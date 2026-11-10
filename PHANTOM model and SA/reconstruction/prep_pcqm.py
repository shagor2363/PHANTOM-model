"""Stage 1 + 2: PCQM records -> anchors, exclusion sectors, band densities, belt infill, reference.

Usage:
    python reconstruction/prep_pcqm.py --out results/reconstruction/geom --seed 42 --check
    python reconstruction/prep_pcqm.py --out results/reconstruction/geom_loto_T3 --exclude-transect 3

Geometry: transect frames from transect_geometry.py (map-digitised layout of Dahdouh-Guebas
et al. 2004), rotated so the parallel transects are vertical; T5 is oblique.  Belt segments
are rotated rectangles (x0, y0, ux, uy, s0, s1, hw); tiles are corridor x band boxes.

Outputs (CSV/JSON in --out):
    records.csv       one row per PCQM quarter (adult present or empty), positions, band, d_cap
    anchors.csv       placed adults (exact distance, random in-quarter angle, belt-constrained)
    sectors.csv       exclusion sectors (px, py, ux, uy, sx, sy, radius) of recorded distances
    band_density.csv  censored PCQM MLE (2026-09-20 B0), Cottam-Curtis and Pollard per band with
                      bootstrap CIs, censoring counts, species shares, lambda_used
    band_density_variants.csv  search-limit sensitivity (R = 20 / 10 m capped; recorded distances uncapped)
    duplicates.csv    flagged duplicate plot blocks of the whole workbook (never deleted)
    core_rects.csv    belt segment rectangles with anchor / infill counts
    tiles.csv         corridor x band tiles with expected counts, eligible area (habitat mask),
                      expected_n_eligible, sampled belt vs extrapolated area (2026-09-20 B1)
    reference.csv     anchors + infill = MTPPR reference per band (scaled marks included)
    allometry.json    height-DBH and crown-DBH fits used for imputation / assembly
    prep_meta.json    counts, settings, frames, segments, corridors, warnings
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recon_config as C  # noqa: E402
import recon_geometry as G  # noqa: E402
import virtual_pcqm as V  # noqa: E402

QUADRAT_RE = re.compile(r'^(\d+)\.(\d+)([a-d])$')


# --- Records -------------------------------------------------------------------------------

def _norm_species(raw):
    s = re.sub(r'\s+', '', str(raw)).lower()
    if s in ('', 'nan', 'none'):
        return None
    if s not in C.SPECIES_MAP:
        raise ValueError(f'unknown species code {raw!r}')
    return C.SPECIES_MAP[s]


JUV97_COLS = list(range(11, 19))     # juvenile 1997 block: 4 x (species, #/m2)
ADULT_MEAS_COLS = [5, 6, 7, 8, 9, 10]  # species, d, GBH, DBH, height, CD
DROP_DUPLICATE_MODES = ('none', 'first', 'second', 'all_suspected')


def _key_value(v):
    """NaN-safe, float-tolerant representation of a workbook cell for block comparison."""
    if v is None:
        return 'nan'
    s = str(v).strip()
    if s.lower() in ('', 'nan', 'none'):
        return 'nan'
    try:
        return f'{float(s):.6g}'
    except ValueError:
        return re.sub(r'\s+', '', s).lower()


def _parse_workbook(xls_path):
    """Every adult-survey row of the workbook (codes 'A' adult, '0' no tree; 'Y' young trees are
    not part of the PCQM adult sample) with its sheet row index and the juvenile 1997 block."""
    df = pd.read_excel(xls_path, sheet_name=C.XLS_SHEET, header=None).iloc[C.XLS_HEADER_ROWS:]
    df.columns = range(df.shape[1])
    parsed = []
    for idx, r in df.iterrows():
        m = QUADRAT_RE.match(str(r[1]).strip())
        if not m:
            continue
        code = str(r[0]).strip()
        if code not in ('A', '0'):
            continue
        parsed.append(dict(row=int(idx), transect=int(m.group(1)), plot=int(m.group(2)), quarter=m.group(3),
                           code=code, raw=r))
    return parsed


def flag_duplicate_plots(parsed):
    """2026-09-20 B0: measurement-level duplicate check of whole PCQM plots (never deletes).

    Two plots are identical when their four quarters carry identical adult measurements
    (species, d, GBH, DBH, height, CD; an empty quarter is 'EMPTY').  Within a transect the
    pair is 'strongly_suspected' when the juvenile 1997 block (cols 11-18) is identical as well,
    otherwise 'suspected'; identical blocks across transects (T5 24/25 vs T6 9/10) are
    'suspected'.  duplicate_order = rank in sheet order within the group of identical plots
    (1 = first).  Plots without any recorded adult are not compared.

    Returns (plot table DataFrame, {(t, plot): dict(status, pair, order)})."""
    blocks = {}
    for r in parsed:
        b = blocks.setdefault((r['transect'], r['plot']), dict(first_row=r['row'], quarters={}))
        b['first_row'] = min(b['first_row'], r['row'])
        qd = b['quarters'].setdefault(r['quarter'], dict(A=None, zero=None))
        if r['code'] == 'A':
            if qd['A'] is not None:
                raise ValueError(f"duplicate adult row {r['transect']}.{r['plot']}{r['quarter']}")
            qd['A'] = r['raw']
        else:
            qd['zero'] = r['raw']
    keys = {}
    for (t, p), b in blocks.items():
        adult, juv, n_present = [], [], 0
        for q in 'abcd':
            qd = b['quarters'].get(q, dict(A=None, zero=None))
            if qd['A'] is not None:
                adult.append((q,) + tuple(_key_value(qd['A'][c]) for c in ADULT_MEAS_COLS))
                n_present += 1
                src = qd['A']
            else:
                adult.append((q, 'EMPTY'))
                src = qd['zero']
            juv.append(tuple(_key_value(src[c]) for c in JUV97_COLS) if src is not None else ('nan',) * len(JUV97_COLS))
        keys[(t, p)] = dict(adult=tuple(adult), juv=tuple(juv), n_present=n_present, first_row=b['first_row'])
    # groups of identical adult blocks (only plots with at least one recorded adult)
    groups = {}
    for tp, k in keys.items():
        if k['n_present'] > 0:
            groups.setdefault(k['adult'], []).append(tp)
    flags = {}
    rows = []
    for _, members in groups.items():
        if len(members) < 2:
            continue
        members = sorted(members, key=lambda tp: keys[tp]['first_row'])
        for i, tp in enumerate(members):
            partners = [m for m in members if m != tp]
            status = 'suspected'
            for m in partners:
                if m[0] == tp[0] and keys[m]['juv'] == keys[tp]['juv']:
                    status = 'strongly_suspected'
            pair = ';'.join(f'{m[0]}.{m[1]}' for m in partners)
            cross = any(m[0] != tp[0] for m in partners)
            flags[tp] = dict(status=status, pair=pair, order=i + 1, cross_transect=cross,
                             group_transects=sorted({m[0] for m in members}))
            rows.append(dict(transect=tp[0], plot=tp[1], plot_id=f'{tp[0]}.{tp[1]}', duplicate_status=status,
                             duplicate_pair=pair, duplicate_order=i + 1, cross_transect=cross,
                             sheet_first_row=keys[tp]['first_row'], n_adults_in_block=keys[tp]['n_present'],
                             juvenile_1997_identical=all(keys[m]['juv'] == keys[tp]['juv'] for m in partners)))
    table = pd.DataFrame(rows, columns=['transect', 'plot', 'plot_id', 'duplicate_status', 'duplicate_pair',
                                        'duplicate_order', 'cross_transect', 'sheet_first_row',
                                        'n_adults_in_block', 'juvenile_1997_identical'])
    return table.sort_values(['transect', 'plot']).reset_index(drop=True), flags


def plots_to_drop(flags, transects, mode):
    """2026-09-20 B0 sensitivity variants.  'first' / 'second': remove the first / second member (sheet
    order) of every STRONGLY suspected pair; 'all_suspected': remove the second member of every flagged
    pair, and, for a cross-transect pair whose first member lies outside the analysed transect set (T5
    24 vs T6 9: T5 is never analysed), the analysed member instead, so that the variant removes both
    suspected T6 plots 9-10 as the plan intends.  Returns the set of (transect, plot) to drop."""
    assert mode in DROP_DUPLICATE_MODES, mode
    drop = set()
    if mode == 'none':
        return drop
    for tp, f in flags.items():
        if tp[0] not in transects:
            continue
        if mode == 'first' and f['status'] == 'strongly_suspected' and f['order'] == 1:
            drop.add(tp)
        elif mode == 'second' and f['status'] == 'strongly_suspected' and f['order'] > 1:
            drop.add(tp)
        elif mode == 'all_suspected':
            if f['order'] > 1:
                drop.add(tp)
            elif f['cross_transect'] and not any(t in transects for t in f['group_transects'] if t != tp[0]):
                drop.add(tp)
    return drop


def load_records(xls_path=C.XLS_PATH, transects=None, drop_duplicates='none', return_meta=False):
    """One row per PCQM quarter of the adult survey (present adults and empty quarters).

    2026-09-20 B0: whole-workbook duplicate flagging (columns duplicate_status / duplicate_pair /
    duplicate_order, nothing deleted), optional removal of one member per flagged pair for the
    sensitivity variants (`drop_duplicates`, see plots_to_drop; the retained record count is reported
    in the meta, never preset), and the list of quarters carrying both a '0' row and an 'A' row
    (treated as adult records).  With return_meta the function returns (rec, meta)."""
    transects = list(transects or C.STAND_TRANSECTS)
    parsed = _parse_workbook(xls_path)
    dup_table, flags = flag_duplicate_plots(parsed)
    drop = plots_to_drop(flags, transects, drop_duplicates)
    rows = []
    quarters = set()
    double_coded = sorted({f"{r['transect']}.{r['plot']}{r['quarter']}" for r in parsed if r['code'] == 'A'} &
                          {f"{r['transect']}.{r['plot']}{r['quarter']}" for r in parsed if r['code'] == '0'},
                          key=lambda s: (int(s.split('.')[0]), int(s.split('.')[1][:-1]), s[-1]))
    n_dropped_adults = 0
    for r0 in parsed:
        t, plot, q, code, r = r0['transect'], r0['plot'], r0['quarter'], r0['code'], r0['raw']
        if t not in transects:
            continue
        if (t, plot) in drop:
            n_dropped_adults += int(code == 'A')
            continue
        quarters.add((t, plot, q))
        if code != 'A':
            continue
        rows.append(dict(transect=t, plot=plot, quarter=q, present=True,
                         species_raw=str(r[5]).strip(), species=_norm_species(r[5]),
                         veg_class=str(r[2]).strip(),
                         d=pd.to_numeric(r[6], errors='coerce'),
                         gbh_cm=pd.to_numeric(r[7], errors='coerce'),
                         dbh_cm=pd.to_numeric(r[8], errors='coerce'),
                         height_m=pd.to_numeric(r[9], errors='coerce'),
                         cd_m=pd.to_numeric(r[10], errors='coerce')))
    adults = pd.DataFrame(rows)
    dup = adults.duplicated(['transect', 'plot', 'quarter'], keep=False)
    if dup.any():
        raise ValueError(f'duplicate adult rows: {adults[dup][["transect", "plot", "quarter"]].values.tolist()}')
    present = set(zip(adults['transect'], adults['plot'], adults['quarter']))
    empty = [dict(transect=t, plot=p, quarter=q, present=False, species_raw='', species=None,
                  veg_class='', d=np.nan, gbh_cm=np.nan, dbh_cm=np.nan, height_m=np.nan, cd_m=np.nan)
             for (t, p, q) in sorted(quarters) if (t, p, q) not in present]
    rec = pd.concat([adults, pd.DataFrame(empty)], ignore_index=True)
    miss = rec['dbh_cm'].isna() & rec['gbh_cm'].notna()
    rec.loc[miss, 'dbh_cm'] = rec.loc[miss, 'gbh_cm'] / math.pi
    rec = rec.sort_values(['transect', 'plot', 'quarter']).reset_index(drop=True)
    rec['record_id'] = [f'{t}.{p}{q}' for t, p, q in zip(rec['transect'], rec['plot'], rec['quarter'])]
    tp = list(zip(rec['transect'], rec['plot']))
    rec['duplicate_status'] = [flags[k]['status'] if k in flags else '' for k in tp]
    rec['duplicate_pair'] = [flags[k]['pair'] if k in flags else '' for k in tp]
    rec['duplicate_order'] = [flags[k]['order'] if k in flags else 0 for k in tp]
    rec['double_coded'] = rec['record_id'].isin(double_coded)
    meta = dict(drop_duplicates=drop_duplicates, transects=transects,
                plots_dropped=sorted(f'{t}.{p}' for t, p in drop), n_plots_dropped=len(drop),
                n_adult_records_dropped=int(n_dropped_adults),
                n_adult_records_retained=int(len(adults)), n_quarters_retained=int(len(rec)),
                duplicate_table=dup_table.to_dict(orient='records'),
                double_coded_quarters=double_coded, n_double_coded_quarters=len(double_coded))
    if return_meta:
        return rec, meta, dup_table
    return rec


def add_positions(rec, frames):
    """PCQM point of each quarter (rotated frame), transect axes, fringe distance, d_cap."""
    rec = rec.copy()
    rec['u'] = [float(G.plot_u(t, k)) for t, k in zip(rec['transect'], rec['plot'])]
    rec['s'] = rec['u']
    px, py = [], []
    for t, u in zip(rec['transect'], rec['u']):
        x, y = G.position(frames[t], u)
        px.append(float(x))
        py.append(float(y))
    rec['px'], rec['py'] = px, py
    for c in ('ux', 'uy', 'nx', 'ny', 'f'):
        rec[c] = rec['transect'].map(lambda t: frames[t][c])
    rec['sx'] = rec['quarter'].map(lambda q: C.QUARTER_SIGNS[q][0])
    rec['sy'] = rec['quarter'].map(lambda q: C.QUARTER_SIGNS[q][1])
    d = rec['d']
    rec['d_cap'] = np.where(rec['present'], np.minimum(d, C.DMAX), C.DMAX)
    rec.loc[rec['present'] & d.isna(), 'd_cap'] = np.nan
    rec['in_stand'] = rec['species'].isin(C.STAND_SPECIES)
    rec['sector_ok'] = rec['present'] & d.notna() & (d <= C.DMAX)
    rec['placeable'] = rec['sector_ok'] & rec['in_stand']
    rec['exclusion_reason'] = ''
    rec.loc[rec['present'] & d.isna(), 'exclusion_reason'] = 'd_missing'
    rec.loc[rec['present'] & (d > C.DMAX), 'exclusion_reason'] = 'd_gt_dmax'
    rec.loc[rec['present'] & ~rec['in_stand'], 'exclusion_reason'] = 'species_not_modelled'
    return rec


# --- Allometry and imputation ---------------------------------------------------------------

def fit_allometry(rec):
    ad = rec[rec['present']]
    out = {'height': {}, 'crown': {}}

    def hfit(sub):
        sub = sub.dropna(subset=['dbh_cm', 'height_m'])
        sub = sub[(sub['dbh_cm'] > 0) & (sub['height_m'] > 0)]
        if len(sub) < 3:
            return None
        X = np.log(sub['dbh_cm'].values)
        Y = np.log(sub['height_m'].values)
        b, a = np.polyfit(X, Y, 1)
        resid = Y - (a + b * X)
        return dict(a=float(a), b=float(b), sd=float(resid.std(ddof=2)) if len(sub) > 2 else 0.0, n=int(len(sub)))

    def cfit(sub):
        sub = sub.dropna(subset=['dbh_cm', 'cd_m'])
        sub = sub[(sub['dbh_cm'] > 0) & (sub['cd_m'] > 0) & (sub['cd_m'] < 40)]
        if len(sub) < 5:
            return None
        c1, c0 = np.polyfit(sub['dbh_cm'].values, sub['cd_m'].values / 2.0, 1)
        return dict(c0=float(c0), c1=float(c1), n=int(len(sub)))

    out['height']['pooled'] = hfit(ad)
    out['crown']['pooled'] = cfit(ad)
    for sp in C.SPECIES_LEVELS:
        sub = ad[ad['species'] == sp]
        out['height'][sp] = hfit(sub) or out['height']['pooled']
        cf = cfit(sub)
        if cf is None or cf['c1'] <= 0 or cf['c0'] < 0:
            cf = out['crown']['pooled']
        out['crown'][sp] = cf
    return out


def impute_marks(rec, allom):
    rec = rec.copy()
    rec['imputed_dbh'] = False
    rec['imputed_h'] = False
    rec['imputed_crown'] = False
    for i in rec.index[rec['present']]:
        sp = rec.at[i, 'species']
        hf = allom['height'].get(sp) or allom['height']['pooled']
        dbh, h = rec.at[i, 'dbh_cm'], rec.at[i, 'height_m']
        if pd.isna(h) and pd.notna(dbh) and hf:
            rec.at[i, 'height_m'] = float(np.exp(hf['a'] + hf['b'] * np.log(dbh)))
            rec.at[i, 'imputed_h'] = True
        elif pd.isna(dbh) and pd.notna(h) and hf and hf['b'] > 0:
            rec.at[i, 'dbh_cm'] = float(np.exp((np.log(h) - hf['a']) / hf['b']))
            rec.at[i, 'imputed_dbh'] = True
        cf = allom['crown'].get(sp) or allom['crown']['pooled']
        cd = rec.at[i, 'cd_m']
        if (pd.isna(cd) or cd <= 0 or cd >= 40) and pd.notna(rec.at[i, 'dbh_cm']) and cf:
            rec.at[i, 'cd_m'] = 2.0 * max(0.3, cf['c0'] + cf['c1'] * rec.at[i, 'dbh_cm'])
            rec.at[i, 'imputed_crown'] = True
    rec['r_stem'] = rec['dbh_cm'] / 200.0
    rec['h_stem'] = rec['height_m']
    rec['r_crown'] = rec['cd_m'] / 2.0
    return rec


# --- Density ----------------------------------------------------------------------------------

def _censored_terms(grp, R, uncap_recorded):
    """2026-09-20 B0: per-quarter terms of the censored PCQM likelihood.

    Uncensored quarter: recorded d <= R (contributes d^2).  Censored: EMPTY quarter (no tree
    recorded; search distance unknown, censored at the assumed R: contributes R^2) and, in the
    capped variant, a recorded d > R (contributes R^2).  With uncap_recorded a recorded d > R
    enters at its recorded distance (uncensored, d^2); empty quarters stay censored at R.  A
    present adult without a recorded d is excluded (missing).  The two censoring cases are kept
    distinct in the counts.  Returns (uncensored flag, contribution, category) arrays."""
    present = grp['present'].values.astype(bool)
    d = grp['d'].values.astype(float)
    has_d = present & np.isfinite(d)
    unc = has_d & (d <= R)
    gt = has_d & (d > R)
    empty = ~present
    missing = present & ~np.isfinite(d)
    contrib = np.zeros(len(grp))
    contrib[unc] = d[unc] ** 2
    contrib[empty] = R ** 2
    if uncap_recorded:
        contrib[gt] = d[gt] ** 2
        unc = unc | gt
    else:
        contrib[gt] = R ** 2
    cat = np.full(len(grp), 'uncensored', dtype=object)
    cat[empty] = 'censored_empty'
    cat[gt] = 'recorded_gt_R_uncapped' if uncap_recorded else 'censored_gt_R'
    cat[missing] = 'missing_d'
    return unc, contrib, cat, missing


def _mle(unc, contrib, missing):
    use = ~missing
    s = float(contrib[use].sum())
    n_unc = int(unc[use].sum())
    return 4.0 * n_unc / (math.pi * s) if s > 0 else np.nan


def band_density(rec, tiles, rng, n_boot=2000, R=None, uncap_recorded=False):
    """Band of a record = band of the tile box containing its PCQM point.

    Per band: Cottam-Curtis and Pollard (checks, on distances capped at R) and, 2026-09-20 B0,
    the censored PCQM maximum-likelihood estimator
        lambda_mle = 4 n_uncensored / (pi * sum_q min(d_q, R)^2)
    with a plot-level bootstrap interval (plots resampled with replacement; rng = the bootstrap
    stream), and the modelled-species share of the recorded adults (Bruguiera removed from the
    total): lambda_mle_species = lambda_mle x share_modelled.  R defaults to C.DMAX (the assumed
    search limit); uncap_recorded implements the 'recorded distances uncapped' variant."""
    R = float(C.DMAX if R is None else R)
    rec = rec.copy()
    rec['tile_id'] = G.tile_of(rec['px'].values, rec['py'].values, tiles)
    rec['band'] = G.band_of_tile(rec['tile_id'].values, tiles)
    if 'd_cap' not in rec.columns or R != C.DMAX:
        d = rec['d']
        rec['d_cap'] = np.where(rec['present'], np.minimum(d, R), R)
        rec.loc[rec['present'] & d.isna(), 'd_cap'] = np.nan
    rows = []
    for k, grp in rec.groupby('band'):
        if k < 0:
            continue
        point_means, point_sumsq = [], []
        plot_terms = []
        for _, g in grp.groupby(['transect', 'plot']):
            dc = g['d_cap'].dropna().values
            if len(dc):
                point_means.append(dc.mean())
                point_sumsq.append((dc ** 2).sum())
            u, c, _, m = _censored_terms(g, R, uncap_recorded)
            plot_terms.append((u, c, m))
        point_means = np.array(point_means)
        n_pts = len(point_means)
        dcap_all = grp['d_cap'].dropna().values
        lam = 1.0 / np.mean(dcap_all) ** 2
        pollard = 4 * (4 * n_pts - 1) / (math.pi * np.sum(point_sumsq)) if n_pts else np.nan
        boots = [1.0 / np.mean(point_means[rng.integers(0, n_pts, n_pts)]) ** 2 for _ in range(n_boot)]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        # censored MLE over the band's quarters, plot-level bootstrap
        unc, contrib, cat, missing = _censored_terms(grp, R, uncap_recorded)
        lam_mle = _mle(unc, contrib, missing)
        n_plots = len(plot_terms)
        boots_mle = []
        for _ in range(n_boot):
            pick = rng.integers(0, n_plots, n_plots)
            bu = np.concatenate([plot_terms[i][0] for i in pick])
            bc = np.concatenate([plot_terms[i][1] for i in pick])
            bm = np.concatenate([plot_terms[i][2] for i in pick])
            boots_mle.append(_mle(bu, bc, bm))
        mlo, mhi = np.nanpercentile(boots_mle, [2.5, 97.5])
        pres = grp[grp['present'] & grp['in_stand']]
        share_modelled = float(len(pres) / grp['present'].sum()) if grp['present'].sum() else np.nan
        row = dict(band=int(k), n_points=int(n_pts), n_quarters=int(len(grp)), n_present=int(grp['present'].sum()),
                   n_present_modelled=int(len(pres)), share_modelled=share_modelled,
                   n_empty=int((~grp['present']).sum()), n_capped=int((grp['d'] > R).sum()),
                   n_d_missing=int((grp['present'] & grp['d'].isna()).sum()), mean_d_cap=float(np.mean(dcap_all)),
                   lambda_cc=float(lam), lambda_lo=float(lo), lambda_hi=float(hi), lambda_pollard=float(pollard),
                   R_search=R, search_limit_assumed=bool(C.SEARCH_LIMIT_ASSUMED),
                   recorded_variant='uncapped_recorded' if uncap_recorded else 'capped_at_R',
                   n_uncensored=int(unc[~missing].sum()), n_censored_empty=int((cat == 'censored_empty').sum()),
                   n_censored_gt_R=int((cat == 'censored_gt_R').sum()),
                   n_recorded_gt_R_uncapped=int((cat == 'recorded_gt_R_uncapped').sum()),
                   sum_min_d2=float(contrib[~missing].sum()),
                   lambda_mle=float(lam_mle), lambda_mle_lo=float(mlo), lambda_mle_hi=float(mhi),
                   lambda_mle_species=float(lam_mle * share_modelled) if np.isfinite(share_modelled) else np.nan,
                   transects=' '.join(str(t) for t in sorted(grp['transect'].unique())))
        for sp in C.SPECIES_LEVELS:
            row[f'n_{sp}'] = int((pres['species'] == sp).sum())
            row[f'p_{sp}'] = float((pres['species'] == sp).mean()) if len(pres) else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values('band').reset_index(drop=True), rec


def band_density_variants(rec, tiles, rng, n_boot=2000):
    """2026-09-20 B0: search-limit sensitivity table: R in C.SEARCH_LIMIT_VARIANTS (capped) and the
    'recorded distances uncapped' variant at the main R (empty quarters censored at R)."""
    parts = []
    for R in C.SEARCH_LIMIT_VARIANTS:
        d, _ = band_density(rec, tiles, rng, n_boot=n_boot, R=R, uncap_recorded=False)
        d.insert(0, 'variant', f'R{R:g}_capped')
        parts.append(d)
    d, _ = band_density(rec, tiles, rng, n_boot=n_boot, R=C.DMAX, uncap_recorded=True)
    d.insert(0, 'variant', f'R{C.DMAX:g}_recorded_uncapped')
    parts.append(d)
    return pd.concat(parts, ignore_index=True)


# --- Sectors and anchors ------------------------------------------------------------------------

def build_sectors(rec):
    """A recorded nearest tree at distance d means no tree closer than d in that quarter."""
    sec = rec[rec['sector_ok']].copy()
    sec['radius'] = sec['d'] - 1e-3
    sec['kind'] = np.where(sec['in_stand'], 'recorded', 'recorded_other_species')
    return sec[['record_id', 'transect', 'plot', 'quarter', 'px', 'py', 'ux', 'uy', 'sx', 'sy', 'radius', 'kind']].reset_index(drop=True)


def _sector_hits(xp, yp, sectors, exclude):
    keep = ~np.asarray(exclude, dtype=bool)
    s = sectors[keep]
    dx = xp - s['px'].values
    dy = yp - s['py'].values
    along = dx * s['ux'].values + dy * s['uy'].values
    perp = dx * s['uy'].values - dy * s['ux'].values
    hit = (s['sx'].values * perp > 0) & (s['sy'].values * along > 0) & (dx * dx + dy * dy < s['radius'].values ** 2)
    return int(hit.sum())


def place_anchors(rec, sectors, rng, max_tries=200):
    """Exact-distance, random-angle placement inside the quarter (transect-local frame),
    belt-constrained, avoiding all other exclusion sectors; unresolvable conflicts flagged."""
    rows = []
    n_constrained = 0
    n_conflict = 0
    sec_ids = sectors['record_id'].values
    for r in rec[rec['placeable']].itertuples(index=False):
        d = float(r.d)
        phi_min = math.acos(min(1.0, (C.BELT_HALF - 0.01) / d)) if d > 0 else 0.0
        if phi_min > 0:
            n_constrained += 1
        own = sec_ids == r.record_id
        ok = False
        best = None
        for _ in range(max_tries):
            phi = rng.uniform(phi_min, math.pi / 2)
            dperp = r.sx * d * math.cos(phi)
            dalong = r.sy * d * math.sin(phi)
            xp = r.px + dalong * r.ux + dperp * r.nx
            yp = r.py + dalong * r.uy + dperp * r.ny
            hits = _sector_hits(xp, yp, sectors, own)
            if best is None or hits < best[0]:
                best = (hits, phi, dalong, dperp, xp, yp)
            if hits == 0:
                ok = True
                break
        if not ok:
            n_conflict += 1
            _, phi, dalong, dperp, xp, yp = best
        rows.append(dict(record_id=r.record_id, transect=r.transect, plot=r.plot, quarter=r.quarter,
                         xp=xp, yp=yp, dalong=dalong, dperp=dperp, phi=phi, angle_constrained=phi_min > 0,
                         sector_conflict=not ok, species=r.species, d=d, dbh_cm=r.dbh_cm,
                         height_m=r.height_m, r_stem=r.r_stem, h_stem=r.h_stem, r_crown=r.r_crown,
                         imputed_dbh=r.imputed_dbh, imputed_h=r.imputed_h, veg_class=r.veg_class))
    anchors = pd.DataFrame(rows)
    anchors['anchor_id'] = np.arange(len(anchors))
    return anchors, dict(n_angle_constrained=n_constrained, n_sector_conflicts=n_conflict)


# --- Infill ----------------------------------------------------------------------------------------

def infill_rects(core_rects, anchors, sectors, dens, frames, rng, max_batches=400):
    """Dart-throwing infill of each belt piece to lambda_k * area, honouring sectors of all
    transects, a hard core against all trees, and pieces already filled (parallel first)."""
    lam = dict(zip(dens['band'], dens['lambda_cc']))
    rows, stats = [], []
    order = list(core_rects[core_rects['clipped'] == 0].index) + list(core_rects[core_rects['clipped'] == 1].index)
    existing = np.column_stack([anchors['xp'].values, anchors['yp'].values]) if len(anchors) else np.zeros((0, 2))
    filled = []
    for idx in order:
        rct = core_rects.loc[idx]
        a_in = anchors[anchors['rect_id'] == rct.rect_id]
        others = core_rects[core_rects['rect_id'].isin(filled)]
        # effective area: piece minus parts inside pieces filled earlier
        mx, my, frac = G.sample_in_rect(rct, rng, 20000)
        overlap = float(G.in_rects(mx, my, others).mean()) if (len(others) and len(mx)) else 0.0
        area_eff = rct.area * (1 - overlap)
        target = max(0, int(round(lam.get(int(rct.band), 0.0) * area_eff)) - len(a_in))
        new = []
        placed = 0
        batches = 0
        while placed < target and batches < max_batches:
            batches += 1
            cx, cy, _ = G.sample_in_rect(rct, rng, max(200, 8 * (target - placed)))
            keep = ~G.in_sector(cx, cy, sectors)
            if len(others):
                keep &= ~G.in_rects(cx, cy, others)
            cx, cy = cx[keep], cy[keep]
            pool = np.vstack([existing, np.array(new).reshape(-1, 2)]) if new else existing
            tree = cKDTree(pool) if len(pool) else None
            for x_, y_ in zip(cx, cy):
                if placed >= target:
                    break
                if tree is not None:
                    dd, _ = tree.query([x_, y_], k=1)
                    if dd < C.HARDCORE:
                        continue
                if new:
                    dn = np.hypot(np.array(new)[:, 0] - x_, np.array(new)[:, 1] - y_).min()
                    if dn < C.HARDCORE:
                        continue
                new.append((x_, y_))
                placed += 1
        if new:
            existing = np.vstack([existing, np.array(new)])
        src = anchors[anchors['band'] == rct.band]
        if len(src) == 0:
            src = anchors
        if placed:
            pick = rng.integers(0, len(src), placed)
            for (x_, y_), j in zip(new, pick):
                sr = src.iloc[j]
                rows.append(dict(transect=int(rct.transect), band=int(rct.band), rect_id=int(rct.rect_id), xp=x_, yp=y_,
                                 species=sr['species'], r_stem=sr['r_stem'], h_stem=sr['h_stem'],
                                 r_crown=sr['r_crown'], donor_anchor_id=int(sr['anchor_id'])))
        stats.append(dict(rect_id=int(rct.rect_id), n_anchors=len(a_in), area_eff=area_eff, overlap_frac=overlap,
                          n_infill_target=target, n_infill_placed=placed, batches=batches))
        filled.append(rct.rect_id)
    cols = ['transect', 'band', 'rect_id', 'xp', 'yp', 'species', 'r_stem', 'h_stem', 'r_crown', 'donor_anchor_id']
    return pd.DataFrame(rows, columns=cols), pd.DataFrame(stats)


# --- Pipeline ----------------------------------------------------------------------------------------

def run(out_dir, seed=C.ANCHOR_SEED, exclude_transect=None, check=False, verbose=True, density_mode=None,
        band_breaks=None, drop_duplicates='none', corridor_transects=None, n_boot=2000):
    """2026-09-20 B0/B1: geometry (corridors, tile and belt segments, tiles) is derived from the
    RETAINED transects only (corridor_transects defaults to them); band_breaks (default
    C.BAND_BREAKS) are persisted in prep_meta['settings']; three RNG streams: rng_anchor (seeded
    by `seed` alone, so the anchor realisation depends only on the seed), rng_boot (bootstrap
    intervals) and rng_geom (Monte-Carlo area integrals); duplicate records are flagged, and
    optionally one member per flagged pair removed (drop_duplicates, sensitivity variants only)."""
    density_mode = density_mode or C.DENSITY_MODE
    assert density_mode in C.DENSITY_MODES, density_mode
    os.makedirs(out_dir, exist_ok=True)
    rng_anchor = np.random.default_rng(seed)
    rng_boot = np.random.default_rng([seed, 1])
    rng_geom = np.random.default_rng([seed, 2])
    band_breaks = [float(b) for b in (C.BAND_BREAKS if band_breaks is None else band_breaks)]
    warnings = []
    # 2026-09-20: exclude_transect may be one id or an iterable of ids (multi-transect hold-out)
    excl = set() if exclude_transect is None else set(np.atleast_1d(exclude_transect).astype(int).tolist())
    transects = [t for t in C.STAND_TRANSECTS if t not in excl]
    corridor_transects = list(corridor_transects or transects)
    assert set(corridor_transects) <= set(transects), 'corridors must come from retained transects'

    # 2026-09-20 B1: frames of the RETAINED transects only -> corridors, segments, tiles
    frames_all = G.transect_frames(corridor_transects)
    frames = {t: frames_all[t] for t in transects}
    corridors = G.corridor_bounds(frames_all)
    rec, rec_meta, dup_table = load_records(transects=transects, drop_duplicates=drop_duplicates, return_meta=True)
    rec = add_positions(rec, frames_all)
    allom = fit_allometry(rec)
    rec = impute_marks(rec, allom)

    # belt segments (surveyed part of each belt) and tile segments (whole corridor, from u = 0)
    belt_segs = {t: G.band_segments(frames_all[t]['tape_length'], u_start=frames_all[t]['tape_offset'], breaks=band_breaks)
                 for t in frames_all}
    tile_segs = {t: G.band_segments(frames_all[t]['tape_length'], breaks=band_breaks) for t in frames_all}
    tiles = G.build_tiles(frames_all, tile_segs, corridors)
    dens, rec = band_density(rec, tiles, rng_boot, n_boot=n_boot)
    dens_var = band_density_variants(rec, tiles, rng_boot, n_boot=n_boot)
    n_unbanded = int((rec['band'] < 0).sum())
    if n_unbanded:
        warnings.append(f'{n_unbanded} PCQM quarters outside all tile boxes (excluded from band densities)')

    sectors = build_sectors(rec)
    anchors, place_stats = place_anchors(rec, sectors, rng_anchor)

    # extend belt ends for overshooting anchors (in tape metres)
    for t in transects:
        a_t = anchors[anchors['transect'] == t]
        if len(a_t):
            along, _ = G.local_coords(frames_all[t], a_t['xp'].values, a_t['yp'].values)
            u = along / frames_all[t]['f']
            lo = float(np.floor(u.min() * 2) / 2 - 0.5)
            hi = float(np.ceil(u.max() * 2) / 2 + 0.5)
            u_min = min(frames_all[t]['tape_offset'], lo)
            u_max = max(frames_all[t]['tape_length'], hi)
            if u_min < frames_all[t]['tape_offset'] or u_max > frames_all[t]['tape_length']:
                warnings.append(f'T{t}: belt extended to u in [{u_min}, {u_max}] for overshooting anchors')
                belt_segs[t] = G.band_segments(frames_all[t]['tape_length'], u_start=frames_all[t]['tape_offset'],
                                               u_min=u_min, u_max=u_max, breaks=band_breaks)
    core_rects = G.build_core_rects(frames, {t: belt_segs[t] for t in frames}, tiles, rng_geom)

    # anchors -> first containing piece (parallel belts first), band = band of that piece
    ri = G.rect_index(anchors['xp'].values, anchors['yp'].values, core_rects)
    anchors['rect_id'] = np.where(ri >= 0, core_rects['rect_id'].values[np.clip(ri, 0, None)], -1)
    anchors['band'] = np.where(ri >= 0, core_rects['band'].values[np.clip(ri, 0, None)], -1)
    anchors['tile_id'] = G.tile_of(anchors['xp'].values, anchors['yp'].values, tiles)
    n_out = int((anchors['rect_id'] < 0).sum())
    if n_out:
        warnings.append(f'{n_out} anchors outside every belt piece (kept as orphans)')

    dens['lambda_used'] = dens['lambda_cc'] * dens['share_modelled']   # Cottam-Curtis species density ('pcqm' mode)
    if density_mode in ('recorded', 'pcqm_mle'):
        # Reference belts = the recorded trees only, in belts as wide as the trees they hold
        # (Pambala recipe), no infill.  'recorded': the belt area is also the density denominator
        # (recorded trees per belt area).  'pcqm_mle' (2026-09-20 B0): the belts remain the
        # REFERENCE geometry, but lambda_used is the censored PCQM MLE x modelled-species share;
        # the recorded-per-belt density is kept alongside as `lambda_recorded_belt`.
        for i, rct in core_rects.iterrows():
            a = anchors[anchors['rect_id'] == rct.rect_id]
            if len(a):
                _, perp = G.local_coords(frames_all[int(rct.transect)], a['xp'].values, a['yp'].values)
                hw = float(np.clip(np.abs(perp).max() + 0.25, C.RECORDED_MIN_HALF, C.BELT_HALF))
            else:
                hw = C.RECORDED_MIN_HALF
            core_rects.at[i, 'hw'] = hw
            core_rects.at[i, 'area'] = (rct.s1 - rct.s0) * 2 * hw
        n_by_band = anchors[anchors['rect_id'] >= 0].groupby('band').size()
        a_by_band = core_rects.groupby('band')['area'].sum()
        dens['lambda_recorded_belt'] = [float(n_by_band.get(b, 0) / a_by_band.get(b, np.nan)) for b in dens['band']]
        if density_mode == 'recorded':
            dens['lambda_used'] = dens['lambda_recorded_belt']
        else:
            dens['lambda_used'] = dens['lambda_mle_species']
        infill = pd.DataFrame(columns=['transect', 'band', 'rect_id', 'xp', 'yp', 'species', 'r_stem', 'h_stem', 'r_crown', 'donor_anchor_id'])
        infill_stats = pd.DataFrame([dict(rect_id=int(r.rect_id), n_anchors=int((anchors['rect_id'] == r.rect_id).sum()),
                                          area_eff=float(r.area), overlap_frac=0.0, n_infill_target=0, n_infill_placed=0, batches=0)
                                     for r in core_rects.itertuples(index=False)])
    else:
        infill, infill_stats = infill_rects(core_rects, anchors, sectors, dens, frames, rng_geom)
    core_rects = core_rects.merge(infill_stats, on='rect_id', how='left').fillna(
        {'n_anchors': 0, 'area_eff': 0.0, 'overlap_frac': 0.0, 'n_infill_target': 0, 'n_infill_placed': 0, 'batches': 0})
    for r in core_rects[core_rects['n_infill_placed'] < core_rects['n_infill_target']].itertuples(index=False):
        warnings.append(f'rect {r.rect_id}: infill shortfall {int(r.n_infill_placed)}/{int(r.n_infill_target)}')

    a_ref = anchors[['anchor_id', 'band', 'transect', 'rect_id', 'xp', 'yp', 'species', 'r_stem', 'h_stem', 'r_crown']].copy()
    a_ref['source'] = 'anchor'
    i_ref = infill[['band', 'transect', 'rect_id', 'xp', 'yp', 'species', 'r_stem', 'h_stem', 'r_crown']].copy()
    i_ref['source'] = 'infill'
    i_ref['anchor_id'] = -1
    ref = pd.concat([a_ref, i_ref], ignore_index=True)
    ref = ref.reset_index(drop=True)
    ref['ref_id'] = np.arange(len(ref))
    ref['tile_id'] = G.tile_of(ref['xp'].values, ref['yp'].values, tiles)
    ref['r_stem_s'] = ref['r_stem'] / C.R_STEM_SCALE
    ref['h_stem_s'] = ref['h_stem'] / C.H_STEM_SCALE
    x, y = G.to_model(ref['xp'].values, ref['yp'].values)
    ref['x'], ref['y'] = x, y

    lam = dict(zip(dens['band'], dens['lambda_used']))
    tiles['lambda'] = tiles['band'].map(lam).fillna(0.0)
    tiles['expected_n'] = tiles['lambda'] * tiles['area']
    tiles['n_ref'] = tiles['band'].map(ref.groupby('band').size()).fillna(0).astype(int)
    tiles['core_area_union'] = tiles['band'].map(
        {k: G.union_area(core_rects[core_rects['band'] == k], rng_geom) for k in core_rects['band'].unique()}).fillna(0.0)
    tiles['has_core'] = [int(any(G.rect_touches_box(r, tl) for r in core_rects.itertuples(index=False))) for tl in tiles.itertuples(index=False)]
    npts = rec[rec['quarter'] == 'a'].groupby(['transect', 'band']).size()
    tiles['n_pcqm_points'] = [int(npts.get((t, k), 0)) for t, k in zip(tiles['transect'], tiles['band'])]
    tiles['n_anchors'] = tiles['tile_id'].map(ref[ref['source'] == 'anchor'].groupby('tile_id').size()).fillna(0).astype(int)
    tiles['n_infill'] = tiles['tile_id'].map(ref[ref['source'] == 'infill'].groupby('tile_id').size()).fillna(0).astype(int)
    tiles['in_domain_frac'] = [_rect_domain_fraction(r, rng_geom) for r in tiles.itertuples(index=False)]
    # 2026-09-20 B1: eligible area (tile minus the treeless cells of the habitat mask), expected
    # count in it, and the sampled (own reference belt) vs extrapolated part of every tile
    elig = G.tile_eligible_areas(tiles)
    tiles = tiles.merge(elig, on='tile_id', how='left')
    tiles['eligible_area'] = np.minimum(tiles['eligible_area'], tiles['area'] * tiles['in_domain_frac'])
    tiles['expected_n_eligible'] = tiles['lambda'] * tiles['eligible_area']
    own = {(int(r.transect), int(r.band)): float(r.area) for r in core_rects.itertuples(index=False) if not getattr(r, 'clipped', 0)}
    tiles['sampled_belt_area'] = [own.get((int(t), int(k)), 0.0) for t, k in zip(tiles['transect'], tiles['band'])]
    tiles['extrapolated_area'] = np.maximum(tiles['eligible_area'] - tiles['sampled_belt_area'], 0.0)

    rec.to_csv(os.path.join(out_dir, 'records.csv'), index=False)
    anchors.to_csv(os.path.join(out_dir, 'anchors.csv'), index=False)
    sectors.to_csv(os.path.join(out_dir, 'sectors.csv'), index=False)
    dens.to_csv(os.path.join(out_dir, 'band_density.csv'), index=False)
    dens_var.to_csv(os.path.join(out_dir, 'band_density_variants.csv'), index=False)
    dup_table.to_csv(os.path.join(out_dir, 'duplicates.csv'), index=False)
    core_rects.to_csv(os.path.join(out_dir, 'core_rects.csv'), index=False)
    tiles.to_csv(os.path.join(out_dir, 'tiles.csv'), index=False)
    ref.to_csv(os.path.join(out_dir, 'reference.csv'), index=False)
    infill.to_csv(os.path.join(out_dir, 'infill.csv'), index=False)
    with open(os.path.join(out_dir, 'allometry.json'), 'w') as f:
        json.dump(allom, f, indent=2)
    x0, x1, y0, y1 = G.domain_bbox_rot()
    dup_meta = {k: v for k, v in rec_meta.items() if k != 'duplicate_table'}
    meta = dict(seed=seed, exclude_transect=(None if exclude_transect is None else sorted(excl)), transects=transects,
                corridor_transects=corridor_transects, density_mode=density_mode, drop_duplicates=drop_duplicates,
                n_quarters=int(len(rec)), n_adults=int(rec['present'].sum()), n_anchors=int(len(anchors)),
                n_infill=int(len(infill)), n_reference=int(len(ref)), n_orphan_reference=int((ref['tile_id'] < 0).sum()),
                n_unplaceable=int((rec['present'] & ~rec['placeable']).sum()), n_empty_quarters=int((~rec['present']).sum()),
                n_other_species=int((rec['present'] & ~rec['in_stand']).sum()), stand_species=C.STAND_SPECIES,
                n_adults_by_transect={str(t): int(n) for t, n in rec[rec['present']].groupby('transect').size().items()},
                n_anchors_by_transect={str(t): int(n) for t, n in anchors.groupby('transect').size().items()},
                duplicates=dup_meta,
                settings=dict(DMAX=C.DMAX, SEARCH_LIMIT_ASSUMED=C.SEARCH_LIMIT_ASSUMED, BELT_HALF=C.BELT_HALF,
                              HARDCORE=C.HARDCORE, BAND_BREAKS=C.BAND_BREAKS, band_breaks=band_breaks,
                              transects=transects, corridors={str(t): list(corridors[t]) for t in corridors},
                              MIN_REMNANT=C.MIN_REMNANT, R_STEM_SCALE=C.R_STEM_SCALE, H_STEM_SCALE=C.H_STEM_SCALE,
                              density_mode=density_mode, COUNT_CONTROL=C.COUNT_CONTROL, drop_duplicates=drop_duplicates,
                              rng_streams=dict(rng_anchor=f'default_rng({seed})', rng_boot=f'default_rng([{seed}, 1])',
                                               rng_geom=f'default_rng([{seed}, 2])'),
                              substrate_npz=os.path.relpath(C.SUBSTRATE_NPZ, C.ROOT)),
                belt_segments={str(t): [list(s) for s in belt_segs[t]] for t in belt_segs},
                tile_segments={str(t): [list(s) for s in tile_segs[t]] for t in tile_segs},
                frames={str(t): frames_all[t] for t in frames_all},
                corridors={str(t): list(corridors[t]) for t in corridors},
                domain_bbox_rot=[x0, x1, y0, y1], warnings=warnings, **place_stats)
    with open(os.path.join(out_dir, 'prep_meta.json'), 'w') as f:
        json.dump(meta, f, indent=2)

    if verbose:
        print(f'transects {transects} (corridors from {corridor_transects}), band breaks {band_breaks}, '
              f'drop_duplicates={drop_duplicates} (plots dropped {rec_meta["plots_dropped"]})')
        print(f'quarters {len(rec)}, adults {int(rec["present"].sum())}, anchors {len(anchors)} '
              f'(angle-constrained {place_stats["n_angle_constrained"]}, conflicts {place_stats["n_sector_conflicts"]}), '
              f'infill {len(infill)}, reference {len(ref)} (orphans outside tiles {int((ref["tile_id"] < 0).sum())})')
        print(f'density mode: {density_mode} (R = {C.DMAX} m, search limit assumed = {C.SEARCH_LIMIT_ASSUMED})')
        cols = ['band', 'n_points', 'n_quarters', 'n_present', 'n_present_modelled', 'n_uncensored', 'n_censored_empty',
                'n_censored_gt_R', 'lambda_mle', 'lambda_mle_lo', 'lambda_mle_hi', 'lambda_cc', 'lambda_pollard',
                'share_modelled', 'lambda_used', 'transects']
        print(dens[cols].round(4).to_string(index=False))
        print(core_rects[['rect_id', 'tile_id', 'transect', 'band', 'hw', 'area', 'area_eff', 'n_anchors', 'n_infill_target', 'n_infill_placed']].round(1).to_string(index=False))
        print(tiles[['tile_id', 'transect', 'band', 'area', 'eligible_area', 'lambda', 'expected_n', 'expected_n_eligible',
                     'sampled_belt_area', 'extrapolated_area', 'n_ref']].round(2).to_string(index=False))
        for w in warnings:
            print('WARNING:', w)
    if check:
        run_checks(rec, anchors, sectors, core_rects, tiles, ref, frames_all, belt_segs, dens, density_mode)
    return dict(rec=rec, anchors=anchors, sectors=sectors, dens=dens, dens_variants=dens_var, duplicates=dup_table,
                core_rects=core_rects, tiles=tiles, ref=ref, infill=infill, frames=frames_all, belt_segments=belt_segs,
                tile_segments=tile_segs, corridors=corridors, meta=meta)


def _rect_domain_fraction(r, rng, n=4000):
    xp = rng.uniform(r.xp_min, r.xp_max, n)
    yp = rng.uniform(r.yp_min, r.yp_max, n)
    x, y = G.to_model(xp, yp)
    return float(G.in_domain(x, y).mean())


def run_checks(rec, anchors, sectors, core_rects, tiles, ref, frames, segments, dens, density_mode=None):
    density_mode = density_mode or C.DENSITY_MODE
    print('--- checks ---')
    for t in sorted(set(rec['transect'])):
        n_q = int((rec['transect'] == t).sum())
        assert n_q == 4 * frames[t]['n_plots'], f'T{t}: {n_q} quarters'
    for t in sorted(set(anchors['transect'])):
        a_t = anchors[anchors['transect'] == t]
        _, perp = G.local_coords(frames[t], a_t['xp'].values, a_t['yp'].values)
        assert (np.abs(perp) < C.BELT_HALF).all(), f'T{t}: anchor outside belt width'
    inf = ref[ref['source'] == 'infill']
    assert not G.in_sector(inf['xp'].values, inf['yp'].values, sectors).any(), 'infill inside a sector'
    pts = ref[['xp', 'yp']].values
    tree = cKDTree(pts)
    pairs = np.array(list(tree.query_pairs(C.HARDCORE - 1e-9)))
    if len(pairs):
        src = ref['source'].values
        bad = [(i, j) for i, j in pairs if src[i] == 'infill' or src[j] == 'infill']
        assert not bad, f'{len(bad)} infill pairs closer than hard core'
        print(f'note: {len(pairs)} anchor-anchor pairs closer than {C.HARDCORE} m (measured data, kept)')
    for i in range(len(tiles)):
        for j in range(i + 1, len(tiles)):
            a, b = tiles.iloc[i], tiles.iloc[j]
            ov = (a.xp_min < b.xp_max and b.xp_min < a.xp_max and a.yp_min < b.yp_max and b.yp_min < a.yp_max)
            assert not ov, f'tiles {a.tile_id} and {b.tile_id} overlap'
    # every reference point sits in a tile of its own band (or in no tile)
    tb = G.band_of_tile(ref['tile_id'].values, tiles)
    bad = ((tb >= 0) & (tb != ref['band'].values)).sum()
    assert bad == 0, f'{bad} reference points banded differently from their tile'
    for k, g in ref[ref['band'] >= 0].groupby('band'):
        area = core_rects.loc[core_rects['band'] == k, 'area_eff'].sum()
        lam = float(dens.loc[dens['band'] == k, 'lambda_used'].iloc[0])
        rel = len(g) / area / lam - 1
        print(f'band {k}: reference density {len(g) / area:.4f} vs lambda {lam:.4f} (rel {rel:+.3f})')
        if density_mode == 'pcqm_mle':
            # 2026-09-20 B0: the reference belts hold the recorded (nearest-tree) sample only; their
            # density is BELOW the distance-based lambda by design (METHODS limitation 3), so the
            # equality check of the filled-belt modes does not apply; the ratio is reported.
            continue
        assert abs(rel) < 0.06, f'band {k} density off by {rel:.3f}'
    n_ok = n_explained = n_bad = 0
    conflict_ids = set(anchors.loc[anchors['sector_conflict'], 'anchor_id'])
    stand = ref.reset_index(drop=True)
    for t in sorted(set(rec['transect'])):
        vs = V.survey_transect(stand, t, frames)
        real = rec[rec['transect'] == t].set_index('record_id')
        for v in vs.itertuples(index=False):
            r = real.loc[f'{t}.{v.plot}{v.quarter}']
            if not r['placeable']:
                continue
            if v.present and abs(v.d - r['d']) < 2e-3:   # sector margin is 1e-3
                n_ok += 1
            else:
                found = int(stand.loc[v.idx, 'anchor_id']) if v.idx >= 0 else -1
                if found in conflict_ids:
                    n_explained += 1
                else:
                    n_bad += 1
                    print(f'  mismatch {t}.{v.plot}{v.quarter}: real d={r["d"]}, virtual present={v.present} d={v.d} (ref idx {v.idx})')
    print(f'virtual survey (all reference trees): {n_ok} anchors reproduced exactly, {n_explained} mismatches '
          f'explained by flagged conflicts, {n_bad} unexplained')
    assert n_bad == 0
    print('all checks passed')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=C.GEOM_DIR)
    ap.add_argument('--seed', type=int, default=C.ANCHOR_SEED)
    # 2026-09-20 B0/B1: several transects may be excluded; band breaks and duplicate handling explicit
    ap.add_argument('--exclude-transect', type=int, nargs='+', default=None)
    ap.add_argument('--band-breaks', type=float, nargs='+', default=None,
                    help=f'fringe-distance breaks (m), last band open; default {C.BAND_BREAKS}')
    ap.add_argument('--drop-duplicates', choices=list(DROP_DUPLICATE_MODES), default='none',
                    help='sensitivity variants only: remove one member of each flagged duplicate plot pair')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--density-mode', choices=list(C.DENSITY_MODES), default=None)
    args = ap.parse_args()
    run(args.out, seed=args.seed, exclude_transect=args.exclude_transect, check=args.check, density_mode=args.density_mode,
        band_breaks=args.band_breaks, drop_duplicates=args.drop_duplicates)


if __name__ == '__main__':
    main()
