"""Per-run summaries and pattern scores (written for every run by transport_treatments.run_grid).

One summary condenses a 365-day run (about 200k propagules, a 40 MB cached CSV) into a
small JSON so that calibration objectives and figures never reload the caches:

  strata      all / boundary / internal propagules
  stages      parents (the stand), first grounding, terminal grounded (established + dead),
              established
  per stage   species zone counts and shares (position zones), 2 cm elevation histograms,
              quantiles, pairwise separation between species (W1 in m, JSD base 2)
  scores      against the 1999 and 1997 position-zoned juvenile surveys: 9-cell Spearman
              rho, peak-zone matches, within-zone species rank, per-zone composition JSD
              (validate_1999_patterns), per-species zone-share JSD ("zonation error")
  signatures  S1 Ceriops the most landward recruiter; S2 Rhizophora the most seaward with its
              peak in zone B; S3 Avicennia shows the largest landward shift of recruits
              relative to adults (S1-S3 are ranked among the CALIBRATED species only);
              S4 (reported, not gated; 2026-09-20): Sonneratia's share of the recruits inside
              the transects and the 10th/90th percentiles of its terminal elevation (observed
              benchmark: < 1 % of the juvenile density, quadrat elevations 1.37-1.66 m)
  guardrails  Z3 mortality gradient, Z6 predation dominant in zone C, landward-clamp share

The observation operator is the one of the validation pipeline: established seedlings
inside the survey transect polygons (T1-T5), zone by position, divided by the ASSUMED
virtual-sampling footprint per zone (n_quadrats x transect_geometry.MODEL_QUADRAT_FOOTPRINT_M2,
25 m2 = the paper's 5 x 5 m quadrat; data audit 2026-09-20, docs/juvenile_data_audit.md);
zone shares are RELATIVE DENSITIES ACROSS ZONES (mean density per zone / sum over zones),
exactly as for the field data, and do not depend on the footprint. Absolute sampled
densities are conditional on the footprint assumption.

Observed data (2026-09-20): the survey records per-quadrat DENSITIES (#/m2) of unknown and
possibly unequal quadrat area; no verified seedling counts exist, so the chi-square goodness
of fit of observed counts (`chi2_tests`) is DROPPED for the observed data: it returns
p = None with a note for every species, and comparisons rest on the density-based zone
shares (JSD 'zonation error', bootstrap TVD in equivalence_test.py / virtual_survey.py).

Usage:
    python cascade_metrics.py --cache results/knockout_geom2/knockout_baseline_seed42.csv \
        --stand results/reconstruction/demo_reconstructed.csv
    python cascade_metrics.py --cache-dir results/knockout_geom2 --stand ... [--force]
"""
import argparse
import csv
import glob
import json
import math
import os
import sys
from collections import Counter
from itertools import combinations

import numpy as np
from scipy.stats import wasserstein_distance

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import tidal_dispersal_model as tdm  # noqa: E402
import transect_geometry as tg  # noqa: E402
import validate_1999_patterns as v99  # noqa: E402

# 2026-09-20: 4 species (Sonneratia); the species list is the model's, not re-declared here.
# CALIBRATED_SPECIES = the three species whose zonation was the calibration target; the
# signatures, the calibration objective and every '9-cell' score are evaluated on them so
# that the numbers stay comparable with the 3-species runs. '_all' variants cover SPECIES.
from tidal_dispersal_model import SPECIES as _MODEL_SPECIES, CALIBRATED_SPECIES as _MODEL_CALIBRATED  # noqa: E402
SPECIES = list(_MODEL_SPECIES)
CALIBRATED_SPECIES = list(_MODEL_CALIBRATED)
ZONES = ['A', 'B', 'C']
STRATA = {'all': None, 'boundary': 'boundary', 'internal': 'internal'}
PAIRS = list(combinations(SPECIES, 2))
MIN_NONZERO_QUADRATS = 5   # below this the observed composition of a species is not testable
Z_EDGES = np.round(np.linspace(1.0, 2.2, 61), 4)   # 2 cm elevation bins
QUANTILES = [0.05, 0.25, 0.5, 0.75, 0.95]
OBS_YEARS = ('1999', '1997')
_NUMERIC = {'x', 'y', 'day', 'survival_days', 'burial_mm', 'first_grounding_t',
            'first_grounding_x', 'first_grounding_y', 'first_grounding_z', 'terminal_z'}


# ---------------------------------------------------------------------------
# Observed data and the stand
# ---------------------------------------------------------------------------

def obs_summary_path(year):
    return os.path.join(SCRIPT_DIR, f'juvenile_zone_species_summary_{year}_t1_5.csv')


def zone_shares_from_densities(dens):
    """{species: {zone: share}} from {(zone, species): density}."""
    out = {}
    for sp in SPECIES:
        tot = sum(dens.get((z, sp), 0.0) for z in ZONES)
        out[sp] = {z: (dens.get((z, sp), 0.0) / tot if tot > 0 else 0.0) for z in ZONES}
    return out


def obs_quadrat_path(year):
    return os.path.join(SCRIPT_DIR, f'juvenile_quadrat_records_{year}_t1_5.csv')


def load_observed_quadrats(year):
    """Per species: the recorded per-quadrat DENSITIES (#/m2) of one survey year, by zone.
    2026-09-20 (data audit): the survey has no verified counts (the former `count_per_quadrat`
    was density x an assumed 2.5 m2 and is gone), so no dispersion factor phi and no
    chi-square are derived here; the densities feed the bootstrap statistics. Returns
      densities            list of (zone, density) per quadrat
      mean_density_zone    {zone: mean recorded density}
      mean_density_all     mean recorded density over ALL quadrats of the year
      n_quadrats           {zone: n}
      n_nonzero_quadrats   occupied quadrats
      testable             n_nonzero_quadrats >= MIN_NONZERO_QUADRATS (a zone profile of
                           fewer occupied quadrats is not compared; Sonneratia: 3 in 1999)
      note                 '' or why the species is not testable
    Rows of species outside SPECIES (Bruguiera) are ignored."""
    path = obs_quadrat_path(year)
    if not os.path.exists(path):
        return None
    per = {sp: [] for sp in SPECIES}
    # 2026-09-21 FIX: 26 of the 335 quadrats (all on the diagonal T5) lie OUTSIDE the model box, so
    # the model can never place a recruit in them. Scoring them against the model biased the target:
    # 2-3 share points for the calibrated species and 39 points for Sonneratia in 1997, whose only
    # zone C record is one of them. The observation is restricted to the domain the model represents.
    import tidal_dispersal_model as _t
    n_out = 0
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            sp = r['species']
            x, y = float(r['x']), float(r['y'])
            if not (_t.X_MIN <= x <= _t.X_MAX and _t.Y_MIN <= y <= _t.Y_MAX):
                n_out += 1
                continue
            if sp in per:
                per[sp].append((r['model_zone'], float(r['density_per_m2'])))
    out = {}
    for sp, rows in per.items():
        d = np.asarray([v for _z, v in rows], dtype=float)
        z = np.asarray([zz for zz, _v in rows])
        n_nonzero = int((d > 0).sum())
        testable = n_nonzero >= MIN_NONZERO_QUADRATS
        out[sp] = {'densities': rows,
                   'mean_density_zone': {zz: (float(d[z == zz].mean()) if (z == zz).any() else 0.0) for zz in ZONES},
                   'mean_density_all': float(d.mean()) if len(d) else float('nan'),
                   'n_quadrats': {zz: int((z == zz).sum()) for zz in ZONES},
                   'n_nonzero_quadrats': n_nonzero,
                   'testable': testable,
                   'note': ('' if testable else
                            f'only {n_nonzero} occupied quadrats (< {MIN_NONZERO_QUADRATS}): zone profile not compared')}
    return out


def mean_density_all_quadrats(year, species=None, n_boot=4000, seed=42):
    """Mean recorded density (#/m2) over ALL quadrats of the survey year per species, with a
    95 % bootstrap interval (quadrats resampled within zones, the survey strata). This is
    the supply-calibration statistic for Sonneratia (2026-09-20): it needs neither quadrat
    areas nor counts. Returns {species: (mean, lo, hi, n_quadrats, n_nonzero)}."""
    oq = load_observed_quadrats(year)
    if oq is None:
        return None
    rng = np.random.default_rng(seed)
    out = {}
    for sp in (species or SPECIES):
        rows = oq[sp]['densities']
        z = np.asarray([zz for zz, _v in rows])
        d = np.asarray([v for _z, v in rows], dtype=float)
        idx = {zz: np.where(z == zz)[0] for zz in ZONES if (z == zz).any()}
        boots = np.empty(n_boot)
        for b in range(n_boot):
            pick = np.concatenate([rng.choice(ix, len(ix)) for ix in idx.values()])
            boots[b] = d[pick].mean()
        out[sp] = (float(d.mean()), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)),
                   int(len(d)), int((d > 0).sum()))
    return out


def sonneratia_mean_density(year, n_boot=4000, seed=42):
    """Sonneratia's mean recorded density over all quadrats of the year with its bootstrap
    interval: (mean, lo, hi, n_quadrats, n_nonzero); None if Sonneratia is not modelled."""
    if 'Sonneratia' not in SPECIES:
        return None
    return mean_density_all_quadrats(year, ['Sonneratia'], n_boot, seed)['Sonneratia']


def load_observed(year):
    """Observed zone x species mean densities, relative densities across zones ('shares'),
    quadrat counts, species fractions, the per-quadrat densities and the mean density over
    all quadrats per species (`mean_density_all_quadrats`) of one survey year."""
    dens, nq = {}, {}
    with open(obs_summary_path(year), newline='') as f:
        for r in csv.DictReader(f):
            if r['species'] in SPECIES:   # Bruguiera rows (not modelled) are skipped
                dens[(r['zone'], r['species'])] = float(r['mean_density'])
                nq[r['zone']] = int(r['n_quadrats'])
    quadrats = load_observed_quadrats(year)
    return {'densities': dens, 'shares': zone_shares_from_densities(dens), 'n_quadrats': nq,
            'species_fraction': species_fraction_from_densities(dens),
            'quadrats': quadrats,
            'mean_density_all_quadrats': ({sp: quadrats[sp]['mean_density_all'] for sp in SPECIES}
                                          if quadrats else None)}


def species_fraction_from_densities(dens):
    """{species: share of the summed species mean densities}: sum over zones of the species'
    mean density per zone, divided by the same sum over all species (zones weighted EQUALLY,
    regardless of their quadrat numbers or areas; it is not a fraction of individuals). The
    observed benchmark of S4 (2026-09-20). For the supply calibration use
    mean_density_all_quadrats instead (quadrat-weighted, no zone weighting)."""
    tot_sp = {sp: sum(dens.get((z, sp), 0.0) for z in ZONES) for sp in SPECIES}
    tot = sum(tot_sp.values())
    return {sp: (tot_sp[sp] / tot if tot > 0 else float('nan')) for sp in SPECIES}


NO_COUNTS_NOTE = 'no verified counts (recorded densities only); use TVD/JSD'


def chi2_gof(obs_zone_counts, sim_zone_weights, phi=1.0, min_expected=1.0, testable=True, note=None):
    """Chi-square goodness of fit of observed zone COUNTS of one species to simulated zone
    proportions. NOT APPLIED TO THE FIELD DATA since 2026-09-20 (the survey records densities
    of unknown quadrat area; no verified counts exist; see chi2_tests). Kept as a generic
    utility for genuinely counted data.

    obs_zone_counts   observed individuals per zone
    sim_zone_weights  simulated recruits per zone in the same sampled area (relative)
    phi               overdispersion factor of the quadrat counts (quasi-likelihood
                      correction: chi2_adj = chi2 / phi)
    Expected counts below min_expected are raised to min_expected (a zone the model leaves
    empty still enters the statistic); df = 2 (three zones).
    testable=False (a species with too few occupied quadrats, 2026-09-20) keeps the
    statistic but reports p = NaN with a note; phi <= 0 is treated as 1.
    """
    from scipy.stats import chi2 as chi2_dist
    n = float(sum(obs_zone_counts.get(z, 0.0) for z in ZONES))
    w = {z: max(float(sim_zone_weights.get(z, 0.0)), 0.0) for z in ZONES}
    tot = sum(w.values())
    if n <= 0 or tot <= 0:
        return None
    phi = float(phi) if phi and phi > 0 else 1.0
    exp = {z: max(n * w[z] / tot, min_expected) for z in ZONES}
    obs = {z: float(obs_zone_counts.get(z, 0.0)) for z in ZONES}
    chi2 = sum((obs[z] - exp[z]) ** 2 / exp[z] for z in ZONES)
    df = len(ZONES) - 1
    adj = chi2 / phi
    out = {'chi2': float(chi2), 'phi': float(phi), 'chi2_adj': float(adj), 'df': int(df),
           'p': float(chi2_dist.sf(adj, df)) if testable else float('nan'),
           'n': n, 'expected': exp, 'observed': obs, 'cells': list(ZONES)}
    if not testable:
        out['note'] = note or 'too few occupied quadrats: p not computed'
    return out


def p_value(test):
    """p of a chi2_gof result, NaN when the test is missing or not computable (JSON None)."""
    if not test or test.get('p') is None:
        return float('nan')
    return float(test['p'])


def chi2_tests(sampled_densities, n_quadrats, obs_quadrats):
    """Per-species observed-vs-simulated chi-square: DROPPED for the observed data
    (2026-09-20). The survey records per-quadrat densities of unknown, possibly unequal
    quadrat area; the counts the test needs were manufactured (density x 2.5 m2) and are no
    longer produced. The signature is kept so every caller still works: each species gets
    {'p': None, 'note': NO_COUNTS_NOTE, 'n_nonzero_quadrats': ..., 'testable': ...} and
    p_value() returns NaN. Use the density-based zone shares (JSD 'zonation_error') and the
    bootstrap TVD (equivalence_test.py, virtual_survey.py) instead."""
    out = {}
    for sp in SPECIES:
        oq = (obs_quadrats or {}).get(sp)
        out[sp] = ({'p': None, 'note': NO_COUNTS_NOTE,
                    'n_nonzero_quadrats': oq.get('n_nonzero_quadrats'),
                    'testable': oq.get('testable')} if oq else None)
    return out


def load_all_observed(years=OBS_YEARS):
    return {str(y): load_observed(y) for y in years if os.path.exists(obs_summary_path(y))}


def zone_of_z(z):
    return 'A' if z >= tdm.ZONE_A_MIN else ('B' if z >= tdm.ZONE_B_MIN else 'C')


def stand_summary(stand_csv):
    """Parent-tree zone shares and elevation distribution per species (stage 'parents')."""
    trees = tdm.load_trees(stand_csv)
    zvals = {sp: [] for sp in SPECIES}
    zones = {sp: Counter() for sp in SPECIES}
    for t in trees:
        if t['species'] not in zvals:
            continue
        zvals[t['species']].append(tdm.elevation(t['x'], t['y']))
        zones[t['species']][tdm.zone_at(t['x'], t['y'])] += 1
    st = _stage(zvals, zones)
    st['csv'] = stand_csv
    st['n_trees'] = len(trees)
    return st


# ---------------------------------------------------------------------------
# Small numerics
# ---------------------------------------------------------------------------

def jsd2(p, q):
    """Jensen-Shannon divergence (base 2, in [0, 1]) of two discrete distributions."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.sum() <= 0 or q.sum() <= 0:
        return float('nan')
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)

    def _kl(a, b):
        mask = a > 0
        return float(np.sum(a[mask] * np.log2(a[mask] / b[mask])))

    return 0.5 * _kl(p, m) + 0.5 * _kl(q, m)


def _float(v):
    if v is None or v == '':
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _separation(zvals):
    """Pairwise W1 (m) and histogram JSD between species elevation samples."""
    out = {}
    w1s = []
    for a, b in PAIRS:
        za, zb = np.asarray(zvals.get(a, [])), np.asarray(zvals.get(b, []))
        if len(za) >= 2 and len(zb) >= 2:
            w1 = float(wasserstein_distance(za, zb))
            ha, _ = np.histogram(za, bins=Z_EDGES)
            hb, _ = np.histogram(zb, bins=Z_EDGES)
            js = jsd2(ha, hb)
        else:
            w1, js = float('nan'), float('nan')
        out[f'{a[:3]}-{b[:3]}'] = {'w1': w1, 'jsd': js}
        w1s.append(w1)
    out['mean_w1'] = float(np.nanmean(w1s)) if np.any(~np.isnan(w1s)) else float('nan')
    return out


def _stage(zvals, zone_counts):
    """Stage record: species zone counts/shares, elevation histograms, quantiles, separation."""
    rec = {'n': {}, 'zone_counts': {}, 'zone_shares': {}, 'hist': {}, 'quantiles': {}, 'mean': {}}
    for sp in SPECIES:
        z = np.asarray(zvals.get(sp, []), dtype=float)
        zc = zone_counts.get(sp, Counter())
        tot = sum(zc.get(zz, 0) for zz in ZONES)
        rec['n'][sp] = int(len(z))
        rec['zone_counts'][sp] = {zz: int(zc.get(zz, 0)) for zz in ZONES}
        rec['zone_shares'][sp] = {zz: (zc.get(zz, 0) / tot if tot > 0 else 0.0) for zz in ZONES}
        if len(z):
            h, _ = np.histogram(z, bins=Z_EDGES)
            rec['hist'][sp] = [int(c) for c in h]
            rec['quantiles'][sp] = [float(q) for q in np.quantile(z, QUANTILES)]
            rec['mean'][sp] = float(z.mean())
        else:
            rec['hist'][sp] = [0] * (len(Z_EDGES) - 1)
            rec['quantiles'][sp] = [float('nan')] * len(QUANTILES)
            rec['mean'][sp] = float('nan')
    rec['separation'] = _separation(zvals)
    return rec


# ---------------------------------------------------------------------------
# Observation operator (validation-pipeline convention)
# ---------------------------------------------------------------------------

def transect_membership(results):
    """{id(result): True/False} for grounded propagules (established or dead at their
    terminal position) inside the survey transect polygons."""
    member = {}
    for r in results:
        if r.get('outcome') in ('established', 'dead'):
            member[id(r)] = tg.assign_to_transect(float(r['x']), float(r['y'])) is not None
    return member


def sampled_densities(established, member, n_quadrats):
    """Zone x species density of established seedlings inside the transect polygons per
    ASSUMED sampled footprint (n_quadrats x transect_geometry.MODEL_QUADRAT_FOOTPRINT_M2;
    transect_geometry.model_zone_sampled_densities convention). Zone shares derived from
    these densities are footprint-free; the absolute values are conditional on the
    footprint assumption (2026-09-20)."""
    counts = Counter()
    for r in established:
        if member.get(id(r), False):
            counts[(r['zone'], r['species'])] += 1
    dens = {}
    for z in ZONES:
        area = n_quadrats.get(z, 0) * tg.MODEL_QUADRAT_FOOTPRINT_M2
        for sp in SPECIES:
            dens[(z, sp)] = counts[(z, sp)] / area if area > 0 else 0.0
    return dens


# ---------------------------------------------------------------------------
# Run summary
# ---------------------------------------------------------------------------

def _summarise_stratum(sub, member, n_quadrats, stand, obs):
    n = len(sub)
    s = {'n': n,
         'by_outcome': dict(Counter(r.get('outcome') for r in sub)),
         'by_settle_method': dict(Counter(r.get('settle_method') for r in sub)),
         'by_species': dict(Counter(r.get('species') for r in sub))}
    s['clamp_frac'] = (s['by_settle_method'].get('boundary_clamp', 0) / n) if n else float('nan')

    grounded = [r for r in sub if r.get('outcome') in ('established', 'dead')]
    established = [r for r in sub if r.get('outcome') == 'established']

    fg_z = {sp: [] for sp in SPECIES}
    fg_zone = {sp: Counter() for sp in SPECIES}
    gr_z = {sp: [] for sp in SPECIES}
    gr_zone = {sp: Counter() for sp in SPECIES}
    for r in grounded:
        sp = r['species']
        if sp not in fg_z:
            continue
        z = _float(r.get('first_grounding_z'))
        if z is not None:
            fg_z[sp].append(z)
            fg_zone[sp][zone_of_z(z)] += 1
        tz = _float(r.get('terminal_z'))
        if tz is not None:
            gr_z[sp].append(tz)
        gr_zone[sp][r['zone']] += 1
    es_z = {sp: [] for sp in SPECIES}
    es_zone = {sp: Counter() for sp in SPECIES}
    for r in established:
        sp = r['species']
        if sp not in es_z:
            continue
        tz = _float(r.get('terminal_z'))
        if tz is not None:
            es_z[sp].append(tz)
        es_zone[sp][r['zone']] += 1

    s['stages'] = {
        'first_grounding': _stage(fg_z, fg_zone),
        'grounded': _stage(gr_z, gr_zone),
        'established': _stage(es_z, es_zone),
    }

    dens = sampled_densities(established, member, n_quadrats)
    shares = zone_shares_from_densities(dens)
    s['sampled_densities'] = {f'{z}_{sp}': dens[(z, sp)] for z in ZONES for sp in SPECIES}
    s['sampled_shares'] = shares
    s['n_sampled'] = int(sum(1 for r in established if member.get(id(r), False)))
    # the stranded (terminal) cohort through the same operator, and the survival that
    # turns it into the recruit distribution
    dens_gr = sampled_densities(grounded, member, n_quadrats)
    s['sampled_densities_grounded'] = {f'{z}_{sp}': dens_gr[(z, sp)] for z in ZONES for sp in SPECIES}
    s['sampled_shares_grounded'] = zone_shares_from_densities(dens_gr)
    s['n_sampled_grounded'] = int(sum(1 for r in grounded if member.get(id(r), False)))
    s['survival_zone'] = {sp: {z: (dens[(z, sp)] / dens_gr[(z, sp)] if dens_gr[(z, sp)] > 0 else None)
                               for z in ZONES} for sp in SPECIES}
    # elevation histograms of the same transect-sampled cohorts (survey footprint)
    s['hist_sampled'] = {}
    for stage_name, cohort in (('grounded', grounded), ('established', established)):
        hs = {}
        for sp in SPECIES:
            zs = [_float(r.get('terminal_z')) for r in cohort
                  if r['species'] == sp and member.get(id(r), False)]
            zs = [z for z in zs if z is not None]
            h, _ = np.histogram(np.asarray(zs, dtype=float), bins=Z_EDGES) if zs else (np.zeros(len(Z_EDGES) - 1), None)
            hs[sp] = [int(c) for c in h]
        s['hist_sampled'][stage_name] = hs

    # landward shift of recruits relative to the parents of this run
    if stand is not None:
        s['shift_sampled'] = {sp: shares[sp]['A'] - stand['zone_shares'][sp]['A'] for sp in SPECIES}
        s['shift_domain'] = {sp: s['stages']['established']['zone_shares'][sp]['A']
                             - stand['zone_shares'][sp]['A'] for sp in SPECIES}
    else:
        s['shift_sampled'] = s['shift_domain'] = None

    # signatures (model-only). 2026-09-20: 4 species (Sonneratia) -- S1-S3 rank the CALIBRATED
    # species only (otherwise S2 fails by construction whenever Sonneratia, the smallest
    # propagule, is the most seaward); S4 reports Sonneratia and is not gated.
    a_share = {sp: shares[sp]['A'] for sp in SPECIES}
    rank = sorted(CALIBRATED_SPECIES, key=lambda sp: a_share[sp], reverse=True)
    rhi_peak = max(ZONES, key=lambda z: shares['Rhizophora'][z])
    s['landward_rank'] = rank
    s['landward_rank_all'] = sorted(SPECIES, key=lambda sp: a_share[sp], reverse=True)
    shift_cal = ({sp: s['shift_sampled'][sp] for sp in CALIBRATED_SPECIES} if s['shift_sampled'] else None)
    s['signatures'] = {
        'S1_ceriops_most_landward': rank[0] == 'Ceriops',
        'S2_rhizophora_most_seaward_peak_B': rank[-1] == 'Rhizophora' and rhi_peak == 'B',
        'S3_avicennia_largest_shift': (
            max(shift_cal, key=shift_cal.get) == 'Avicennia' if shift_cal else None),
    }
    s['signatures'].update(_sonneratia_signature(established, member, dens))

    # pattern scores against each survey year
    s['scores'] = {}
    for yr, ob in obs.items():
        share_jsd = {sp: jsd2([shares[sp][z] for z in ZONES], [ob['shares'][sp][z] for z in ZONES])
                     for sp in SPECIES}
        comp = v99.zone_composition_jsd(dens, ob['densities'])
        peak = v99.species_peak_zone_match(dens, ob['densities'])
        rank_z = v99.within_zone_species_rank(dens, ob['densities'])
        share_jsd_gr = {sp: jsd2([s['sampled_shares_grounded'][sp][z] for z in ZONES],
                                 [ob['shares'][sp][z] for z in ZONES]) for sp in SPECIES}
        # zonation error = mean JSD over the calibrated species (the calibration objective,
        # comparable with the 3-species runs); '_all' = mean over every modelled species
        # (NaN-safe: a species absent from a survey has no observed shares)
        s['scores'][yr] = {
            'chi2': chi2_tests(s['sampled_densities'], n_quadrats, ob.get('quadrats')),
            'chi2_grounded': chi2_tests(s['sampled_densities_grounded'], n_quadrats, ob.get('quadrats')),
            'zonation_error_grounded': _mean_over(share_jsd_gr, CALIBRATED_SPECIES),
            'zonation_error_grounded_all': _mean_over(share_jsd_gr, SPECIES),
            'spearman_rho': float(v99.zone_spearman_rho(dens, ob['densities'])),
            'spearman_rho_all': float(v99.zone_spearman_rho(dens, ob['densities'], species=SPECIES)),
            'peak_matches': int(peak['matches']),
            'peak_details': {sp: d for sp, d in peak.get('details', {}).items()},
            'within_zone_rank_agreeing': int(rank_z['n_agreeing_zones']),
            'composition_jsd': {z: float(comp['per_zone'][z]['jsd']) for z in ZONES},
            'composition_jsd_mean': float(comp['mean_jsd']),
            'share_jsd': share_jsd,
            'zonation_error': _mean_over(share_jsd, CALIBRATED_SPECIES),
            'zonation_error_all': _mean_over(share_jsd, SPECIES),
            'species_fraction_obs': ob.get('species_fraction'),
        }

    # mortality guardrails (2026-09-20: over the CALIBRATED species; Sonneratia's uncalibrated
    # seedlings die at ~100 % in zones B and C and would fail Z3 by construction; the pooled
    # 4-species numbers are kept under 'pct_all')
    sub_cal = [r for r in sub if r['species'] in CALIBRATED_SPECIES]
    mort = v99.zone_mortality_gradient(sub_cal)
    cause = v99.mortality_cause_analysis(sub_cal)
    mort_all = v99.zone_mortality_gradient(sub)
    s['mortality'] = {
        'pct': mort['mortality'], 'Z3_pass': bool(mort['pass']),
        'ordering_correct': bool(mort['ordering_correct']),
        'cause_C': {k: cause['C'][k] for k in ('siltation', 'predation', 'washed_out')},
        'Z6_pass': bool(cause['C']['predation'] > cause['C']['siltation']),
        'pct_all': mort_all['mortality'], 'Z3_pass_all': bool(mort_all['pass']),
        'species': list(CALIBRATED_SPECIES),
    }
    return s


def _mean_over(d, species):
    """NaN-ignoring mean of d[sp] over the given species (NaN if none is finite)."""
    vals = [d[sp] for sp in species if sp in d and d[sp] is not None and not math.isnan(d[sp])]
    return float(np.mean(vals)) if vals else float('nan')


def _sonneratia_signature(established, member, dens, sp='Sonneratia'):
    """S4 (reported, not gated; 2026-09-20): Sonneratia's share of the recruits and the
    10th/90th percentiles of its terminal elevation.
      S4_sonneratia_recruit_fraction  density-weighted share under the zone-share operator:
                                      sum over zones of Sonneratia's sampled density / the same
                                      sum over all species (zones weighted equally, like the
                                      observed 'species_fraction' benchmark)
      S4_sonneratia_count_fraction    plain count share of the established recruits inside the
                                      transect polygons
      S4_sonneratia_elev_p10_p90      [p10, p90] of the terminal elevation (m) of Sonneratia
                                      recruits inside the transects; '_domain' = whole domain
    Observed benchmark: < 1 % of the juvenile density; quadrat elevations 1.37-1.66 m.
    All values NaN/None when the species is not modelled or has no recruits."""
    if sp not in SPECIES:
        return {}
    tot_sp = sum(dens.get((z, sp), 0.0) for z in ZONES)
    tot = sum(dens.get((z, s2), 0.0) for z in ZONES for s2 in SPECIES)
    n_in = sum(1 for r in established if member.get(id(r), False))
    n_sp = sum(1 for r in established if member.get(id(r), False) and r['species'] == sp)
    z_in = [_float(r.get('terminal_z')) for r in established if r['species'] == sp and member.get(id(r), False)]
    z_in = [z for z in z_in if z is not None]
    z_dom = [_float(r.get('terminal_z')) for r in established if r['species'] == sp]
    z_dom = [z for z in z_dom if z is not None]
    nan2 = [float('nan'), float('nan')]
    return {
        'S4_sonneratia_recruit_fraction': (tot_sp / tot) if tot > 0 else float('nan'),
        'S4_sonneratia_count_fraction': (n_sp / n_in) if n_in > 0 else float('nan'),
        'S4_sonneratia_n_sampled': int(n_sp),
        'S4_sonneratia_elev_p10_p90': [float(q) for q in np.percentile(z_in, [10, 90])] if z_in else nan2,
        'S4_sonneratia_elev_p10_p90_domain': [float(q) for q in np.percentile(z_dom, [10, 90])] if z_dom else nan2,
    }


def summarise_run(results, stand_csv=None, obs=None, strata=STRATA, meta=None):
    """Summary dict of one run; JSON-serialisable."""
    obs = obs if obs is not None else load_all_observed()
    stand = stand_summary(stand_csv) if stand_csv else None
    n_quadrats = tg.load_observed_quadrats_per_zone()
    member = transect_membership(results)
    out = {'meta': meta or {}, 'n_results': len(results), 'n_quadrats': n_quadrats,
           'stand': stand, 'strata': {}}
    for name, flt in strata.items():
        sub = results if flt is None else [r for r in results if r.get('source') == flt]
        out['strata'][name] = _summarise_stratum(sub, member, n_quadrats, stand, obs)
    return _clean(out)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if (isinstance(o, float) and math.isnan(o)) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


# ---------------------------------------------------------------------------
# Calibration objective
# ---------------------------------------------------------------------------

def score_summary(summary, stratum='all', w1997=0.5, penalty=0.2, clamp_max=0.05):
    """Objective (lower is better) and its components for one run summary."""
    s = summary['strata'][stratum]
    e99 = s['scores'].get('1999', {}).get('zonation_error')
    e97 = s['scores'].get('1997', {}).get('zonation_error')
    sig = s['signatures']
    comps = {
        'zonation_error_1999': e99,
        'zonation_error_1997': e97,
        # 2026-09-20: the objective stays on the calibrated species; the all-species means are reported only
        'zonation_error_1999_all': s['scores'].get('1999', {}).get('zonation_error_all'),
        'zonation_error_1997_all': s['scores'].get('1997', {}).get('zonation_error_all'),
        'S4_sonneratia_recruit_fraction': sig.get('S4_sonneratia_recruit_fraction'),
        'S1S2_ok': bool(sig['S1_ceriops_most_landward'] and sig['S2_rhizophora_most_seaward_peak_B']),
        'S3_ok': bool(sig['S3_avicennia_largest_shift']),
        'Z3_ok': bool(s['mortality']['Z3_pass']),
        'clamp_ok': bool(s['clamp_frac'] is not None and s['clamp_frac'] <= clamp_max),
        'clamp_frac': s['clamp_frac'],
        'spearman_1999': s['scores'].get('1999', {}).get('spearman_rho'),
        'mean_w1_established': s['stages']['established']['separation']['mean_w1'],
    }
    obj = (e99 if e99 is not None else 1.0) + w1997 * (e97 if e97 is not None else 1.0)
    obj += penalty * sum(not comps[k] for k in ('S1S2_ok', 'S3_ok', 'Z3_ok', 'clamp_ok'))
    comps['objective'] = float(obj)
    return comps


def compact_line(summary, stratum='all'):
    s = summary['strata'][stratum]
    sh = s['sampled_shares']
    parts = [f"{sp[:3]} A/B/C={sh[sp]['A']:.2f}/{sh[sp]['B']:.2f}/{sh[sp]['C']:.2f}" for sp in SPECIES if sp in sh]
    sc = score_summary(summary, stratum)
    chi = s['scores'].get('1999', {}).get('chi2') or {}
    # 2026-09-20: p is None for the observed data (no verified counts); printed as n/a
    pv = ' '.join(f"p{sp[:3]}={p_value(chi.get(sp)):.2f}" if not math.isnan(p_value(chi.get(sp))) else f"p{sp[:3]}=n/a"
                  for sp in SPECIES)
    e99 = sc['zonation_error_1999'] if sc['zonation_error_1999'] is not None else float('nan')
    e97 = sc['zonation_error_1997'] if sc['zonation_error_1997'] is not None else float('nan')
    return (', '.join(parts) + f" | {pv} | err99={e99:.3f} err97={e97:.3f}"
            f" rho99={sc['spearman_1999']:.2f} S1S2={int(sc['S1S2_ok'])} S3={int(sc['S3_ok'])}"
            f" Z3={int(sc['Z3_ok'])} clamp={sc['clamp_frac']:.3f} W1={sc['mean_w1_established']:.3f}"
            f" obj={sc['objective']:.3f}")


# ---------------------------------------------------------------------------
# Cache IO
# ---------------------------------------------------------------------------

def load_results_csv(path):
    results = []
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            for k in _NUMERIC:
                v = row.get(k, '')
                if v != '':
                    try:
                        row[k] = float(v)
                    except (TypeError, ValueError):
                        pass
            results.append(row)
    return results


def summary_path_for(cache_path):
    return os.path.splitext(cache_path)[0] + '_summary.json'


def summarise_cache(cache_path, stand_csv, out_path=None, force=False, meta=None, obs=None):
    out_path = out_path or summary_path_for(cache_path)
    if os.path.exists(out_path) and not force:
        with open(out_path) as f:
            return json.load(f)
    results = load_results_csv(cache_path)
    meta = dict(meta or {})
    meta.setdefault('cache', os.path.basename(cache_path))
    summary = summarise_run(results, stand_csv, obs=obs, meta=meta)
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=1)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cache', help='one cached run CSV')
    ap.add_argument('--cache-dir', help='summarise every *.csv run cache in this directory')
    ap.add_argument('--stand', required=True, help='tree inventory CSV used by the run(s)')
    ap.add_argument('--out', default=None, help='summary JSON path (single cache only)')
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    paths = [args.cache] if args.cache else sorted(
        p for p in glob.glob(os.path.join(args.cache_dir, '*.csv')) if not p.endswith('_trees.csv'))
    obs = load_all_observed()
    for p in paths:
        summ = summarise_cache(p, args.stand, args.out if args.cache else None, force=args.force, obs=obs)
        print(os.path.basename(p), '|', compact_line(summ))


if __name__ == '__main__':
    main()
