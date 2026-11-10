"""Simulated recruitment versus observed ADULT zonation on the withheld transects (2026-09-20).

Second evaluation pattern of the model (the first is the juvenile zonation). The SAME simulated
recruits (established propagules of the replicate runs) are surveyed with the point-centred
quarter protocol of the adult survey: at every PCQM sample point of the withheld transects T3
and T1 (and, separately, of the orthogonal T5, never used by any reconstruction), the nearest
recruit in each quarter within the search limit R is recorded (`reconstruction/virtual_pcqm.
survey_transect`, marks ignored), quarters are zoned by the position of their plot point, and
the share of each species' nearest-individual records per zone is compared with the same
quantity computed from the real adult records of the same points and quarters. Sampling support
is therefore matched by construction (same points, quarters, search limit, empty-quarter rule);
the quantity is frequency-based (share of a species' nearest-individual records per zone) and
is NOT the relative density used for the juvenile comparison.

Uncertainty: adults by bootstrap over plots within zones; recruits over the replicate runs
(and, when present, the reconstruction-ensemble replicates) with plots resampled the same way.
Reported: zone shares, the recruit-minus-adult difference per zone with its 95 % interval
(direction and size as found), the total variation distance (TVD) per species with its 95 %
interval and the per-zone JSD, as measures of spatial correspondence with the standing forest.
No pass/fail threshold is applied. Agreement is not simulated adult development.

Secondary operator (table only): recruit counts inside the 10 m transect belt polygons.

Usage:
    python evaluate_recruits_vs_adults.py [--caches "results/knockout_cascade/cascade_D1T1S1E1_seed*.csv"]
                                          [--stand-caches "results/knockout_cascade/stand*_seed*.csv"]
                                          [--withheld 3 1] [--n-boot 4000]
Outputs:
    results/evaluation/recruits_vs_adults_tests.csv, recruits_vs_adults_zone_composition.csv,
    recruits_vs_adults_belt_counts.csv, results/figures/fig_recruits_vs_adults.png
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
sys.path.insert(0, os.path.join(SCRIPT_DIR, 'reconstruction'))
import recon_config as C  # noqa: E402
import recon_geometry as G  # noqa: E402
import transect_geometry as tg  # noqa: E402
import prep_pcqm  # noqa: E402
import virtual_pcqm as V  # noqa: E402
from tidal_dispersal_model import SPECIES, zone_at  # noqa: E402

SPECIES = list(SPECIES)
ZONES = ['C', 'B', 'A']
ZONE_FACE = {'C': '#d0e3f3', 'B': '#adb27c', 'A': '#4b4e36'}
ZONE_LABEL = {'C': 'Zone C (seaward)', 'B': 'Zone B', 'A': 'Zone A (landward)'}
BAR_EDGE = '#333333'
SP_COLORS = {'Avicennia': '#E07A5F', 'Rhizophora': '#3D8B7A', 'Ceriops': '#7B2D8E', 'Sonneratia': '#3B6CC0'}

DEFAULT_CACHES = os.path.join(SCRIPT_DIR, 'results', 'knockout_cascade', 'cascade_D1T1S1E1_seed*.csv')
DEFAULT_STAND_CACHES = os.path.join(SCRIPT_DIR, 'results', 'knockout_cascade', 'stand*_seed*.csv')


# --- data ------------------------------------------------------------------------------------

def load_real(transects, frames):
    """Real PCQM records of the given transects with positions, zone and in-domain flag."""
    rec = prep_pcqm.load_records(transects=transects)
    rec = prep_pcqm.add_positions(rec, frames)
    other = rec['present'] & ~rec['species'].isin(C.STAND_SPECIES)   # Bruguiera: absent from the stand
    rec.loc[other, ['present', 'd', 'species']] = [False, np.nan, None]
    rec.loc[other, 'd_cap'] = C.DMAX
    x, y = G.to_model(rec['px'].values, rec['py'].values)
    rec['x'], rec['y'] = x, y
    rec['in_domain'] = G.in_domain(x, y)
    rec['zone'] = [zone_at(a, b) for a, b in zip(x, y)]
    return rec[rec['in_domain']].reset_index(drop=True)


def survey(stand, transects, frames, real):
    """Virtual PCQM survey of a stand along the transects, zoned like the real records."""
    parts = []
    zone_by = {(int(t), int(p)): z for t, p, z in zip(real['transect'], real['plot'], real['zone'])}
    for t in transects:
        v = V.survey_transect(stand, t, frames)
        v = v[[(int(t), int(p)) in zone_by for p in v['plot']]].copy()
        v['zone'] = [zone_by[(int(t), int(p))] for p in v['plot']]
        parts.append(v)
    return pd.concat(parts, ignore_index=True)


def species_zone_counts(rec):
    pres = rec[rec['present'] == True]
    return pd.crosstab(pres['species'], pres['zone']).reindex(index=SPECIES, columns=ZONES, fill_value=0).values.astype(float)


def shares(counts):
    tot = counts.sum(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(tot > 0, counts / np.where(tot > 0, tot, 1), np.nan)


def bootstrap_counts(recs, rng, n_boot):
    """Resample plots within zones (each plot = 4 quarters) for one record set; returns
    (n_boot, n_species, n_zones) counts."""
    plots = recs.drop_duplicates(['transect', 'plot'])[['transect', 'plot', 'zone']].reset_index(drop=True)
    row_of = {(int(t), int(p)): i for i, (t, p) in enumerate(zip(plots['transect'], plots['plot']))}
    mat = np.zeros((len(plots), len(SPECIES)))
    pres = recs[recs['present'] == True]
    sp_idx = {sp: i for i, sp in enumerate(SPECIES)}
    for t, p, sp in zip(pres['transect'], pres['plot'], pres['species']):
        if sp in sp_idx:
            mat[row_of[(int(t), int(p))], sp_idx[sp]] += 1
    by_zone = {z: np.flatnonzero(plots['zone'].values == z) for z in ZONES}
    out = np.zeros((n_boot, len(SPECIES), len(ZONES)))
    for j, z in enumerate(ZONES):
        rows = by_zone[z]
        if len(rows) == 0:
            continue
        pick = rows[rng.integers(0, len(rows), (n_boot, len(rows)))]
        out[:, :, j] = mat[pick].sum(axis=1)
    return out


def evaluate_set(name, transects, real, stands, frames, rng, n_boot, margin):
    """Observed vs simulated species zone shares for one transect set and a list of stands
    (replicates).  Returns (rows, profiles, summary rows, vs dict for the figure)."""
    rr = real[real['transect'].isin(transects)]
    obs_counts = species_zone_counts(rr)
    obs_share = shares(obs_counts)
    obs_boot = shares(bootstrap_counts(rr, rng, n_boot).reshape(-1, len(ZONES))).reshape(n_boot, len(SPECIES), len(ZONES))
    virts, sim_shares, sim_boots, summ = [], [], [], []
    for k, (label, stand) in enumerate(stands):
        vv = survey(stand, transects, frames, rr)
        virts.append(vv)
        sim_shares.append(shares(species_zone_counts(vv)))
        nb = max(1, n_boot // len(stands))
        sim_boots.append(shares(bootstrap_counts(vv, rng, nb).reshape(-1, len(ZONES))).reshape(nb, len(SPECIES), len(ZONES)))
        for t in transects:
            m = V.compare_records(rr[rr['transect'] == t], vv[vv['transect'] == t], {})
            m.update(set=name, transect=t, stand=label)
            summ.append(m)
    sim_share = np.nanmedian(np.array(sim_shares), axis=0)
    sim_lo = np.nanmin(np.array(sim_shares), axis=0)
    sim_hi = np.nanmax(np.array(sim_shares), axis=0)
    sim_boot = np.concatenate(sim_boots, axis=0)
    nb = min(len(sim_boot), n_boot)
    rows = []
    vs = {'obs': {}, 'sim': {}, 'obs_boot': {}, 'sim_boot': {}}
    for i, sp in enumerate(SPECIES):
        ob, sb = obs_boot[:nb, i, :], sim_boot[:nb, i, :]
        ok = ~np.isnan(ob).any(axis=1) & ~np.isnan(sb).any(axis=1)
        tvd_draws = 0.5 * np.abs(ob[ok] - sb[ok]).sum(axis=1)
        tvd = float(0.5 * np.nansum(np.abs(obs_share[i] - sim_share[i]))) if not np.isnan(obs_share[i]).any() else np.nan
        bound = float(np.percentile(tvd_draws, 95)) if ok.sum() > 20 else np.nan
        jsd_sp = V.jsd(obs_counts[i], np.nan_to_num(sim_share[i])) if obs_counts[i].sum() > 0 else np.nan
        vs['obs'][sp] = (np.nan_to_num(obs_share[i]),)
        vs['sim'][sp] = (np.nan_to_num(sim_share[i]),)
        vs['obs_boot'][sp] = np.nan_to_num(ob)
        vs['sim_boot'][sp] = np.nan_to_num(sb)
        for j, z in enumerate(ZONES):
            olo, ohi = (np.nanpercentile(ob[:, j], [2.5, 97.5]) if np.isfinite(ob[:, j]).sum() > 20 else (np.nan, np.nan))
            rows.append(dict(set=name, transects=' '.join(f'T{t}' for t in transects), species=sp, zone=z,
                             n_obs_adults=int(obs_counts[i, j]), n_obs_species=int(obs_counts[i].sum()),
                             obs_share=obs_share[i, j], obs_lo=olo, obs_hi=ohi,
                             sim_share=sim_share[i, j], sim_lo=sim_lo[i, j], sim_hi=sim_hi[i, j],
                             tvd=tvd, tvd_upper95=bound, equivalent_at_margin=(bound < margin / 100.0) if np.isfinite(bound) else False,
                             jsd_species=jsd_sp, margin_points=margin, n_stands=len(stands), n_boot=int(ok.sum())))
    # per-zone composition (share of the species among the adults of a zone)
    zone_rows = []
    for j, z in enumerate(ZONES):
        o = obs_counts[:, j]
        s = np.nanmean([np.nan_to_num(species_zone_counts(v)[:, j]) for v in virts], axis=0)
        zone_rows.append(dict(set=name, zone=z, n_obs=int(o.sum()), n_sim_mean=float(s.sum()),
                              jsd_zone=V.jsd(o, s) if o.sum() > 0 and s.sum() > 0 else np.nan,
                              **{f'obs_{sp}': (o[i] / o.sum() if o.sum() else np.nan) for i, sp in enumerate(SPECIES)},
                              **{f'sim_{sp}': (s[i] / s.sum() if s.sum() else np.nan) for i, sp in enumerate(SPECIES)}))
    profiles = []
    for t in transects:
        pr = V.profile_by_plot(rr[rr['transect'] == t], levels=SPECIES)
        pr['transect'], pr['kind'], pr['stand'] = t, 'real', 'field'
        profiles.append(pr.reset_index().rename(columns={'index': 'plot'}))
        for (label, _), vv in zip(stands, virts):
            pv = V.profile_by_plot(vv[vv['transect'] == t], levels=SPECIES)
            pv['transect'], pv['kind'], pv['stand'] = t, 'virtual', label
            profiles.append(pv.reset_index().rename(columns={'index': 'plot'}))
    return rows, zone_rows, pd.concat(profiles, ignore_index=True), pd.DataFrame(summ), vs


def stacked_panel(ax, vs, title, ylabel):
    w, gap = 0.40, 0.03
    x = 1.3 * np.arange(len(SPECIES))
    for i, sp in enumerate(SPECIES):
        for k, (key, boot, name) in enumerate((('obs', vs['obs_boot'][sp], 'Observed'),
                                               ('sim', vs['sim_boot'][sp], 'Simulated'))):
            est = vs[key][sp][0]
            xc = x[i] + (k - 0.5) * (w + gap)
            bottom = 0.0
            for j, z in enumerate(ZONES):
                ax.bar(xc, 100 * est[j], w, bottom=100 * bottom, facecolor=ZONE_FACE[z], edgecolor=BAR_EDGE, lw=0.7, zorder=2)
                bottom += est[j]
            if len(boot):
                cum = np.cumsum(boot, axis=1)[:, :2]
                for b, xb in enumerate((xc - w / 4, xc + w / 4)):
                    lo, hi = np.percentile(cum[:, b], [2.5, 97.5])
                    if hi - lo < 1e-9:
                        continue
                    ax.plot([xb, xb], [100 * lo, 100 * hi], color='k', lw=2.4, zorder=6, solid_capstyle='butt', clip_on=False)
                    for yv in (lo, hi):
                        ax.plot([xb - 0.07, xb + 0.07], [100 * yv] * 2, color='k', lw=2.4, zorder=6, clip_on=False)
            ax.text(xc, -4.5, name, ha='center', va='top', fontsize=7.5, color='#333333')
    ax.set_xticks(x)
    ax.set_xticklabels([f'$\\it{{{sp}}}$' for sp in SPECIES], fontsize=10)
    ax.tick_params(axis='x', pad=22, length=0)
    ax.set_xlim(x[0] - 0.65, x[-1] + 0.65)
    ax.set_ylim(-2.5, 101)
    ax.spines['bottom'].set_position(('data', -2.5))
    ax.set_ylabel(ylabel, fontsize=9.5)
    ax.set_title(title, fontsize=10.5, loc='left', fontweight='bold')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def load_recruits(path):
    """Established recruits of one run as a 'stand' for the virtual PCQM survey (rotated frame)."""
    df = pd.read_csv(path, usecols=['species', 'x', 'y', 'outcome'])
    df = df[df['outcome'] == 'established']
    xp, yp = G.to_rot(df['x'].values, df['y'].values)
    return pd.DataFrame({'xp': xp, 'yp': yp, 'species': df['species'].values,
                         'r_stem': 0.0, 'h_stem': 0.0, 'x': df['x'].values, 'y': df['y'].values})


def belt_counts(recruits, transects):
    """Secondary operator: recruits per species x zone inside the 10 m belt polygons."""
    tid = tg.assign_to_transect_batch(recruits['x'].values, recruits['y'].values, ids=list(transects))
    sub = recruits[np.isin(tid, list(transects))].copy()
    sub['zone'] = [zone_at(a, b) for a, b in zip(sub['x'], sub['y'])]
    return pd.crosstab(sub['species'], sub['zone']).reindex(index=SPECIES, columns=ZONES, fill_value=0)


def difference_panel(ax, vs, sp, title, show_ylabel):
    """Recruit share minus adult share per zone with its 95 % interval (direction as found)."""
    ob, sb = vs['obs_boot'][sp], vs['sim_boot'][sp]
    n = min(len(ob), len(sb))
    diff = 100 * (sb[:n] - ob[:n])
    est = 100 * (vs['sim'][sp][0] - vs['obs'][sp][0])
    x = np.arange(len(ZONES))
    for j, z in enumerate(ZONES):
        ax.bar(x[j], est[j], 0.6, facecolor=ZONE_FACE[z], edgecolor=BAR_EDGE, lw=0.7, zorder=2)
        lo, hi = np.percentile(diff[:, j], [2.5, 97.5])
        for col, lw, zo in (('white', 3.2, 5), ('k', 1.3, 6)):
            ax.plot([x[j], x[j]], [lo, hi], color=col, lw=lw, zorder=zo, solid_capstyle='butt')
            for yv in (lo, hi):
                ax.plot([x[j] - 0.08, x[j] + 0.08], [yv, yv], color=col, lw=lw, zorder=zo)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(['C\n(seaward)', 'B', 'A\n(landward)'], fontsize=8)
    ax.set_title(title, fontsize=9.5, style='italic')
    if show_ylabel:
        ax.set_ylabel('Simulated recruits minus observed adults\n(share points of the species\' records)', fontsize=8.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.grid(axis='y', color='0.9', zorder=0)


def figure(vs_sets, withheld, path):
    fig = plt.figure(figsize=(10.5, 16))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.75], hspace=0.75)
    held = ' + '.join(f'T{t}' for t in withheld)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    stacked_panel(axA, vs_sets['heldout'], f'A. Withheld transects {held}: observed adults vs simulated recruits (PCQM operator)',
                  "Share of the species' nearest-individual records (%)")
    stacked_panel(axB, vs_sets['T5'], 'B. Orthogonal transect T5 (used by no reconstruction): observed adults vs simulated recruits',
                  "Share of the species' nearest-individual records (%)")
    for ax in (axA, axB):   # relabel the pair: adults are the observation, recruits the simulation
        for t in ax.texts:
            if t.get_text() == 'Observed':
                t.set_text('Adults\n(observed)')
            elif t.get_text() == 'Simulated':
                t.set_text('Recruits\n(simulated)')
    handles = [Patch(facecolor=ZONE_FACE[z], edgecolor=BAR_EDGE, lw=0.7, label=ZONE_LABEL[z]) for z in ZONES]
    handles.append(plt.Line2D([0], [0], color='k', lw=2.4, label='95 % interval of the zone boundary'))
    axB.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=4, frameon=False, fontsize=9)
    gsC = gs[2].subgridspec(1, len(SPECIES), wspace=0.3)
    axes_c = []
    for i, sp in enumerate(SPECIES):
        ax = fig.add_subplot(gsC[i], sharey=axes_c[0] if axes_c else None)
        difference_panel(ax, vs_sets['heldout'], sp, sp, show_ylabel=(i == 0))
        axes_c.append(ax)
    axes_c[0].set_ylim(-100, 100)
    for ax in axes_c[1:]:
        plt.setp(ax.get_yticklabels(), visible=False)
    pos = axes_c[0].get_position()
    fig.text(0.125, pos.y1 + 0.045, f'C. Withheld transects {held}: recruit share minus adult share per zone (95 % interval)',
             ha='left', fontsize=10.5, fontweight='bold')
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--caches', default=DEFAULT_CACHES)
    ap.add_argument('--stand-caches', default=DEFAULT_STAND_CACHES)
    # 2026-09-22: the withheld transect is T2 (the T1/T3 designs are superseded). T5 is blind in
    # every design and is handled separately inside the script.
    ap.add_argument('--withheld', type=int, nargs='+', default=[2])
    ap.add_argument('--n-boot', type=int, default=4000)
    ap.add_argument('--out-dir', default=os.path.join(SCRIPT_DIR, 'results', 'evaluation'))
    ap.add_argument('--fig', default=os.path.join(SCRIPT_DIR, 'results', 'figures', 'fig_recruits_vs_adults.png'))
    args = ap.parse_args()
    withheld = sorted(args.withheld)
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.fig), exist_ok=True)
    paths = sorted(glob.glob(args.caches)) + sorted(glob.glob(args.stand_caches))
    assert paths, f'no run caches match {args.caches}'
    rng = np.random.default_rng(42)
    frames = G.transect_frames(C.STAND_TRANSECTS + [5])
    real = load_real(withheld + [5], frames)
    runs = [(os.path.basename(p).replace('.csv', ''), load_recruits(p)) for p in paths]
    print(f'{len(runs)} replicate runs; search limit R = {C.DMAX} m'
          f'{" (assumed)" if getattr(C, "SEARCH_LIMIT_ASSUMED", False) else ""}; '
          f'real records in domain: {real.groupby("transect").size().to_dict()}')
    # 2026-09-24 (user): a POOLED set, the withheld transect(s) and T5 together, one observed-adults vs
    # modelled-recruits comparison over all blind plots; the two separate sets are kept for the supplement
    sets = {'heldout': withheld, 'T5': [5], 'pooled': withheld + [5]}
    rows, zone_rows, belts, vs_sets = [], [], [], {}
    for name, ts in sets.items():
        r, zr, _, _, vs = evaluate_set(f'{name}:recruits', ts, real, runs, frames, rng, args.n_boot, margin=100.0)
        for row in r:
            row['operator'] = 'PCQM nearest-individual records'
            row['comparison'] = 'simulated recruits vs observed adults'
            row.pop('equivalent_at_margin', None)
            row.pop('margin_points', None)
            # two-sided 95 % interval of the TVD and the per-zone recruit-minus-adult difference
            i = SPECIES.index(row['species'])
            j = ZONES.index(row['zone'])
            ob, sb = vs['obs_boot'][row['species']], vs['sim_boot'][row['species']]
            n = min(len(ob), len(sb))
            tvd = 0.5 * np.abs(ob[:n] - sb[:n]).sum(axis=1)
            d = sb[:n, j] - ob[:n, j]
            row.update(tvd_lo=float(np.percentile(tvd, 2.5)), tvd_hi=float(np.percentile(tvd, 97.5)),
                       diff_recruit_minus_adult=float(vs['sim'][row['species']][0][j] - vs['obs'][row['species']][0][j]),
                       diff_lo=float(np.percentile(d, 2.5)), diff_hi=float(np.percentile(d, 97.5)))
        rows += r
        zone_rows += zr
        vs_sets[name] = vs
        for label, rec in runs:
            bc = belt_counts(rec, ts)
            for sp in SPECIES:
                belts.append(dict(set=name, run=label, species=sp, **{f'n_{z}': int(bc.loc[sp, z]) for z in ZONES}))
    tests = pd.DataFrame(rows)
    tests.to_csv(os.path.join(args.out_dir, 'recruits_vs_adults_tests.csv'), index=False)
    pd.DataFrame(zone_rows).to_csv(os.path.join(args.out_dir, 'recruits_vs_adults_zone_composition.csv'), index=False)
    pd.DataFrame(belts).to_csv(os.path.join(args.out_dir, 'recruits_vs_adults_belt_counts.csv'), index=False)
    figure({k: v for k, v in vs_sets.items() if k != 'pooled'}, withheld, args.fig)
    pd.set_option('display.width', 220)
    print('\n=== share of the species\' nearest-individual records per zone: observed adults [95 % CI] vs simulated recruits [replicate range]; TVD [95 % CI] ===')
    for name in tests['set'].unique():
        t = tests[tests['set'] == name]
        print(f'--- {name}')
        for sp in SPECIES:
            r = t[t['species'] == sp]
            if not len(r):
                continue
            cells = '  '.join(f"{z}: {100*x.obs_share:4.0f} [{100*x.obs_lo:3.0f}-{100*x.obs_hi:3.0f}] vs {100*x.sim_share:4.0f} [{100*x.sim_lo:3.0f}-{100*x.sim_hi:3.0f}] (diff {100*x.diff_recruit_minus_adult:+4.0f} [{100*x.diff_lo:+3.0f},{100*x.diff_hi:+3.0f}])"
                              for z, x in zip(ZONES, r.itertuples()))
            r0 = r.iloc[0]
            print(f"  {sp:<11} adults n={r0.n_obs_species:3d}  {cells}  | TVD {100*r0.tvd:4.1f} [{100*r0.tvd_lo:4.1f}, {100*r0.tvd_hi:4.1f}]; JSD {r0.jsd_species:.3f}")
    print(f'\nSaved {args.fig}')


if __name__ == '__main__':
    main()
