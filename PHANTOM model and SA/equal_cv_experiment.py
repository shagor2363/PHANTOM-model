"""Equal-contrast test of the three transport traits (secondary, mechanistic analysis).

The knockout sweeps (figure_transport_knockout.py, fig_mechanisms.png) retain a fraction alpha of each
trait's CALIBRATED species contrast. Those contrasts differ in size (population CV across the modelled
species; with the three calibrated species: draft 0.55, mobility 0.13, lateral spread 0.34), so
'alpha = 0' removes very different amounts of relative variation and the endpoint effects (-73 % draft,
-36 % mobility) rank contributions of the calibrated contrasts, not the mechanisms. Here every trait is
given the SAME relative contrast c (target CV), keeping the calibrated mean and the calibrated pattern
of species differences:

    T_i(c) = mean(T) + (c / CV_cal) * (T_i,cal - mean(T))

one trait at a time, the other two at their calibrated values. Levels: c = 0 (all species equal; the
existing knockout_*_eq_mean caches), CV_mob/3, 2 CV_mob/3, CV_mob (the largest contrast within every
trait's calibrated range: mobility at this level IS the baseline), then the extrapolated levels up to
the ceiling CV_MAX (a trait at its own CV_cal is the baseline; the other traits are then outside their
calibrated ranges: to be labelled as extrapolation). Seeds are paired: the same seed set runs every
trait x level. Stand and constants as the final runs.

2026-09-20: 4 species (Sonneratia). Target = mean of the calibrated species; Sonneratia (uncalibrated,
set by rule) blends toward the same target but does not define it: CV_cal per trait, the mean and the
ceiling (CV_ZERO, CV_MAX: where the smallest calibrated value reaches zero, Avicennia's draft) are
computed over CALIBRATED_SPECIES only, and T_i(c) = mean_cal + (c / CV_cal) (T_i,cal - mean_cal) is
applied to all four species, Sonneratia's values floored at 0. Sonneratia's draft (0.01 m, below the
calibrated mean) therefore saturates at 0 above c = CV_SAT['draft'] (about 0.61; its mobility and
lateral spread equal Avicennia's and never saturate). All numbers are derived from the constants at
import time (printed); the level grid is that of the 3-species design.

Usage:
  python equal_cv_experiment.py run [--seeds 42,49,56] [--workers 10] [--levels common|all]
  python equal_cv_experiment.py figure [--seeds 42,49,56]
Outputs: results/equal_cv/knockout_eqcv_{trait}_c{ccc}_seed{seed}.csv, results/equal_cv/equal_cv_w1.csv,
         results/figures/fig_equal_cv.png
"""
import argparse
import os
import sys
from itertools import combinations

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
STAND = os.path.join(SCRIPT_DIR, 'results', 'reconstruction', 'stand3T_T4T3T1T6_single', 'demo_reconstructed_stand3T.csv')
os.environ.setdefault('KNOCKOUT_TREES', STAND)
import transport_treatments as tt  # noqa: E402  (reads the model constants at import)

RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results', 'equal_cv')
KNOCKOUT_DIR = os.path.join(SCRIPT_DIR, 'results', 'knockout_final')
SPECIES = tt.SPECIES_ORDER   # 2026-09-20: 4 species (Sonneratia): the vector Eq. S1 is applied to
CALIBRATED = tt.CALIBRATED_SPECIES   # the species that define the target mean, CV_cal and the ceiling

TRAITS = {  # trait -> (tdm attribute, calibrated values, label)
    'draft': ('DRAFT', tt.DRAFT_VALS, 'Draft'),
    'mobility': ('MOBILITY', tt.MOBILITY_VALS, 'Mobility'),
    'spread': ('SIGMA_LATERAL', tt.SIGMA_LAT_VALS, 'Lateral spread'),
}


def cal_mean(vals):
    """Target mean of a trait: over the calibrated species only (2026-09-20)."""
    return float(np.mean([vals[sp] for sp in CALIBRATED]))


def cv(vals):
    """Population CV of a trait over the CALIBRATED species (2026-09-20: Sonneratia excluded)."""
    a = np.array([vals[sp] for sp in CALIBRATED], dtype=float)
    return float(a.std() / a.mean())


CV_CAL = {t: cv(v[1]) for t, v in TRAITS.items()}
CV_MOB = min(CV_CAL.values())   # the smallest calibrated contrast: the largest level inside EVERY trait's range
CV_MOB_TRAIT = min(CV_CAL, key=CV_CAL.get)   # which trait that is (mobility)
CV_DRAFT = CV_CAL['draft']


def cv_at_zero(trait, species=None):
    """The target CV at which the smallest value of `trait` among `species` (default: the calibrated
    species) reaches zero under Eq. S1: mean_cal + (c / CV_cal) (T_min - mean_cal) = 0  ->
    c = CV_cal * mean_cal / (mean_cal - T_min). A trait value cannot be negative, so over the calibrated
    species this is the trait's own ceiling (2026-09-20: computed, no longer typed)."""
    _, vals, _ = TRAITS[trait]
    species = CALIBRATED if species is None else species
    mean = cal_mean(vals)
    t_min = float(min(vals[sp] for sp in species))
    return float('inf') if mean - t_min <= 0 else CV_CAL[trait] * mean / (mean - t_min)


CV_ZERO = {t: cv_at_zero(t) for t in TRAITS}            # over the calibrated species (Avicennia's draft)
CV_ZERO_TRAIT = min(CV_ZERO, key=CV_ZERO.get)
# Sonneratia (uncalibrated, set by rule) blends toward the calibrated target but does not define the
# ceiling: where one of its values would turn negative it is floored at 0 (saturation level per trait).
UNCALIBRATED = [sp for sp in SPECIES if sp not in CALIBRATED]
CV_SAT = {t: cv_at_zero(t, UNCALIBRATED) for t in TRAITS} if UNCALIBRATED else {}
# The typed grid of the 3-species design (pilot of 2026-09-14, 39 runs: levels 0, CV_mob/3, 2CV_mob/3,
# CV_mob, 0.30, CV_draft; densified grid, user's choice: 0.16, 0.20, 0.25, 0.35, 0.40, 0.45, 0.50 (63 runs);
# stress-test extension, user: draft must be pushed beyond its calibrated contrast too, otherwise the
# design is asymmetric: 0.60, 0.65, 0.70, 0.75 (36 runs); 0.75 is the ceiling because Avicennia's draft
# reaches zero at CV 0.76, mean kept at 0.108 m).
_GRID = [0.25, 0.45, 0.65, 0.75]   # 2026-09-23: reduced from 12 to 4 extrapolated levels (plus the CV_cal anchors) to halve runtime; the full grid is in equal_cv_experiment.py.bak_fullgrid
# Ceiling = the last grid level below the CV at which a calibrated value reaches zero (0.75 with the
# calibrated constants; 2026-09-20: computed from the constants, over the calibrated species only).
CV_MAX = max(c for c in _GRID if c < min(CV_ZERO.values()) - 1e-9)
COMMON_LEVELS = [0.0, CV_MOB / 3, 2 * CV_MOB / 3, CV_MOB]     # within every calibrated range (0 .. CV_mob)
# every trait's own CV_cal is a level (it reuses the baseline cache); the grid is cut at the ceiling
EXTRA_LEVELS = sorted({round(c, 6) for c in _GRID + list(CV_CAL.values())
                       if CV_MOB + 1e-9 < c <= CV_MAX + 1e-9})
print(f'equal_cv_experiment: {len(SPECIES)} species {SPECIES}, target/CV over {CALIBRATED}; CV_cal '
      + ', '.join(f'{t} {v:.3f}' for t, v in CV_CAL.items())
      + '; calibrated value reaches zero at c = '
      + ', '.join(f'{t} {v:.3f}' for t, v in CV_ZERO.items())
      + f' -> CV_MAX = {CV_MAX:.2f} (set by {CV_ZERO_TRAIT})'
      + ''.join(f"; NOTE {sp}'s {t} saturates at 0 above c = {CV_SAT[t]:.2f}"
                for sp in UNCALIBRATED for t in TRAITS if CV_SAT[t] <= CV_MAX + 1e-9))


def level_tag(c):
    return f'c{int(round(c * 1000)):03d}'


def scaled_values(trait, c):
    """Eq. S1 for every modelled species around the calibrated-species mean; a complete dict (copy of
    the calibrated one) so the override carries every species key the model expects. 2026-09-20:
    target = mean of the calibrated species; Sonneratia (uncalibrated, set by rule) blends toward the
    same target but does not define it, and its values are floored at 0 (saturation, CV_SAT)."""
    attr, vals, _ = TRAITS[trait]
    mean = cal_mean(vals)
    k = c / CV_CAL[trait]
    out = dict(vals)
    out.update({sp: mean + k * (vals[sp] - mean) for sp in CALIBRATED})
    out.update({sp: max(0.0, mean + k * (vals[sp] - mean)) for sp in UNCALIBRATED})
    return out


def treatment_key(trait, c):
    """Cache key; reuses existing knockout caches where the scaled values equal them exactly."""
    if c == 0.0:
        return f'{trait}_eq_mean', KNOCKOUT_DIR
    if abs(c - CV_CAL[trait]) < 1e-9:
        return 'baseline', KNOCKOUT_DIR
    return f'eqcv_{trait}_{level_tag(c)}', RESULTS_DIR


def register_treatments(levels):
    plan = []
    for trait in TRAITS:
        for c in levels:
            key, d = treatment_key(trait, c)
            if d == RESULTS_DIR:
                attr = TRAITS[trait][0]
                tt.TREATMENTS.setdefault(key, {attr: scaled_values(trait, c)})
            plan.append((trait, c, key, d))
    return plan


def cmd_run(args):
    seeds = [int(s) for s in args.seeds.split(',')]
    levels = COMMON_LEVELS + (EXTRA_LEVELS if args.levels == 'all' else [])
    os.makedirs(RESULTS_DIR, exist_ok=True)
    plan = register_treatments(levels)
    print('stand:', os.environ['KNOCKOUT_TREES'])
    print('calibrated CV:', {t: round(v, 3) for t, v in CV_CAL.items()})
    for trait, c, key, d in plan:
        vals = scaled_values(trait, c)
        print(f'  {trait:9s} c={c:.4f} -> {key:26s} ', {tt.SPECIES_SHORT[s]: round(vals[s], 3) for s in SPECIES},
              '(existing cache)' if d == KNOCKOUT_DIR else '')
    # only treatments that still lack a cache for at least one seed: run_treatments would otherwise load
    # every existing 45 MB cache into memory before submitting the pending runs
    new_keys = sorted({key for _, _, key, d in plan if d == RESULTS_DIR
                       and not all(os.path.exists(os.path.join(d, f'knockout_{key}_seed{s}.csv')) for s in seeds)})
    print(f'{len(new_keys)} treatments without a complete cache x {len(seeds)} seeds')
    tt.run_treatments(new_keys, seeds, args.days, RESULTS_DIR, max_workers=args.workers)
    print('run complete')


def pairwise_w1(df, species=None):
    """Mean pairwise W1 of terminal elevation over the pairs of `species` (default: all modelled
    species, 6 pairs with four; pairs without data are left out). 2026-09-20: pass
    tt.CALIBRATED_SPECIES for `w1_calibrated`, the 3-pair mean comparable with the earlier runs."""
    from scipy.stats import wasserstein_distance
    species = SPECIES if species is None else species
    zs = {sp: df.loc[df.species == sp, 'terminal_z'].dropna().values for sp in species}
    vals = [wasserstein_distance(zs[a], zs[b]) for a, b in combinations(species, 2)
            if len(zs[a]) > 1 and len(zs[b]) > 1]
    return float(np.mean(vals)) if vals else np.nan


def cmd_figure(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    seeds = [int(s) for s in args.seeds.split(',')]
    metric = args.metric   # 2026-09-23: 'w1_calibrated' (default) or 'w1' (all four species)
    levels = COMMON_LEVELS + EXTRA_LEVELS
    plan = register_treatments(levels)
    rows = []
    cache = {}
    for trait, c, key, d in plan:
        for seed in seeds:
            path = os.path.join(d, f'knockout_{key}_seed{seed}.csv')
            if not os.path.exists(path):
                print('  missing', os.path.relpath(path, SCRIPT_DIR))
                continue
            if path not in cache:
                df = pd.read_csv(path, usecols=['species', 'source', 'outcome', 'terminal_z'])
                df = df[df.outcome.isin(['established', 'dead'])]
                strata_df = {'boundary': df[df.source == 'boundary'], 'internal': df[df.source == 'internal'],
                             'all': df}
                # 2026-09-20: both metrics per run, all pairs (w1) and calibrated pairs (w1_calibrated)
                cache[path] = {s: (pairwise_w1(d), pairwise_w1(d, tt.CALIBRATED_SPECIES))
                               for s, d in strata_df.items()}
            for stratum, (w, w_cal) in cache[path].items():
                rows.append({'trait': trait, 'target_cv': c, 'treatment': key, 'seed': seed,
                             'stratum': stratum, 'w1': w, 'w1_calibrated': w_cal,
                             'extrapolated': c > CV_CAL[trait] + 1e-9})
    tab = pd.DataFrame(rows)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_csv = os.path.join(RESULTS_DIR, 'equal_cv_w1.csv')
    tab.to_csv(out_csv, index=False, float_format='%.5f')
    print('Saved:', out_csv)

    colors = {'draft': '#2A9D8F', 'mobility': '#3B4CC0', 'spread': '#E69F00'}
    markers = {'draft': 'o', 'mobility': 's', 'spread': '^'}
    strata = [('boundary', 'A. Common boundary source'), ('internal', 'B. Species-specific parent sources'),
              ('all', 'C. Full source configuration')]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True, sharey=True)
    for ax, (stratum, title) in zip(axes, strata):
        # one line per trait (the panels are per trait, not per species): the status of the x range
        # differs by trait, since each trait leaves its calibrated contrast at its own CV_cal
        for trait, (_, _, label) in TRAITS.items():
            t = tab[(tab.trait == trait) & (tab.stratum == stratum)]
            # increment above the trait's own zero-contrast run (other two traits calibrated),
            # per seed, so the paired design carries through
            t = t.copy()
            base = t[t.target_cv == 0.0].set_index('seed')[metric]
            t['dw1'] = t[metric].values - base.reindex(t.seed).values
            g = t.groupby('target_cv').dw1
            x = np.array(sorted(g.groups))
            mean = g.mean().loc[x].values
            lo = g.min().loc[x].values
            hi = g.max().loc[x].values
            inside = x <= CV_CAL[trait] + 1e-9
            ax.fill_between(x, lo, hi, color=colors[trait], alpha=0.15, lw=0)
            ax.plot(x, mean, "-", color=colors[trait], lw=1.8)
            ax.plot(x, mean, ls='none', marker=markers[trait], color=colors[trait], ms=5.5, label=label)
            # 2026-09-23 (user: the alpha sweep says draft matters most, this figure says mobility): mark
            # the trait's OWN calibrated contrast on its curve, so the reader sees both readings at once,
            # the slope (effect per unit relative contrast) and the value at the contrast the species
            # actually have (= the endpoint drop of the alpha sweep, same runs, same metric)
            # 2026-09-23 22:40: removed again at the user's request ('Remove the calibrated contrast of
            # the trait'); kept behind --mark-calibrated, the reconciliation goes in the text instead
            if args.mark_calibrated:
                i_cal = int(np.argmin(abs(x - CV_CAL[trait])))
                ax.plot(x[i_cal], mean[i_cal], ls='none', marker=markers[trait], ms=13, mfc='white',
                        mec=colors[trait], mew=2.2, zorder=5)
        ax.set_title(title, fontsize=11, fontweight='bold', loc='left')
        ax.set_xlabel('Species specific difference in trait', fontsize=10)
        ticks = np.arange(0, CV_MAX + 1e-9, 0.10)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=8.5)
        ax.axhline(0, color='0.5', lw=0.6)
        ax.set_xlim(-0.01, CV_MAX + 0.02)
        ax.grid(True, axis='y', lw=0.4, alpha=0.5)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel('Elevation difference between species (m)', fontsize=10)
    h, l = axes[0].get_legend_handles_labels()
    if args.mark_calibrated:
        from matplotlib.lines import Line2D
        h.append(Line2D([], [], ls='none', marker='o', ms=11, mfc='white', mec='0.3', mew=2.0))
        l.append('Calibrated contrast of the trait (as in the model)')
    fig.legend(h, l, loc='lower center', ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    out = args.out or os.path.join(SCRIPT_DIR, 'results', 'figures',
                                   'fig_equal_cv.png' if metric == 'w1_calibrated' else 'fig_equal_cv_4species.png')
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved:', out)

    # increment above the trait's own zero-contrast run, per paired seed, into the CSV and the console
    idx = list(zip(tab.stratum, tab.trait, tab.seed))
    base = tab[tab.target_cv == 0.0].set_index(['stratum', 'trait', 'seed'])
    tab['dw1'] = tab.w1.values - base.w1.reindex(idx).values
    tab['dw1_calibrated'] = tab.w1_calibrated.values - base.w1_calibrated.reindex(idx).values   # 2026-09-20
    tab.to_csv(out_csv, index=False, float_format='%.5f')
    piv = tab.groupby(['stratum', 'trait', 'target_cv']).dw1.mean().unstack('target_cv')
    print('\nSeparation added above zero contrast (m), mean over paired seeds:')
    print(piv.round(4).to_string())
    top = tab[(abs(tab.target_cv - CV_MOB) < 1e-9)].groupby(['stratum', 'trait']).dw1.agg(['mean', 'min', 'max'])
    print(f'\nAt matched CV = {CV_MOB:.3f}:')
    print(top.round(4).to_string())
    piv_c = tab.groupby(['stratum', 'trait', 'target_cv']).dw1_calibrated.mean().unstack('target_cv')
    print('\nSame, calibrated pairs only (w1_calibrated, comparable with the 3-species runs):')
    print(piv_c.round(4).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--seeds', default='42,49,56')
    r.add_argument('--workers', type=int, default=10)
    r.add_argument('--days', type=int, default=365)
    r.add_argument('--levels', choices=['common', 'all'], default='all')
    r.set_defaults(func=cmd_run)
    f = sub.add_parser('figure')
    f.add_argument('--seeds', default='42,49,56')
    # 2026-09-23: the mechanism figure is drawn on the CALIBRATED species (the pairs the calibration
    # targets); 'w1' adds the Sonneratia pairs, whose assumed draft (0.02 m) is multiplied 5x when
    # equalised and dominates the parent-tree panel (sign flip); that version is a supplement
    f.add_argument('--metric', default='w1_calibrated', choices=['w1', 'w1_calibrated'])
    f.add_argument('--out', default=None)
    f.add_argument('--mark-calibrated', action='store_true',
                   help="hollow marker at each trait's calibrated contrast (off by default, user 2026-09-23)")
    f.set_defaults(func=cmd_figure)
    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
