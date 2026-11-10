"""How the vegetation-drag profile changes simulated zonation (2026-09-23/24; replaces the categorical bar figure
fig_drag_sensitivity.png, which the user rejected: 'I was expecting how increasing or decreasing drag would
affect zonation').

Main figure (user's specification, 2026-09-24): six panels, rows = boundary-origin recruits / parent-tree-origin
recruits (both tracked inside the same full-model runs, zone shares of the established recruits over the WHOLE
domain), columns = the calibrated species; one 100 % stacked bar per drag profile (zones C / B / A in the zone
palette), x = 'Drag profile' with named scenarios Strong -> Calibrated -> Intermediate -> Weak -> Uniform and,
after a gap, Reversed; the calibrated bar outlined; shared y 'Share of recruits (%)'; no observed reference
lines (the observed-simulated comparison belongs to the evaluation figures). The C / B / A retention factors
of each profile go in the caption (table printed to the console and written to
results/drag_sensitivity/drag_profiles_table.csv; smaller retention = stronger drag).
Supplement (--supplement): the same bars for all recruits, whole domain and inside the survey belts (1999
operator), the sampling-area change that was removed from the main figure; drawn for ALL FOUR modelled species
(Sonneratia added 2026-09-28 at the user's request; the main figure stays with the calibrated species).

Mean over seeds 42/49/56 (seed range < 1 point). Caches: results/drag_sensitivity (the calibrated profile =
the pipeline baseline cascade_D1T1S1E1, same seeds).

    python figure_drag_zonation.py                 # compute results/drag_sensitivity/drag_zonation_shares.csv + draw
    python figure_drag_zonation.py --from-csv      # redraw
    python figure_drag_zonation.py --from-csv --supplement
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import cascade_metrics as cm  # noqa: E402
import tidal_dispersal_model as tdm  # noqa: E402

DRAG_DIR = os.path.join(SCRIPT_DIR, 'results', 'drag_sensitivity')
CASCADE_DIR = os.path.join(SCRIPT_DIR, 'results', 'knockout_cascade')
SEEDS = (42, 49, 56)
# 2026-09-28 (user: 'we also need sonneratia in this figure', on the supplement): the CSV holds all four
# modelled species; the main figure keeps the calibrated species (user's specification of 2026-09-24), the
# supplement draws all four. Sonneratia's fitted supply gives ~85 propagules per run (a handful inside the
# survey belts), so its seed ribbon is wide; the n per run is printed for the caption.
SPECIES = list(tdm.SPECIES)
CALIB_SPECIES = list(tdm.CALIBRATED_SPECIES)
ZONES = ['C', 'B', 'A']
ZONE_COL = {'C': '#d0e3f3', 'B': '#adb27c', 'A': '#4b4e36'}   # user's zone palette (2026-09-17)
ZONE_LAB = {'C': 'Zone C (seaward)', 'B': 'Zone B', 'A': 'Zone A (landward)'}
# key, short label, (F_C, F_B, F_A), cache pattern ({seed})
PROFILES = [
    ('steep', 'Strong', (0.2, 0.6, 1.0), os.path.join(DRAG_DIR, 'calib_drag_steep_C0.2_B0.6_A1.0_seed{seed}.csv')),
    ('calibrated', 'Calibrated', (0.4, 0.7, 1.0), os.path.join(CASCADE_DIR, 'cascade_D1T1S1E1_seed{seed}.csv')),
    ('shallow', 'Intermediate', (0.6, 0.8, 1.0), os.path.join(DRAG_DIR, 'calib_drag_shallow_1_C0.6_B0.8_A1.0_seed{seed}.csv')),
    ('near flat', 'Weak', (0.8, 0.9, 1.0), os.path.join(DRAG_DIR, 'calib_drag_shallow_2_C0.8_B0.9_A1.0_seed{seed}.csv')),
    ('no gradient', 'No drag', (1.0, 1.0, 1.0), os.path.join(DRAG_DIR, 'calib_drag_flat_C1.0_B1.0_A1.0_seed{seed}.csv')),
]
GAP = 0.0
# scenario labels of the five ordered profiles (user, 2026-09-24): ticks D1-D5, D2 = baseline; the retention
# factors go into a compact key under the panels, not on the axis
SCENARIO_TICKS = ['D1', 'D2\n(baseline)', 'D3', 'D4', 'D5']
SCENARIO_KEY = [['D1', '0.2 / 0.6 / 1.0'], ['D2 — baseline', '0.4 / 0.7 / 1.0'], ['D3', '0.6 / 0.8 / 1.0'],
                ['D4', '0.8 / 0.9 / 1.0'], ['D5 — no drag', '1.0 / 1.0 / 1.0']]


def compute():
    rows = []
    for key, label, (fc, fb, fa), pat in PROFILES:
        for seed in SEEDS:
            path = pat.format(seed=seed)
            if not os.path.exists(path):
                print('  missing', os.path.relpath(path, SCRIPT_DIR))
                continue
            d = pd.read_csv(path, usecols=['species', 'source', 'outcome', 'zone'])
            d = d[d.species.isin(SPECIES)]
            # 2026-09-24 (user: 'We need stranded version'): two stages, recruits = established only,
            # stranded = every grounded propagule (established + dead), i.e. the deposition before the filter
            for stage, e in (('recruits', d[d.outcome == 'established']),
                             ('stranded', d[d.outcome.isin(['established', 'dead'])])):
                for stratum in ('boundary', 'internal', 'all'):
                    s = e if stratum == 'all' else e[e.source == stratum]
                    for sp in SPECIES:
                        z = s.loc[s.species == sp, 'zone']
                        n = len(z)
                        for zone in ZONES:
                            rows.append({'profile': key, 'F_C': fc, 'F_B': fb, 'F_A': fa, 'seed': seed, 'stage': stage,
                                         'operator': 'domain', 'stratum': stratum, 'species': sp, 'zone': zone,
                                         'share': float((z == zone).mean()) if n else np.nan, 'n': n})
            # belt operator (calibration target) from the run summary, both stages
            sm = json.load(open(cm.summary_path_for(path)))['strata']['all']
            for stage, k_sh, k_n in (('recruits', 'sampled_shares', 'n_sampled'),
                                     ('stranded', 'sampled_shares_grounded', 'n_sampled_grounded')):
                for sp in SPECIES:
                    for zone in ZONES:
                        rows.append({'profile': key, 'F_C': fc, 'F_B': fb, 'F_A': fa, 'seed': seed, 'stage': stage,
                                     'operator': 'belt', 'stratum': 'all', 'species': sp, 'zone': zone,
                                     'share': sm[k_sh][sp][zone], 'n': sm[k_n]})
            print('  read', key, seed, flush=True)
    return pd.DataFrame(rows)


def profile_table():
    t = pd.DataFrame([{'profile': lab, 'F_C': c, 'F_B': b, 'F_A': a} for _, lab, (c, b, a), _ in PROFILES])
    return t


def _stacked_panels(tab, panels, out, figsize, species=None):
    """panels = list of (operator, stratum, row label); one row of species per panel spec."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    species = list(species) if species else CALIB_SPECIES
    nrow = len(panels)
    fig, axes = plt.subplots(nrow, len(species), figsize=figsize, sharex=True, sharey=True, squeeze=False)
    x = np.arange(len(PROFILES), dtype=float)
    x[-1] += GAP                      # the reversed bar stands apart
    width = 0.72
    for i, (op, stratum, rlab) in enumerate(panels):
        for j, sp in enumerate(species):
            ax = axes[i, j]
            bottom = np.zeros(len(PROFILES))
            for zone in ZONES:
                vals = []
                for key, _, _, _ in PROFILES:
                    sel = tab[(tab.operator == op) & (tab.stratum == stratum) & (tab.species == sp)
                              & (tab.zone == zone) & (tab.profile == key)]
                    vals.append(100 * sel.share.mean())
                vals = np.array(vals)
                ax.bar(x, vals, width, bottom=bottom, color=ZONE_COL[zone], edgecolor='0.35', lw=0.5)
                bottom += vals
            # outline the calibrated bar
            k = [p[0] for p in PROFILES].index('calibrated')
            ax.bar([x[k]], [100], width, bottom=0, fill=False, edgecolor='black', lw=2.0, zorder=5)
            ax.axvline(x[-1] - GAP / 2 - 0.5, color='0.6', lw=0.8, ls=':')
            if i == 0:
                ax.set_title(sp, fontsize=12, fontweight='bold', style='italic')
            if j == 0:
                ax.set_ylabel(f'{rlab}\nShare of recruits (%)', fontsize=10)
            if i == nrow - 1:
                ax.set_xlabel('Drag profile', fontsize=10)
            ax.set_xticks(x)
            ax.set_xticklabels([p[1] for p in PROFILES], fontsize=8.5, rotation=30, ha='right')
            ax.set_ylim(0, 100)
            ax.set_yticks([0, 25, 50, 75, 100])
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
    handles = [Patch(facecolor=ZONE_COL[z], edgecolor='0.35', lw=0.5, label=ZONE_LAB[z]) for z in ZONES]
    handles.append(Patch(facecolor='none', edgecolor='black', lw=2.0, label='Calibrated profile'))
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=9.5,
               bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved:', out)


STAGE_TEXT = {'recruits': ('recruits', 'Recruit share (%)'), 'stranded': ('propagules', 'Stranded share (%)')}


def _line_panels(tab, panels, out, figsize, stage='recruits', species=None):
    """Line design (user, 2026-09-24: 'previous design was better, I dont like stack bar ... use drag values';
    then 'Remove calibrated profile, reverse profile. And think about a smart way to present x axis'): one line
    per zone (C / B / A share of the species' recruits) across the five ordered profiles; x = movement retained
    in zone C, and under each tick a two-row table C / B (zone A = 1.0 in every profile); mean over seeds,
    ribbon = seed range. No calibrated band, no reversed profile, no observed lines."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import transforms
    line_col = {'C': '#7fa8cf', 'B': '#adb27c', 'A': '#4b4e36'}   # zone C line darkened so it reads on white
    tab = tab[tab.stage == stage] if 'stage' in tab.columns else tab
    noun, ylab = STAGE_TEXT[stage]
    species = list(species) if species else CALIB_SPECIES
    nrow = len(panels)
    fig, axes = plt.subplots(nrow, len(species), figsize=figsize, sharex=True, sharey=True, squeeze=False)
    grad = [p for p in PROFILES if p[0] != 'reversed']
    keys = [p[0] for p in grad]
    xg = np.array([p[2][0] for p in grad])          # movement retained in zone C
    fb = [p[2][1] for p in grad]                     # ... and in zone B
    for i, (op, stratum, rlab) in enumerate(panels):
        for j, sp in enumerate(species):
            ax = axes[i, j]
            for zone in ZONES:
                col, fill = line_col[zone], ZONE_COL[zone]
                sel = tab[(tab.operator == op) & (tab.stratum == stratum) & (tab.species == sp) & (tab.zone == zone)]
                g = sel[sel.profile.isin(keys)].groupby('profile').share
                mean = 100 * g.mean().reindex(keys).values
                lo = 100 * g.min().reindex(keys).values
                hi = 100 * g.max().reindex(keys).values
                if sp in CALIB_SPECIES:      # user, 2026-09-28: no envelope for Sonneratia (tens of propagules per run)
                    ax.fill_between(xg, lo, hi, color=fill, alpha=0.35, lw=0)
                ax.plot(xg, mean, '-', color=col, lw=2.0)
                ax.plot(xg, mean, 'o', color=fill, mec='0.25', mew=0.7, ms=7)
            if i == 0:
                ax.set_title(sp, fontsize=12, fontweight='bold', style='italic')
            if j == 0:
                ax.set_ylabel(rlab.replace('recruits', noun), fontsize=10.5, fontweight='bold')
            ax.set_xticks(xg)
            ax.set_xticklabels([])          # user, 2026-09-24 01:40: no D1-D5 labels, only the direction words
            ax.tick_params(axis='x', length=3)
            if i == nrow - 1:
                # direction of the axis in words (user, 2026-09-24: 'moving right means less drag'):
                # 'Strong drag' at the left end, 'No drag' at the right end, an arrow between them
                ax.text(0.0, -0.08, 'Strong drag', transform=ax.transAxes, ha='left', va='center', fontsize=10)
                ax.text(1.0, -0.08, 'No drag', transform=ax.transAxes, ha='right', va='center', fontsize=10)
                ax.annotate('', xy=(0.70, -0.08), xytext=(0.30, -0.08), xycoords='axes fraction',
                            arrowprops=dict(arrowstyle='-|>', lw=1.4, color='0.25', mutation_scale=14))
            ax.set_xlim(0.1, 1.08)
            ax.set_ylim(0, 100)
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.grid(True, axis='y', lw=0.4, alpha=0.5)
            for s_ in ('top', 'right'):
                ax.spines[s_].set_visible(False)
    handles = [Line2D([], [], color=line_col[z], lw=2, marker='o', mfc=ZONE_COL[z], mec='0.25', mew=0.7, ms=7,
                      label=ZONE_LAB[z]) for z in ZONES]
    # colour legend (zones) only; the scenario key with the retention factors was removed at the user's request
    # (2026-09-24 01:30, 'Remove the profiles and their numbers'); the factors stay in the caption / drag_profiles_table.csv
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, -0.005))
    fig.supylabel(ylab, fontsize=11, x=0.012)
    fig.tight_layout(rect=(0.02, 0.06, 1, 1))
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved:', out)


MAIN_PANELS = [('domain', 'boundary', 'Boundary-origin recruits'),
               ('domain', 'internal', 'Parent-tree-origin recruits')]
SUPPL_PANELS = [('domain', 'all', 'All recruits, whole domain'),
                ('belt', 'all', 'All recruits, survey belts')]


def draw_main(tab, out, bars=False, stage='recruits', species=None):
    """Main layout (boundary-origin / parent-tree-origin rows). All four species since 2026-09-28 (user: the pooled
    supplement 'does not make sense' as the place for Sonneratia, the boundary concentration in zone C under strong
    drag being what should be seen); --calibrated-only restores the three-species version of 2026-09-24."""
    species = list(species) if species else SPECIES
    figsize = (11.5 * len(species) / 3, 7.4)
    if bars:
        _stacked_panels(tab[tab.stage == stage], MAIN_PANELS, out, figsize, species=species)
    else:
        _line_panels(tab, MAIN_PANELS, out, figsize, stage=stage, species=species)


def draw_supplement(tab, out, bars=False, stage='recruits'):
    """All four modelled species, Sonneratia included (user, 2026-09-28); one column per species."""
    figsize = (11.5 * len(SPECIES) / 3, 7.4)
    if bars:
        _stacked_panels(tab[tab.stage == stage], SUPPL_PANELS, out, figsize, species=SPECIES)
    else:
        _line_panels(tab, SUPPL_PANELS, out, figsize, stage=stage, species=SPECIES)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--csv', default=os.path.join(DRAG_DIR, 'drag_zonation_shares.csv'))
    ap.add_argument('--out', default=os.path.join(SCRIPT_DIR, 'results', 'figures', 'fig_drag_zonation.png'))
    ap.add_argument('--from-csv', action='store_true')
    ap.add_argument('--supplement', action='store_true', help='also write the all-recruits figure (whole domain + belts)')
    ap.add_argument('--bars', action='store_true', help='100 %% stacked bars instead of lines (rejected by the user 2026-09-24)')
    ap.add_argument('--stage', default='recruits', choices=['recruits', 'stranded'],
                    help='recruits = established only; stranded = all grounded propagules (established + dead)')
    ap.add_argument('--calibrated-only', action='store_true', help='main figure with the three calibrated species only')
    args = ap.parse_args()
    if args.from_csv:
        tab = pd.read_csv(args.csv)
    else:
        tab = compute()
        tab.to_csv(args.csv, index=False, float_format='%.4f')
        print('Saved:', args.csv)
    if 'stage' not in tab.columns:
        raise SystemExit('CSV predates the stage column: rerun without --from-csv')
    out = args.out if args.stage == 'recruits' else args.out.replace('.png', '_stranded.png')
    draw_main(tab, out, bars=args.bars, stage=args.stage, species=CALIB_SPECIES if args.calibrated_only else SPECIES)
    if args.supplement:
        draw_supplement(tab, out.replace('.png', '_all_SUPPL.png'), bars=args.bars, stage=args.stage)
    pt = profile_table()
    pt.to_csv(os.path.join(DRAG_DIR, 'drag_profiles_table.csv'), index=False)
    print('\nDrag profiles (movement retained in zones C / B / A; smaller = stronger drag):')
    print(pt.to_string(index=False))
    piv = (tab[tab.stage == args.stage].groupby(['operator', 'stratum', 'species', 'zone', 'profile']).share.mean()
           .unstack('profile')[[p[0] for p in PROFILES]])
    piv.columns = [p[1] for p in PROFILES]
    print('\nZone shares, mean over seeds:')
    print((100 * piv).round(0).astype('Int64').to_string())   # blank = no propagule of that stratum (Sonneratia)
    # sample sizes per species and run (whole domain, all sources) for the caption: Sonneratia's fitted supply
    # is tiny, so its shares rest on tens of propagules per run
    nn = (tab[(tab.stage == args.stage) & (tab.operator == 'domain') & (tab.stratum == 'all') & (tab.zone == 'C')]
          .groupby(['species', 'profile']).n.agg(['min', 'max']))
    print(f'\nPropagules per run behind the whole-domain shares ({args.stage}, all sources; min-max over seeds):')
    print(nn.to_string())


if __name__ == '__main__':
    main()
