"""Evaluation figure, user's final three-panel specification (2026-09-24 08:10): A. the whole 1999 juvenile
survey observed vs simulated (calibration year), B. change since the 1997 survey for the observed 1999 survey and for
the model, C. observed adults on the withheld transects (T2, T5 pooled) vs modelled recruits.

Derived from figure_evaluation_three_patterns_v2.py: (2026-09-24; user: 'Panel C does not seem
to be that convincing, but I think it is a matter of presentation').

Panels A and B are unchanged (figure_evaluation_three_patterns.py). Panel C, withheld ADULT zonation vs modelled
RECRUIT zonation on the withheld transect (T2) and on T5 (never used by any reconstruction), now uses the grammar of
panel A and shows the observed uncertainty, which the bar version lacked:

  C1  per species and transect ONE PAIR of stacked composition bars, observed adults (hatched) | modelled recruits,
      zones C / B / A stacked, with the 95 % interval of each zone boundary on BOTH bars (observed: multinomial
      bootstrap of the n nearest-individual records; model: replicate range), n under each pair; where the model
      places no recruits at the plots (Sonneratia on T2) the bar is empty and says so.
  C2  parity plot of the same numbers, x = observed adult share of a zone, y = modelled recruit share, one point
      per species x transect x zone (fill = zone colour, marker = transect), horizontal whisker = observed 95 %
      interval, vertical = model range, 1:1 line, grey band = the a priori +-15-point margin of the 1999
      goodness-of-fit figure. Structural zeros (zone C on T5, which has no zone C plots) are left out.

Both sides of panel C use the same nearest-individual (PCQM) operator at the same plots. It is a comparison across
life stages, not a prediction of the adult stand.

Usage: python figure_evaluation_three_patterns_v2.py [--tests results/evaluation/recruits_vs_adults_tests.csv]
           [--caches 'results/knockout_cascade/cascade_D1T1S1E1_seed*.csv'] [--out ...]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

TRANSECT_MARKER = {'T2': 'o', 'T5': 's'}


def bracket(ax, xb, lo, hi, lw=2.4):
    if hi - lo < 1e-9:
        return
    ax.plot([xb, xb], [lo, hi], color='k', lw=lw, zorder=6, solid_capstyle='butt', clip_on=False)
    for yv in (lo, hi):
        ax.plot([xb - 0.06, xb + 0.06], [yv] * 2, color='k', lw=lw, zorder=6, clip_on=False)


def panel_c1(ax, tests, sets, set_label, SPECIES, ZONES, F, rng, n_boot=4000, pooled=False):
    w, gap, pair_gap, sp_gap = (0.40, 0.03, 0.0, 0.47) if pooled else (0.34, 0.03, 0.22, 0.55)
    pair_w = 2 * w + gap
    xs_species = []
    x = 0.0
    for i, sp in enumerate(SPECIES):
        x_start = x
        for k, st in enumerate(sets):
            d = tests[(tests['set'] == st) & (tests.species == sp)].set_index('zone').reindex(ZONES)
            n = int(d.n_obs_species.iloc[0]) if len(d) and not np.isnan(d.n_obs_species.iloc[0]) else 0
            counts = d.n_obs_adults.fillna(0).values.astype(int)
            x_obs, x_sim = x + w / 2, x + w / 2 + w + gap
            # observed adults: stacked shares, multinomial bootstrap of the records for the boundary brackets
            if n > 0:
                p = counts / counts.sum()
                bottom = 0.0
                for j, z in enumerate(ZONES):
                    ax.bar(x_obs, 100 * p[j], w, bottom=100 * bottom, facecolor=F.ZONE_FACE[z],
                           edgecolor='white' if z == 'A' else F.BAR_EDGE, lw=0.7, hatch='///', zorder=2)
                    bottom += p[j]
                boot = rng.multinomial(n, p, size=n_boot) / n
                cum = np.cumsum(boot, axis=1)[:, :2]
                for b, xb in enumerate((x_obs - w / 4, x_obs + w / 4)):
                    lo, hi = np.percentile(cum[:, b], [2.5, 97.5])
                    bracket(ax, xb, 100 * lo, 100 * hi)
            # modelled recruits: stacked mean shares, replicate range as the brackets
            if d.sim_share.notna().all() and d.sim_share.sum() > 0:
                s = d.sim_share.values
                s = s / s.sum()          # mean replicate shares need not sum exactly to 1: renormalise for the stack
                bottom = 0.0
                for j, z in enumerate(ZONES):
                    ax.bar(x_sim, 100 * s[j], w, bottom=100 * bottom, facecolor=F.ZONE_FACE[z], edgecolor=F.BAR_EDGE,
                           lw=0.7, zorder=2)
                    bottom += s[j]
                lo_c, hi_c = d.sim_lo['C'], d.sim_hi['C']
                lo_a, hi_a = d.sim_lo['A'], d.sim_hi['A']
                bracket(ax, x_sim - w / 4, 100 * lo_c, 100 * hi_c)
                bracket(ax, x_sim + w / 4, 100 * (1 - hi_a), 100 * (1 - lo_a))
            else:
                ax.bar(x_sim, 100, w, facecolor='white', edgecolor=F.BAR_EDGE, lw=0.7, ls=(0, (2, 2)), zorder=2)
                ax.text(x_sim, 50, 'no recruits\nat these plots', ha='center', va='center', fontsize=7.2,
                        color='#555555', rotation=90)
            if pooled:
                ax.text(x_obs, -4.5, 'Observed' + chr(10) + f'n = {n}', ha='center', va='top', fontsize=8.5, color='#333333', linespacing=1.2)
                ax.text(x_sim, -4.5, 'Simulated', ha='center', va='top', fontsize=8.5, color='#333333')
            else:
                ax.text(x + pair_w / 2, -4.5, f'{set_label[st]}' + chr(10) + f'n = {n}', ha='center', va='top', fontsize=8.2,
                        color='#333333', linespacing=1.2)
            x += pair_w + pair_gap
        xs_species.append((x_start + x - pair_gap) / 2)
        x += sp_gap
    ax.set_xticks(xs_species)
    ax.set_xticklabels([f'$\\it{{{sp}}}$' for sp in SPECIES], fontsize=10.5)
    ax.tick_params(axis='x', pad=30, length=0)
    ax.set_xlim(-0.25, x - sp_gap + 0.1)
    ax.set_ylim(-2.5, 101)
    ax.spines['bottom'].set_position(('data', -2.5))
    ax.set_ylabel("Share of the species' nearest-individual\nrecords per zone (%)", fontsize=9.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def panel_c2(ax, tests, sets, set_label, SPECIES, ZONES, F, margin=15.0):
    xx = np.array([0, 100])
    ax.fill_between(xx, xx - margin, xx + margin, color='#ececec', lw=0, zorder=0)
    ax.plot(xx, xx, color='k', lw=0.9, zorder=1)
    pts = []
    for st in sets:
        for sp in SPECIES:
            d = tests[(tests['set'] == st) & (tests.species == sp)].set_index('zone').reindex(ZONES)
            if d.sim_share.isna().any() or d.sim_share.sum() == 0:
                continue
            for z in ZONES:
                o, s = 100 * d.obs_share[z], 100 * d.sim_share[z]
                if set_label[st] == 'T5' and z == 'C':
                    continue          # structural zero: T5 has no zone C plots in the model frame
                ol, oh = 100 * d.obs_lo[z], 100 * d.obs_hi[z]
                sl, sh = 100 * d.sim_lo[z], 100 * d.sim_hi[z]
                ax.plot([ol, oh], [s, s], color='0.65', lw=0.8, zorder=2)
                ax.plot([o, o], [sl, sh], color='0.65', lw=0.8, zorder=2)
                ax.plot(o, s, marker=TRANSECT_MARKER[set_label[st]], ms=8, mfc=F.ZONE_FACE[z], mec='k', mew=0.8,
                        ls='none', zorder=4)
                pts.append((sp, set_label[st], z, o, s))
    ax.set_xlim(-3, 103)
    ax.set_ylim(-3, 103)
    ax.set_aspect('equal')
    ax.set_xlabel('Observed adult share (%)', fontsize=9.5)
    ax.set_ylabel('Modelled recruit share (%)', fontsize=9.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    diff = np.array([s - o for _, _, _, o, s in pts])
    return pts, diff


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tests', default=os.path.join(SCRIPT_DIR, 'results', 'evaluation', 'recruits_vs_adults_tests.csv'))
    ap.add_argument('--caches', default=os.path.join(SCRIPT_DIR, 'results', 'knockout_cascade', 'cascade_D1T1S1E1_seed*.csv'))
    ap.add_argument('--survey-out-dir', default=os.path.join(SCRIPT_DIR, 'results', 'evaluation'))
    ap.add_argument('--n-boot', type=int, default=4000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--margin', type=float, default=15.0)
    ap.add_argument('--out', default=os.path.join(SCRIPT_DIR, 'results', 'figures', 'fig_evaluation_1999_change_adults.png'))
    args = ap.parse_args()

    import virtual_survey as v
    v.CACHE_GLOB = os.path.abspath(args.caches)
    v.OUT_DIR = os.path.abspath(args.survey_out_dir)
    from virtual_survey import virtual_survey, SPECIES, ZONES  # noqa: E402
    import figure_evaluation_stacked as F  # noqa: E402

    vs99 = virtual_survey('1999', args.n_boot, args.seed)
    vs97 = virtual_survey('1997', args.n_boot, args.seed)
    tests = pd.read_csv(args.tests)
    rng = np.random.default_rng(args.seed)

    fig = plt.figure(figsize=(6.2 + 2.2 * len(SPECIES), 14.6))
    gs = fig.add_gridspec(3, 6, height_ratios=[1.15, 0.85, 1.05], hspace=0.72, wspace=0.10, top=0.955, bottom=0.05)

    # --- A: 1997 withheld juvenile composition (unchanged)
    axA = fig.add_subplot(gs[0, 1:5])
    # 2026-09-24 11:05 (user: 'Panel A is against 1997 juveniles not 1999'): the withheld 1997 survey, all quadrats
    F.stacked_panel(axA, vs97, 'A. 1997 withheld juvenile survey (all quadrats): observed vs simulated zonation')
    h = [Patch(facecolor=F.ZONE_FACE[z], edgecolor=F.BAR_EDGE, lw=0.7) for z in ZONES] + [plt.Line2D([], [], color='k', lw=2.4)]
    axA.legend(h, [F.ZONE_LABEL[z] for z in ZONES] + ['95 % CI of the zone boundary'], loc='upper center',
               bbox_to_anchor=(0.5, -0.17), ncol=4, fontsize=9, frameon=False, columnspacing=2.0)

    # --- B: change since 1997 (unchanged)
    subB = gs[1, :].subgridspec(1, len(SPECIES), wspace=0.12)
    axesB = [fig.add_subplot(subB[0, i]) for i in range(len(SPECIES))]
    for i, sp in enumerate(SPECIES):
        F.change_panel(axesB[i], sp, vs97, vs99, show_ylabel=(i == 0))
        axesB[i].set_title(f'$\\it{{{sp}}}$', fontsize=10.5, pad=4)
        if i > 0 and np.allclose(axesB[i].get_ylim(), axesB[0].get_ylim()):
            axesB[i].tick_params(axis='y', labelleft=False)
    bb = axesB[0].get_position()
    fig.text(bb.x0, bb.y1 + 0.028, 'B. Change since the 1997 survey: observed 1999 survey vs model',
             fontsize=10.5, fontweight='bold', ha='left', va='bottom')
    hB = [Patch(facecolor=F.CHG_LEGEND_FACE, edgecolor=F.BAR_EDGE, lw=0.7, **F.CHG_STYLE[k]) for k in ('obs', 'model')]
    hB.append(plt.Line2D([], [], color='k', lw=1.3))
    b0, b1 = axesB[0].get_position(), axesB[-1].get_position()
    fig.legend(hB, [F.CHG_LABEL['obs'], F.CHG_LABEL['model'], '95 % CI of the change'], loc='upper center',
               bbox_to_anchor=(0.5 * (b0.x0 + b1.x1), b0.y0 - 0.20 * b0.height), ncol=3, fontsize=9,
               frameon=False, columnspacing=2.0, handletextpad=0.6)

    # --- C: observed adults on the withheld transects (pooled T2 + T5) vs modelled recruits, ONE panel
    # (user, 2026-09-24: 'not separating T2 and T5 ... one panel is enough')
    if 'pooled:recruits' not in set(tests['set']):
        raise SystemExit('tests table has no pooled set: rerun evaluate_recruits_vs_adults.py (2026-09-24 version)')
    held = tests.loc[tests['set'] == 'heldout:recruits', 'transects']
    held_lab = str(held.iloc[0]) if len(held) else 'withheld'
    axC = fig.add_subplot(gs[2, 1:5])
    panel_c1(axC, tests, ['pooled:recruits'], {'pooled:recruits': ''}, SPECIES, ZONES, F, rng, n_boot=args.n_boot,
             pooled=True)
    axC.set_title(f'C. Observed adults on the withheld transects ({held_lab}, T5) vs modelled recruits '
                  '(nearest-individual survey, same plots)', fontsize=10.5, loc='left', fontweight='bold')
    hC = [Patch(facecolor=F.CHG_LEGEND_FACE, edgecolor=F.BAR_EDGE, lw=0.7, hatch='///'),
          Patch(facecolor=F.CHG_LEGEND_FACE, edgecolor=F.BAR_EDGE, lw=0.7), Line2D([], [], color='k', lw=2.4)]
    axC.legend(hC, ['Observed adults (withheld from the reconstruction)', 'Modelled recruits',
                    '95 % CI of the zone boundary (observed) / replicate range (model)'], loc='upper center',
               bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=9, frameon=False, columnspacing=2.0, handletextpad=0.6)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved:', args.out)
    print('Panel C, pooled withheld transects: observed adult shares vs modelled recruit shares (%), TVD [95 % upper]:')
    d = tests[tests['set'] == 'pooled:recruits']
    for sp in SPECIES:
        r = d[d.species == sp].set_index('zone').reindex(ZONES)
        print(f'  {sp:<11} n {int(r.n_obs_species.iloc[0]):>3}  obs ' + '/'.join(f'{100 * v:.0f}' for v in r.obs_share)
              + '  sim ' + '/'.join('nan' if np.isnan(v) else f'{100 * v:.0f}' for v in r.sim_share)
              + f'  TVD {100 * r.tvd.iloc[0]:.0f} [{100 * r.tvd_upper95.iloc[0]:.0f}]')


if __name__ == '__main__':
    main()
