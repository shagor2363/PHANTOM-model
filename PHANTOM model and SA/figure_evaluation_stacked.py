"""Calibration and evaluation figure with paired 100 % stacked composition bars.

Two figures.

Calibration figure (fig_evaluation_stacked.png): 1999 calibration target, per species a pair of stacked
   bars (observed | simulated) showing the relative density across zones C / B / A (share of the species' mean quadrat density). A stacked
   composition is fixed by its two internal boundaries (share of C, share of C + B), so the uncertainty is
   drawn there: a black bracket spanning the 95 % bootstrap interval of each cumulative share (quadrats
   resampled within zones; for the model also one random replicate run per draw, virtual_survey.py). A
   boundary with no bracket is fixed at 0. The parity form of the same comparison (the goodness-of-fit
   figure of the paper) is figure_goodness_of_fit_1999.py.

Evaluation figure (fig_zonation_change_1997.png):
A. 1997 withheld juvenile evaluation: stacked composition bars as above for the 1997 survey.
B. Change from the 1997 baseline, per species and zone: the observed change (1999 survey - 1997 survey)
   next to the model's deviation from 1997 (model - 1997), diverging bars about zero with 95 % intervals
   (bootstrap of the difference). Bars on the same side of zero with similar length = the model departs
   from 1997 in the same direction and by the same amount as the 1999 survey itself did, i.e. its
   distance from 1997 is the between-survey variability.

Usage: python figure_evaluation_stacked.py [--n-boot 10000]
Output: results/figures/fig_evaluation_stacked.png and results/figures/fig_zonation_change_1997.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from virtual_survey import virtual_survey, SPECIES, ZONES  # noqa: E402


def tvd(a, b):
    """Total variation distance between two share vectors."""
    return 0.5 * np.abs(np.asarray(a) - np.asarray(b)).sum(axis=-1)


# seaward -> landward gradient: pale blue, soft khaki-green, dark olive (user's palette, 2026-09-17;
# adjacent-pair separation checked, normal vision dE 18.7, deutan 18.0); whiskers stay black
ZONE_FACE = {'C': '#d0e3f3', 'B': '#adb27c', 'A': '#4b4e36'}
ZONE_LABEL = {'C': 'Zone C (seaward)', 'B': 'Zone B', 'A': 'Zone A (landward)'}
BAR_EDGE = '#333333'
ZONE_TICK = {'C': 'C\n(seaward)', 'B': 'B', 'A': 'A\n(landward)'}
# change panel: bars carry the zone colour; the observed change is hatched, the model solid
CHG_STYLE = {'obs': dict(hatch='///'), 'model': dict()}
CHG_LABEL = {'obs': 'Observed change (1999 survey \u2212 1997 survey)', 'model': 'Model \u2212 1997 survey'}
CHG_LEGEND_FACE = '#c9c9c9'   # neutral fill for the legend patches (colour on the bars = zone)


def stacked_panel(ax, vs, title, show_ylabel=True):
    """Paired stacked bars per species with boundary brackets: the observed and simulated bars of a species
    stand close together, species groups are separated by a wider gap."""
    w = 0.40          # bar width
    gap = 0.03        # gap inside the observed | simulated pair
    x = 1.3 * np.arange(len(SPECIES))   # species centres (gap between species groups = 1.3 - 2w - gap)
    for i, sp in enumerate(SPECIES):
        for k, (key, boot, name) in enumerate((('obs', vs['obs_boot'][sp], 'Observed'),
                                               ('sim', vs['sim_boot'][sp], 'Simulated'))):
            est = vs[key][sp][0]
            xc = x[i] + (k - 0.5) * (w + gap)
            bottom = 0.0
            for j, z in enumerate(ZONES):
                ax.bar(xc, 100 * est[j], w, bottom=100 * bottom, facecolor=ZONE_FACE[z], edgecolor=BAR_EDGE, lw=0.7,
                       zorder=2)
                bottom += est[j]
            # brackets at the two internal boundaries: cumulative shares C and C+B
            cum = np.cumsum(boot, axis=1)[:, :2]
            # the C/B boundary bracket on the left third of the bar, the B/A bracket on the right third,
            # so the two stay apart when the boundaries are close
            for b, xb in enumerate((xc - w / 4, xc + w / 4)):
                lo, hi = np.percentile(cum[:, b], [2.5, 97.5])
                if hi - lo < 1e-9:
                    continue
                ax.plot([xb, xb], [100 * lo, 100 * hi], color='k', lw=2.4, zorder=6, solid_capstyle='butt', clip_on=False)
                ax.plot([xb - 0.07, xb + 0.07], [100 * lo] * 2, color='k', lw=2.4, zorder=6, clip_on=False)
                ax.plot([xb - 0.07, xb + 0.07], [100 * hi] * 2, color='k', lw=2.4, zorder=6, clip_on=False)
            ax.text(xc, -4.5, name, ha='center', va='top', fontsize=8.5, color='#333333')
    ax.set_xticks(x)
    ax.set_xticklabels([f'$\\it{{{sp}}}$' for sp in SPECIES], fontsize=10.5)
    ax.tick_params(axis='x', pad=22, length=0)
    ax.set_xlim(x[0] - 0.65, x[-1] + 0.65)
    ax.set_ylim(-2.5, 101)
    ax.spines['bottom'].set_position(('data', -2.5))
    if show_ylabel:
        ax.set_ylabel("Relative density across zones (%)\n(share of the species' mean quadrat density)", fontsize=9.5)
    ax.set_title(title, fontsize=10.5, loc='left', fontweight='bold')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def _whisker(ax, x, lo, hi, cap=0.07):
    """95 % interval as a black line with a white halo (stays visible on the dark model bar)."""
    for col, lw, z in (('white', 3.2, 5), ('k', 1.3, 6)):
        ax.plot([x, x], [lo, hi], color=col, lw=lw, zorder=z, solid_capstyle='butt', clip_on=False)
        for yv in (lo, hi):
            ax.plot([x - cap, x + cap], [yv, yv], color=col, lw=lw, zorder=z, clip_on=False)


def change_panel(ax, sp, vs97, vs99, show_ylabel=True):
    """One species: per zone the observed 1999 - 1997 change next to the model - 1997 deviation
    (percentage points), diverging bars about zero with 95 % bootstrap intervals of the difference."""
    x = 1.15 * np.arange(len(ZONES))   # zone centres; the observed | model pair sits close together
    w = 0.36
    gap = 0.03
    e97 = np.asarray(vs97['obs'][sp][0])
    b97 = vs97['obs_boot'][sp]
    series = {'obs': (np.asarray(vs99['obs'][sp][0]) - e97, vs99['obs_boot'][sp] - b97),
              'model': (np.asarray(vs99['sim'][sp][0]) - e97, vs99['sim_boot'][sp] - b97)}
    reach = 55.0   # the y range of the calibrated species (unchanged); widened only where a whisker needs it
    for k, key in enumerate(('obs', 'model')):
        est, boot = series[key]
        xs = x + (k - 0.5) * (w + gap)
        for xi, e, z in zip(xs, est, ZONES):
            # the hatch takes the edge colour: white on the dark zone A fill, dark grey elsewhere
            edge = 'white' if (key == 'obs' and z == 'A') else BAR_EDGE
            ax.bar(xi, 100 * e, w, facecolor=ZONE_FACE[z], edgecolor=edge, lw=0.7, zorder=2, **CHG_STYLE[key])
        lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
        for xi, l, h in zip(xs, lo, hi):
            if h - l > 1e-9:
                _whisker(ax, xi, 100 * l, 100 * h)
                ext = 100 * max(abs(l), abs(h))
                if ext > 55:
                    reach = max(reach, 10 * np.ceil((ext + 3) / 10))
    ax.axhline(0, color='k', lw=0.9, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([ZONE_TICK[z] for z in ZONES], fontsize=9.5)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(x[0] - 0.7, x[-1] + 0.7)
    ax.set_ylim(-reach, reach)   # 2026-09-21: Sonneratia's intervals (11 / 3 occupied quadrats) exceed +-55
    ax.grid(True, axis='y', lw=0.4, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    if show_ylabel:
        ax.set_ylabel("Change in the relative density across zones since 1997\n(percentage points of the species' mean quadrat density)", fontsize=9.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    return series


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n-boot', type=int, default=10000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    if args.out is None:
        args.out = os.path.join(SCRIPT_DIR, 'results', 'figures', 'fig_evaluation_stacked.png')
    vs99 = virtual_survey('1999', args.n_boot, args.seed)
    vs97 = virtual_survey('1997', args.n_boot, args.seed)

    rows = []
    for sp in SPECIES:
        d = tvd(vs97['obs_boot'][sp], vs97['sim_boot'][sp])
        d_obs = tvd(vs99['obs_boot'][sp], vs97['obs_boot'][sp])
        rows.append({'species': sp, 'model': float(tvd(vs97['obs'][sp][0], vs97['sim'][sp][0])),
                     'model_ub95': float(np.percentile(d, 95)),
                     'change': float(tvd(vs99['obs'][sp][0], vs97['obs'][sp][0])),
                     'change_lo': float(np.percentile(d_obs, 2.5)), 'change_hi': float(np.percentile(d_obs, 97.5))})

    for r in rows:   # discrepancy numbers for the caption (no longer drawn)
        print(f"{r['species']:<11} D model vs 1997 {100 * r['model']:.1f} (95 % ub {100 * r['model_ub95']:.1f})   "
              f"observed change 1999 vs 1997 {100 * r['change']:.1f} [{100 * r['change_lo']:.1f}, {100 * r['change_hi']:.1f}]")

    def zone_legend(ax):
        """Legend of a stacked panel: the three zones, then the bracket as it appears on the bars."""
        h = [Patch(facecolor=ZONE_FACE[z], edgecolor=BAR_EDGE, lw=0.7) for z in ZONES]
        h.append(plt.Line2D([], [], color='k', lw=2.4))
        lab = [ZONE_LABEL[z] for z in ZONES] + ['95 % CI of the zone boundary']
        ax.legend(h, lab, loc='upper center', bbox_to_anchor=(0.5, -0.20), ncol=4, fontsize=9, frameon=False,
                  columnspacing=2.0)

    def species_row(fig, gs_row, draw, title, offset):
        """One sub-panel per species in one gridspec row; draw(ax, i, sp) fills each; panel title above the row."""
        sub = gs_row.subgridspec(1, len(SPECIES), wspace=0.12)      # 2026-09-20: 4 species (Sonneratia)
        axes = [fig.add_subplot(sub[0, i]) for i in range(len(SPECIES))]
        for i, sp in enumerate(SPECIES):
            draw(axes[i], i, sp)
            # 2026-09-21: a panel whose y range differs from the first (Sonneratia, wide intervals) keeps its labels
            if i > 0 and np.allclose(axes[i].get_ylim(), axes[0].get_ylim()):
                axes[i].tick_params(axis='y', labelleft=False)
        bb = axes[0].get_position()
        fig.text(bb.x0, bb.y1 + offset, title, fontsize=10.5, fontweight='bold', ha='left', va='bottom')
        return axes

    # ---- calibration figure: the 1999 stacked composition alone
    fig, axA = plt.subplots(figsize=(8.0 + 2.0 * (len(SPECIES) - 3), 5.4))   # 2026-09-20: 4 species (Sonneratia)
    stacked_panel(axA, vs99, '1999 calibration target: observed vs simulated zonation')
    zone_legend(axA)
    fig.savefig(args.out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved:', args.out)

    # ---- evaluation figure: A. 1997 stacked composition, B. change from 1997 (observed 1999 vs model)
    fig2 = plt.figure(figsize=(11.5 + 2.5 * (len(SPECIES) - 3), 8.6))
    gs2 = fig2.add_gridspec(2, 6, height_ratios=[1.35, 0.95], hspace=0.62, wspace=0.10, top=0.95, bottom=0.06)
    axB = fig2.add_subplot(gs2[0, 1:5])
    stacked_panel(axB, vs97, 'A. 1997 withheld juvenile evaluation')
    zone_legend(axB)
    series_all = {}

    def draw_change(ax, i, sp):
        series_all[sp] = change_panel(ax, sp, vs97, vs99, show_ylabel=(i == 0))
        ax.set_title(f'$\\it{{{sp}}}$', fontsize=10.5, pad=4)
    axD = species_row(fig2, gs2[1, :], draw_change, 'B. Change since the 1997 survey: observed (1999) vs model', 0.034)
    hD = [Patch(facecolor=CHG_LEGEND_FACE, edgecolor=BAR_EDGE, lw=0.7, **CHG_STYLE[k]) for k in ('obs', 'model')]
    hD.append(plt.Line2D([], [], color='k', lw=1.3))
    # legend centred under the whole species row (2026-09-20: the row has one panel per species)
    b0, b1 = axD[0].get_position(), axD[-1].get_position()
    fig2.legend(hD, [CHG_LABEL['obs'], CHG_LABEL['model'], '95 % CI of the change'],
                loc='upper center', bbox_to_anchor=(0.5 * (b0.x0 + b1.x1), b0.y0 - 0.22 * b0.height), ncol=3,
                fontsize=9, frameon=False, columnspacing=2.0, handletextpad=0.6)
    out_d = os.path.join(os.path.dirname(args.out), 'fig_zonation_change_1997.png')
    fig2.savefig(out_d, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig2)
    print('Saved:', out_d)

    # numbers for the caption of the evaluation figure
    obs_chg = np.concatenate([100 * series_all[sp]['obs'][0] for sp in SPECIES])
    mod_chg = np.concatenate([100 * series_all[sp]['model'][0] for sp in SPECIES])
    same = np.sum((np.sign(obs_chg) == np.sign(mod_chg)) | (np.abs(obs_chg) < 2) & (np.abs(mod_chg) < 2))
    print('Change from 1997 (pp), observed vs model per species x zone (C, B, A):')
    for i, sp in enumerate(SPECIES):
        print(f"  {sp:<11} obs {np.round(obs_chg[3*i:3*i+3], 1)}   model {np.round(mod_chg[3*i:3*i+3], 1)}")
    print(f'  same direction (or both within 2 pp of zero): {same}/{len(obs_chg)}; Pearson r = {np.corrcoef(obs_chg, mod_chg)[0, 1]:.2f}; '
          f'mean |obs - model| = {np.mean(np.abs(obs_chg - mod_chg)):.1f} pp')


if __name__ == '__main__':
    main()
