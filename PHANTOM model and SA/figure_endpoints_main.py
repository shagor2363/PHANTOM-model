"""Main-text endpoints figure (user's specification, 2026-09-24): one question, 'how does equalising the transport
traits change stranding under each source configuration?', answered by comparing the two panels of each row.

3 x 2 grid. Rows: Boundary only / Parent trees only / Both sources (short horizontal row labels). Columns:
Species-specific transport traits / Equalised transport traits (headings written once). Each panel: four solid
species curves (kernel density of the stranding elevation of all grounded propagules, established + dead), no
fills, no panel titles; zone boundaries dashed with 'C: seaward / B: middle / A: landward' across the top; panel
letters A-F; one legend with full species names; one shared y label 'Probability density (m^-1)' and one shared
x label 'Stranding elevation (m)'. Data: results/knockout_final baseline and
all_transport_eq, seeds 42/49/56 pooled; Sonneratia from the shape runs, caption must say so).

    python figure_endpoints_main.py [--recruits-only] [--ymax 7] [--out ...]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
STAND = os.path.join(SCRIPT_DIR, 'results', 'reconstruction', 'stand3T_T4T3T1T6_single', 'demo_reconstructed_stand3T.csv')
os.environ.setdefault('KNOCKOUT_TREES', STAND)
import transport_treatments as tt  # noqa: E402

KNOCKOUT_DIR = os.path.join(SCRIPT_DIR, 'results', 'knockout_final')
SON_SHAPE_DIR = os.path.join(SCRIPT_DIR, 'results', 'knockout_son_shape')
SEEDS = (42, 49, 56)
TREATMENTS = [('baseline', ''), ('all_transport_eq', '')]
SPECIES = list(tt.SPECIES_ORDER)
Z = tt.Z_GRID
YMAX = 7.0   # shared density ceiling of every panel
ROWS = [('boundary', 'Boundary only'), ('internal', 'Parent trees only'), ('all', 'Both sources')]
COLS = [('baseline', 'Species-specific transport traits'), ('all_transport_eq', 'Equalised transport traits')]
ZONE_TOP = ((1.2, 'C: seaward'), (1.6, 'B: middle'), (2.0, 'A: landward'))


def load_terminal(outcomes):
    """{(treatment, species, stratum): terminal_z array} pooled over seeds; Sonneratia from the shape runs."""
    out = {}
    for tr, _ in TREATMENTS:
        frames = []
        for seed in SEEDS:
            d = pd.read_csv(os.path.join(KNOCKOUT_DIR, f'knockout_{tr}_seed{seed}.csv'),
                            usecols=['species', 'source', 'outcome', 'terminal_z'])
            d = d[d.species != 'Sonneratia']
            s = pd.read_csv(os.path.join(SON_SHAPE_DIR, f'knockout_{tr}_seed{seed}.csv'),
                            usecols=['species', 'source', 'outcome', 'terminal_z'])
            frames += [d, s[s.species == 'Sonneratia']]
        d = pd.concat(frames)
        d = d[d.outcome.isin(outcomes)].dropna(subset=['terminal_z'])
        for sp in SPECIES:
            for stratum in ('boundary', 'internal', 'all'):
                sel = d[(d.species == sp) & ((d.source == stratum) if stratum != 'all' else True)]
                out[(tr, sp, stratum)] = sel.terminal_z.values
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--recruits-only', action='store_true')
    ap.add_argument('--ymax', type=float, default=YMAX, help='shared density ceiling (0 = data-driven)')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scipy.stats import gaussian_kde

    outcomes = ('established',) if args.recruits_only else ('established', 'dead')
    out = args.out or os.path.join(SCRIPT_DIR, 'results', 'figures',
                                   'fig_endpoints_main_recruits.png' if args.recruits_only else 'fig_endpoints_main.png')
    term = load_terminal(outcomes)

    fig, axes = plt.subplots(3, 2, figsize=(10.5, 9.6), sharex=True, sharey=True,
                             gridspec_kw=dict(left=0.11, right=0.98, top=0.90, bottom=0.11, hspace=0.18, wspace=0.10))
    letters = iter('ABCDEF')
    for i, (stratum, rlab) in enumerate(ROWS):
        for j, (tr, _) in enumerate(COLS):
            ax = axes[i, j]
            # light species fill (10 %) under every curve, opaque outlines drawn above all fills (user, 2026-09-24)
            dens = {}
            for sp in SPECIES:
                data = term[(tr, sp, stratum)]
                if len(data) < 2:
                    continue
                dens[sp] = gaussian_kde(data)(Z)
                ax.fill_between(Z, dens[sp], color=tt.SPECIES_COLORS[sp], alpha=0.10, lw=0, zorder=1)
            for sp, d in dens.items():
                ax.plot(Z, d, color=tt.SPECIES_COLORS[sp], lw=2.0, zorder=3)
            for zb in tt.ZONE_BOUNDARIES:
                ax.axvline(zb, color='#888888', ls='--', lw=0.8, alpha=0.6)
            ax.text(0.015, 0.97, next(letters), transform=ax.transAxes, ha='left', va='top',
                    fontsize=12, fontweight='bold')
            if i == 0:
                for xz, lab in ZONE_TOP:
                    ax.text(xz, 1.02, lab, transform=ax.get_xaxis_transform(), ha='center', va='bottom',
                            fontsize=8.5, color='0.35')
            ax.set_xlim(1.0, 2.2)
            ax.set_xticks([1.0, 1.4, 1.8, 2.2])
            ax.grid(True, axis='y', lw=0.4, alpha=0.4)
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
        # short horizontal row label, left of the row
        # rotated row label just left of the tick labels (user, 2026-09-24: no white space column)
        axes[i, 0].annotate(rlab, xy=(0, 0.5), xycoords='axes fraction', xytext=(-34, 0), textcoords='offset points',
                            rotation=90, ha='center', va='center', fontsize=10.5, fontweight='bold', color='0.2')
    for j, (_, clab) in enumerate(COLS):
        axes[0, j].annotate(clab, xy=(0.5, 1.0), xycoords='axes fraction', xytext=(0, 22), textcoords='offset points',
                            ha='center', va='bottom', fontsize=11.5, fontweight='bold', color='0.2')
    if args.ymax and args.ymax > 0:
        axes[0, 0].set_ylim(0, args.ymax)
    handles = [Line2D([], [], color=tt.SPECIES_COLORS[sp], lw=2.2, label=sp) for sp in SPECIES]
    leg = fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=10.5,
                     bbox_to_anchor=(0.5, 0.0))
    for t in leg.get_texts():
        t.set_fontstyle('italic')
    fig.supxlabel('Stranding elevation (m)', fontsize=11, y=0.055)
    fig.supylabel('Probability density (m$^{-1}$)', fontsize=11, x=0.015)
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved:', out)


if __name__ == '__main__':
    main()
