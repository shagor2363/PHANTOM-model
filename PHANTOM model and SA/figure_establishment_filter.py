"""The establishment filter as a transformation: from stranded propagules to recruits.

One panel per species: the elevation distribution of stranded propagules (terminal
position, established + dead; grey) and of the recruits that survived the 90-day
establishment window (species colour), absolute numbers per 4 cm bin, mean over seeds,
inside the survey transect polygons (the observation footprint of every other figure;
--whole-domain switches to the full box). The share of stranded propagules that recruit is given under the species name.

No mortality cause is attributed: the figure shows what the filter does along the
gradient, not which coded process does it.

Usage: python figure_establishment_filter.py [--summaries GLOB ...] [--out PATH]
                                             [--stratum all|internal|boundary] [--sampled]
"""
import argparse
import glob
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import cascade_metrics as cm  # noqa: E402
import tidal_dispersal_model as tdm  # noqa: E402

SPECIES = tdm.SPECIES   # 2026-09-20: 4 species (Sonneratia); summaries written before lack its keys -> empty
SPECIES_COLORS = {'Avicennia': '#E07A5F', 'Rhizophora': '#3D8B7A', 'Ceriops': '#7B2D8E', 'Sonneratia': '#3B6CC0'}
DEFAULT_GLOB = os.path.join(SCRIPT_DIR, 'results', 'knockout_cascade', 'cascade_D1T1S1E1_seed*_summary.json')
Z_MID = 0.5 * (cm.Z_EDGES[:-1] + cm.Z_EDGES[1:])
REBIN = 2          # model: 2 cm bins -> 4 cm for display


def load(paths):
    out = []
    for p in paths:
        with open(p) as f:
            out.append(json.load(f))
    return out


def _rebin(h, k=REBIN):
    h = np.asarray(h, dtype=float)
    n = len(h) // k * k
    return h[:n].reshape(-1, k).sum(axis=1)


def mean_hist(summaries, stratum, stage, sp, sampled=False):
    """Mean over seeds of the elevation histogram of a stage; sampled=True restricts the
    cohort to the survey transect polygons, otherwise the whole domain."""
    empty = np.zeros(len(cm.Z_EDGES) - 1)     # species absent from a summary (2026-09-20: 4 species) -> no propagules
    if sampled:
        return np.mean([_rebin(s['strata'][stratum]['hist_sampled'][stage].get(sp, empty)) for s in summaries], axis=0)
    return np.mean([_rebin(s['strata'][stratum]['stages'][stage]['hist'].get(sp, empty)) for s in summaries], axis=0)


SON_SHAPE_DIR = os.path.join(SCRIPT_DIR, 'results', 'knockout_son_shape')


def cache_hist(cache_paths, stratum, sp):
    """(grounded, established) elevation histograms of one species, mean over the given run caches, on the
    summaries' 2 cm bins re-binned like mean_hist. 2026-09-24: used for Sonneratia from the supply-amplified
    shape runs (results/knockout_son_shape), whose baseline cohort (~85 propagules per run) is too small to
    show a distribution; the amplification changes the sample size only (propagules do not interact)."""
    import pandas as pd
    hg, he = [], []
    for q in cache_paths:
        df = pd.read_csv(q, usecols=['species', 'source', 'outcome', 'terminal_z'])
        df = df[(df.species == sp) & df.outcome.isin(['established', 'dead'])].dropna(subset=['terminal_z'])
        if stratum != 'all':
            df = df[df.source == stratum]
        hg.append(np.histogram(df.terminal_z, bins=cm.Z_EDGES)[0])
        he.append(np.histogram(df.loc[df.outcome == 'established', 'terminal_z'], bins=cm.Z_EDGES)[0])
    return _rebin(np.mean(hg, axis=0)), _rebin(np.mean(he, axis=0))


def band_table(cache_paths):
    """Per-species boundary-influx and sand-band statistics from the run caches, mean over runs:
    rows = (boundary influx per spring day, boundary share of the stranded cohort, share of boundary
    propagules ending on the sand band, mean terminal elevation of boundary propagules, propagules
    dying on the band boundary / internal). Returns (row labels, {species: [cell strings]})."""
    import pandas as pd
    import tidal_dispersal_model as tdm
    acc = {sp: {'bshare': [], 'band_b': [], 'z_b': [], 'die_b': [], 'die_i': []} for sp in SPECIES}
    for p in cache_paths:
        df = pd.read_csv(p, usecols=['species', 'source', 'outcome', 'cause', 'terminal_z'])
        df = df[df.outcome.isin(['established', 'dead'])]
        for sp in SPECIES:
            d = df[df.species == sp]
            b = d[d.source == 'boundary']
            i = d[d.source == 'internal']
            acc[sp]['bshare'].append(len(b) / max(len(d), 1))
            acc[sp]['band_b'].append((b.cause == 'substrate').mean() if len(b) else np.nan)
            acc[sp]['z_b'].append(b.terminal_z.mean() if len(b) else np.nan)
            acc[sp]['die_b'].append((b.cause == 'substrate').sum())
            acc[sp]['die_i'].append((i.cause == 'substrate').sum())
    labels = ["share of the species' stranded cohort that came from the boundary",
              'share of its boundary propagules ending on the sand band']
    cells = {}
    for sp in SPECIES:
        a = acc[sp]
        band_b = np.asarray(a['band_b'], dtype=float)
        cells[sp] = [f'{100 * np.mean(a["bshare"]):.0f} %',
                     f'{100 * np.mean(band_b[np.isfinite(band_b)]):.0f} %' if np.isfinite(band_b).any() else 'n/a']
    return labels, cells


def plot(summaries, out, stratum='all', sampled=True, table=None, son_shape=None):
    x = _rebin(Z_MID) / REBIN
    width = (cm.Z_EDGES[1] - cm.Z_EDGES[0]) * REBIN
    n_sp = len(SPECIES)                       # 2026-09-20: 4 species (Sonneratia)
    # 2026-09-24 (user): 2 x 2 layout, panel letters, larger headings; the strip made everything too small
    ncol = 2
    nrow = int(np.ceil(n_sp / ncol))
    fig_w = 11.0
    if table is None:
        fig, axes = plt.subplots(nrow, ncol, figsize=(fig_w, 4.6 * nrow), sharex=True, sharey=True)
        axes = list(np.ravel(axes))
    else:
        fig = plt.figure(figsize=(fig_w, 4.6 * nrow + 1.9))
        gs = fig.add_gridspec(nrow + 1, ncol, height_ratios=[4.6] * nrow + [1.5], hspace=0.55, wspace=0.16)
        axes = [fig.add_subplot(gs[r, c]) for r in range(nrow) for c in range(ncol)][:n_sp]
        ax_t = fig.add_subplot(gs[nrow, :])
        ax_t.set_axis_off()
    for c, sp in enumerate(SPECIES):
        ax = axes[c]
        if son_shape and sp == 'Sonneratia' and not sampled:
            h_gr, h_es = cache_hist(son_shape, stratum, sp)
        else:
            h_gr = mean_hist(summaries, stratum, 'grounded', sp, sampled)
            h_es = mean_hist(summaries, stratum, 'established', sp, sampled)
        # 2026-09-24 (user: the Sonneratia counts, ~6 per bin against thousands, 'reviewers will flag it'):
        # every panel shows the share of the species' stranded cohort per 4 cm bin, so the panels are
        # comparable whatever the species' supply; the recruit share in the title is unchanged
        pct = 100 * h_es.sum() / max(h_gr.sum(), 1)
        tot = max(h_gr.sum(), 1)
        n_runs = len(son_shape) if (son_shape and sp == 'Sonneratia' and not sampled) else len(summaries)
        # counts are printed to the console for the caption (the panel annotation was removed at the user's request, 2026-09-24)
        print(f'  {sp:<11} stranded {h_gr.sum():9.1f}  recruits {h_es.sum():9.1f}  ({pct:.0f} %), mean per run over {n_runs} runs')
        h_gr_raw, h_es_raw = h_gr.copy(), h_es.copy()
        h_gr, h_es = 100 * h_gr / tot, 100 * h_es / tot
        ax.bar(x, h_gr, width=width, color='#d9d9d9', edgecolor='#999999', lw=0.4,
               label='All stranded propagules')
        ax.bar(x, h_es, width=width, color=SPECIES_COLORS[sp], edgecolor='none', alpha=0.85,
               label='Established recruits (subset of the stranded propagules)')
        for zb in (1.4, 1.8):
            ax.axvline(zb, color='#888888', ls='--', lw=0.8)
        ax.set_xlim(1.0, 2.2)
        ax.set_ylim(bottom=0)                 # 2026-09-20: a species without propagules keeps a 0-based axis
        ax.set_xticks([1.0, 1.4, 1.8, 2.2])
        ax.set_xlabel('Elevation (m)' if c >= n_sp - ncol else '', fontsize=10)
        ax.set_ylabel('Share of stranded cohort per bin (%)' if c % ncol == 0 else '', fontsize=10)
        ax.set_title(f'$\it{{{sp}}}$' + chr(10) + f'Overall establishment: {pct:.0f}%', fontsize=12,
                     color=SPECIES_COLORS[sp], pad=34)
        ax.text(-0.06, 1.24, 'ABCDEFGH'[c], transform=ax.transAxes, fontsize=13, fontweight='bold', ha='left', va='bottom')
        ax.tick_params(labelsize=9.5)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        for z, (lo, hi) in {'C': (1.0, 1.4), 'B': (1.4, 1.8), 'A': (1.8, 2.2)}.items():
            # 2026-09-25 (user): survival of the stranded cohort within the zone (recruits / stranded in the zone's bins)
            inz = (x >= lo) & (x < hi) if z != 'A' else (x >= lo)
            s_z, r_z = h_gr_raw[inz].sum(), h_es_raw[inz].sum()
            surv = f'{100 * r_z / s_z:.0f}% established' if s_z > 0 else 'no strandings'
            # annotation strip ABOVE the plotting area (user, 2026-09-25): zone names grey, percentages in the species colour
            ax.text(0.5 * (lo + hi), 1.085, f'Zone {z}', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=9, color='#666666')
            ax.text(0.5 * (lo + hi), 1.02, surv, transform=ax.get_xaxis_transform(), ha='center',
                    va='bottom', fontsize=9, color=SPECIES_COLORS[sp], fontweight='bold')
            print(f'    {sp:<11} zone {z}: stranded {s_z:9.0f}  recruits {r_z:9.0f}  survival {100 * r_z / max(s_z, 1):5.1f} %')
    ymax = max(a.get_ylim()[1] for a in axes)
    for a in axes:
        a.set_ylim(0, ymax)                   # one shared scale: shares are comparable across species
    handles, labels = axes[0].get_legend_handles_labels()
    if table is None:
        fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.01))
        fig.tight_layout(rect=(0, 0.05, 1, 1))
    else:
        ax_t.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.30), ncol=2, frameon=False,
                    fontsize=10)
        row_labels, cells = table
        tab = ax_t.table(cellText=[[cells[sp][r] for sp in SPECIES] for r in range(len(row_labels))],
                         rowLabels=row_labels, colLabels=[f'$\it{{{sp}}}$' for sp in SPECIES],
                         cellLoc='center', rowLoc='left', loc='center', colWidths=[0.14] * n_sp, bbox=[0.40, 0.0, 0.145 * n_sp, 0.62])
        tab.auto_set_font_size(False)
        tab.set_fontsize(9)
        for (r, c), cell in tab.get_celld().items():
            cell.set_edgecolor('#bbbbbb')
            cell.set_linewidth(0.6)
            if r == 0:
                cell.set_text_props(color=SPECIES_COLORS[SPECIES[c]] if c >= 0 else 'k', fontweight='bold')
                cell.set_facecolor('#f4f4f4')
        fig.subplots_adjust(top=0.93, bottom=0.03, left=0.08, right=0.98)
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('Saved:', out)


def print_table(summaries, stratum='all', year='1999'):
    obs = cm.load_observed(year)
    tests_gr = cm.chi2_tests({k: float(np.mean([s['strata'][stratum]['sampled_densities_grounded'][k] for s in summaries]))
                              for k in summaries[0]['strata'][stratum]['sampled_densities_grounded']},
                             summaries[0]['n_quadrats'], obs['quadrats'])
    tests_es = cm.chi2_tests({k: float(np.mean([s['strata'][stratum]['sampled_densities'][k] for s in summaries]))
                              for k in summaries[0]['strata'][stratum]['sampled_densities']},
                             summaries[0]['n_quadrats'], obs['quadrats'])
    def _mean(vals):
        vals = np.asarray(vals, dtype=float)
        return float(np.mean(vals[np.isfinite(vals)])) if np.isfinite(vals).any() else np.nan

    def _p(tests, sp):
        t = tests.get(sp)
        return t['p'] if t and t.get('p') is not None else np.nan

    for sp in SPECIES:      # 2026-09-20: 4 species; a species missing from the summaries or the survey prints nan
        sh_gr = {z: _mean([s['strata'][stratum]['sampled_shares_grounded'].get(sp, {}).get(z, np.nan) for s in summaries])
                 for z in 'CBA'}
        sh_es = {z: _mean([s['strata'][stratum]['sampled_shares'].get(sp, {}).get(z, np.nan) for s in summaries])
                 for z in 'CBA'}
        surv = {z: _mean([(s['strata'][stratum]['survival_zone'].get(sp) or {}).get(z) or np.nan for s in summaries])
                for z in 'CBA'}
        obs_sh = obs['shares'].get(sp, {})
        print(f"{sp:<11} stranded C/B/A {' / '.join(f'{100 * sh_gr[z]:4.0f}' for z in 'CBA')}  "
              f"recruited {' / '.join(f'{100 * sh_es[z]:4.0f}' for z in 'CBA')}  "
              f"observed {' / '.join(f'{100 * obs_sh.get(z, np.nan):4.0f}' for z in 'CBA')}  "
              f"survival C/B/A {' / '.join(f'{100 * surv[z]:3.0f}%' for z in 'CBA')}  "
              f"p stranded {_p(tests_gr, sp):.3g}  p recruited {_p(tests_es, sp):.3g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--summaries', nargs='*', default=None)
    ap.add_argument('--out', default=os.path.join(SCRIPT_DIR, 'results', 'figures', 'fig_establishment_filter.png'))
    ap.add_argument('--stratum', default='all', choices=['all', 'internal', 'boundary'])
    ap.add_argument('--transects', action='store_true',
                    help='restrict to the survey transect polygons (default since 2026-09-06: whole domain, so that '
                         'the zone C strandings that cannot root are visible; the transects only reach 1.27 m elevation)')
    ap.add_argument('--whole-domain', action='store_true', help='(default now; kept for compatibility)')
    ap.add_argument('--year', default='1999', help='survey year for the printed table')
    ap.add_argument('--son-shape-dir', default=SON_SHAPE_DIR,
                    help="Sonneratia's histograms from the supply-amplified shape runs in this directory "
                         "(default results/knockout_son_shape; '' = use the baseline summaries)")
    ap.add_argument('--no-table', action='store_true',
                    help='omit the boundary-influx / sand-band table under the panels (computed from the run caches)')
    args = ap.parse_args()
    paths = []
    for pat in (args.summaries or [DEFAULT_GLOB]):
        paths.extend(sorted(glob.glob(pat)) or ([pat] if os.path.exists(pat) else []))
    if not paths:
        sys.exit('no summaries found')
    summaries = load(paths)
    if 'sampled_densities_grounded' in summaries[0]['strata']['all']:
        print_table(summaries, args.stratum, args.year)
    table = None
    if not args.no_table:
        caches = [p.replace('_summary.json', '.csv') for p in paths]
        caches = [c for c in caches if os.path.exists(c)]
        if caches:
            table = band_table(caches)
            for r, lab in enumerate(table[0]):
                print(f'{lab:<68} ' + '  '.join(f'{table[1][sp][r]:>12}' for sp in SPECIES))
    son_shape = sorted(glob.glob(os.path.join(args.son_shape_dir, 'knockout_baseline_seed*.csv'))) if args.son_shape_dir else []
    if son_shape:
        print(f'  Sonneratia panel from {len(son_shape)} shape-run caches in {args.son_shape_dir} (distribution only)')
    plot(summaries, args.out, args.stratum, sampled=args.transects, table=table, son_shape=son_shape)


if __name__ == '__main__':
    main()
