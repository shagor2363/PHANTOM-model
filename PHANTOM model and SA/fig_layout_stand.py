"""The paper's transect map next to the reconstructed (fringe-clipped) stand in the model frame.

Left : the map panel of Dahdouh-Guebas et al. (2004, Fig. 6a), cropped to its frame, with the
       paper's legend below it.
Right: the reconstructed adult stand drawn as in the reconstruction validation map
       (species colours of validate.py, measured trees with a black edge, transect belts,
       surveyed envelope dashed / dotted, zone boundaries), with the area seaward of the forest
       edge (line through the transect seaward ends, T4 at the start of its surveyed part; no
       synthetic tree there) shaded in the paper's sand colour, and its own species legend. Drawn with the sea at the
       bottom (model orientation, user's choice); note that this view is the printed map's mirror
       image (a proper rotation with the sea at the bottom would put T6 on the left).

2026-09-21 (user: 'rotate the dahdouh guebas panel so both panels are in the same direction, not
flipped'): --paper-orient turns panel A so that it, too, has the sea at the bottom. DEFAULT
'flipud' = mirror top-bottom, which is the transform the model frame itself uses
(transect_geometry._px_to_model: scale, mirror, rotate), so both panels show the land-sea gradient
in the same direction AND the transects in the same left-to-right order (T4 left ... T6 right).
'rot180' (a proper rotation, tried first at the user's wording) also puts the sea at the bottom but
reverses the alongshore order with respect to the model panel, which reads as mirrored - rejected
by the user. 'none' leaves the printed orientation (sea at the top).

2026-09-20 (user): the model runs on the stand reconstructed from THREE parallel transects (T4, T2,
T6, T3); T1 is withheld from the reconstruction (blind tests of recruits and adults) and is
drawn as red dashed belts without measured trees. The five-transect stand is kept as
fig_layout_stand_5transects.png (--holdout with no ids reproduces it).

Usage: python fig_layout_stand.py [--page PAGE] [--stand CSV] [--out DIR] [--holdout 3 1]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reconstruction'))
import recon_config as C  # noqa: E402
import recon_geometry as G  # noqa: E402
import transect_geometry as tg  # noqa: E402
from postprocess import load_geom  # noqa: E402
from make_substrate_map import (forest_fringe_line, PAPER_COLORS, LEGEND_ENTRIES, SPECIES_KEY,  # noqa: E402
                                load_panel_a, CROP)
from tidal_dispersal_model import (X_MIN, X_MAX, Y_MIN, Y_MAX, zone_boundary_lines,  # noqa: E402
                                   SUBSTRATE_MAP_FILE, SUBSTRATE_SUITABILITY)
from matplotlib.colors import ListedColormap  # noqa: E402

DEFAULT_PAGE = os.path.join(C.FIG_DIR, 'dahdouh_guebas_2004_page09_300dpi.png')
STAND3T_DIR = os.path.join(C.OUT_ROOT, 'stand3T_T4T3T1T6_single')   # 2026-09-22: four-transect ensemble T4, T3, T1, T6; T2 withheld (the model's stand); single emergent-count realisations
# The three-transect version (T4, T2, T6) is kept at 'stand3T_T4T2T6_single': it could not reproduce the 1999 zonation under any
# transport setting because T1 and T3 carry 80 % of the measured Rhizophora seed mass of zone B (verification record, 2026-09-21).
DEFAULT_STAND = None     # resolved in main(): the representative ensemble member of STAND3T_DIR
DEFAULT_GEOM = None      # resolved in main(): that member's geom directory
DEFAULT_HOLDOUT = [2]   # 2026-09-22: T2 is the withheld transect
# 2026-09-21: the band is filled with a LIGHTENED version of the paper's sand yellow (the printed
# '#fdf109' at 35 % over white), painted SOLID and above the tile hatching, so the domain background
# stays white and the band is the only tinted area of the panel (user: 'domain background should be
# white. Now sand bank cant be distinguished'; the pale cream tried before was invisible under the
# hatch). Lighter than the printed yellow because the user found that one too close to Sonneratia
# ('make the sand band color in the panel more lighter'): dE76 41 from white, 39 from Sonneratia.
SAND = PAPER_COLORS['Sand (incl. sandy stretch)']
SAND_FILL = '#fefaa9'
SAND_ZORDER = 0.6          # above the hatched tile remainders (0.5), below the white reference belts (0.8)
# Withheld transects: deep violet instead of red (2026-09-21, user: 'make the withheld color not
# exactly red, because of daltonism'). Vienot-simulated separation from every other ink of the panel
# (black retained belts, the four species, the sand band): dE76 >= 58 normal vision, >= 30 deuteranopia,
# >= 50 protanopia; red HOLDOUT_COLOR fell to 19 (vs Ceriops ochre) and 27 (vs black). The dashed line
# style carries the distinction on its own for achromatopsia.
HOLDOUT_COLOR = '#6A3D9A'
# 2026-09-21 (user: 'use the same color for all 4 species according to dahdouh-guebas panel'): the
# trees of the model panel carry the species colours of the published legend, so one colour means
# one species in both panels and the shared legend below panel A applies to the stand as well.
# validate.SP_COLORS (the reconstruction's own palette) is unchanged and still used by every other
# reconstruction figure.
SP_COLORS = {sp: PAPER_COLORS[SPECIES_KEY[sp]] for sp in C.SPECIES_LEVELS}
# deepened from the printed swatch (#d2dc34) at the user's request (2026-09-21, same message), so the
# Sonneratia trees stand out against the sand band: dE76 62 from the band wash, 53 from Rhizophora.
SP_COLORS['Sonneratia'] = '#b6bf24'


def fringe_line():
    """Forest edge (clip line) sampled across the box; see make_substrate_map.forest_fringe_line."""
    xs = np.linspace(X_MIN, X_MAX, 400)
    return xs, forest_fringe_line(xs)


# The reconstruction design used to be drawn over the stand: the hatched EXTRAPOLATED remainder of
# every corridor x band tile with its dashed outline, the band labels ('0-60 m', '60-120 m',
# '>= 120 m') and the white SAMPLED reference belts. All of it was removed on 2026-09-21 at the
# user's request ('White background is enough', 'also remove 0-60 m and the other two', 'remove
# sampled portion and its legend'); the design is described in the caption / METHODS instead, and
# reconstruction/validate.py still draws it.


# The orthogonal transect T5 (evaluation only, never used in any reconstruction) was drawn here as a
# blue diagonal; removed 2026-09-21 at the user's request ('remove the T5 diagonal line and its
# legend'); re-added 2026-09-24 as a dotted grey line on both panels and REVERTED the same day
# ('Revert back, remove the T5 line from both panels'). The helpers t5_centre_line / draw_t5_on_map
# stay unused. It is the hold-out transect of evaluate_recruits_vs_adults.py.


# 2026-09-24 (user: 'Add T5 on both ... panels ... in a subtle way because it is withheld'): the orthogonal
# transect T5 is drawn as a thin dotted grey centre line (no belt, no colour) on both panels, labelled 'T5'
# in small grey type, with one legend entry. It was never used by any reconstruction (evaluation only).
T5_COLOR = '0.35'
T5_STYLE = dict(color=T5_COLOR, ls=(0, (1.5, 2.5)), lw=1.0, zorder=4.5)


def t5_centre_line(n=400):
    """(x, y) of T5's tape from u = 0 to its tape length in the model frame, clipped to the box."""
    us = np.linspace(0.0, tg.transect_frame(5)['tape_length'], n)
    pts = np.array([tg.transect_position(5, u) for u in us])
    inside = (pts[:, 0] >= X_MIN) & (pts[:, 0] <= X_MAX) & (pts[:, 1] >= Y_MIN) & (pts[:, 1] <= Y_MAX)
    return pts[inside]


def model_to_page_px(X, Y):
    """Inverse of transect_geometry._px_to_model: model metres -> pixels of the 300 dpi page render."""
    import math
    c, s_ = math.cos(tg._tilt_rad), math.sin(tg._tilt_rad)
    xr, yr = np.asarray(X, float) - tg.MAP_OFFSET[0], np.asarray(Y, float) - tg.MAP_OFFSET[1]
    xm, ym = xr * c + yr * s_, -xr * s_ + yr * c
    return xm * tg.MAP_SCALE_PX_PER_M + tg.MAP_ORIGIN_PX[0], ym * tg.MAP_SCALE_PX_PER_M + tg.MAP_ORIGIN_PX[1]


def draw_t5_on_map(ax_a, orient):
    """T5 on the printed map: model line -> page pixels -> crop pixels -> the panel's orientation."""
    x0, y0, x1, y1 = CROP
    W, H = x1 - x0, y1 - y0
    us = np.linspace(0.0, tg.transect_frame(5)['tape_length'], 600)
    pts = np.array([tg.transect_position(5, u) for u in us])
    px, py = model_to_page_px(pts[:, 0], pts[:, 1])
    cx, cy = px - x0, py - y0
    keep = (cx >= 0) & (cx <= W - 1) & (cy >= 0) & (cy <= H - 1)
    cx, cy = cx[keep], cy[keep]
    if orient == 'flipud':
        cy = H - 1 - cy
    elif orient == 'rot180':
        cx, cy = W - 1 - cx, H - 1 - cy
    ax_a.plot(cx, cy, **T5_STYLE)
    # label at the landward end of the drawn part, nudged inward
    # label at the seaward end of the drawn part (the water side), nudged inside the frame
    import matplotlib.patheffects as pe
    dy = -34 if orient != 'none' else 34
    ax_a.text(cx[0] - 14, cy[0] + dy, 'T5', fontsize=8, color=T5_COLOR, ha='right', va='center', zorder=5,
              path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])


def draw_stand(ax, stand, meta, holdout=(), geom_dir=None):
    """Right panel in the style of validate.fig_map_model, plus the shaded clipped strip."""
    # schematic: the treeless band of the substrate map (seaward of the forest edge, following the mapped
    # sand-bank boundary west of T3 so that it narrows towards T4); no synthetic tree stands there; the small water / terrestrial
    # patches inside the forest are masked in the stand as well but not drawn (user's choice)
    d = np.load(SUBSTRATE_MAP_FILE)
    band = d['treeless']
    sand_img = ax.imshow(np.ma.masked_where(~band, np.ones_like(band, dtype=float)), origin='lower',
                         extent=(X_MIN, X_MAX, Y_MIN, Y_MAX), cmap=ListedColormap([SAND_FILL]),
                         vmin=0.5, vmax=1.5, interpolation='nearest', zorder=SAND_ZORDER)
    sand_img.set_zorder(SAND_ZORDER)
    polys = tg.build_transect_polygons(tg.ALL_TRANSECT_IDS)
    from tidal_dispersal_model import MATURITY
    repro = pd.Series([r >= MATURITY.get(sp, 1.0) for r, sp in zip(stand['r_stem'], stand['species'])],
                      index=stand.index)
    for src, alpha, z in (('mtppr', 0.6, 1), ('infill', 0.7, 2), ('anchor', 1.0, 3)):
        d = stand[stand['source'] == src]
        for sp in C.SPECIES_LEVELS:
            for is_parent in (False, True):
                dd = d[(d['species'] == sp) & (repro[d.index] == is_parent)]
                if not len(dd):
                    continue
                size = 10 + 320 * dd['r_stem']
                if is_parent:
                    ax.scatter(dd['x'], dd['y'], s=size, c=SP_COLORS[sp], alpha=alpha,
                               edgecolors='k' if src == 'anchor' else 'none', linewidths=0.5, zorder=z + 3,
                               label=f'{sp} parent tree' if src == 'mtppr' else None)
                else:
                    ax.scatter(dd['x'], dd['y'], s=size, facecolors='none', edgecolors=SP_COLORS[sp], alpha=0.8,
                               linewidths=0.7, zorder=z)
    n_par = int(repro.sum())
    for t, p in polys.items():
        cs = np.array(p['corners'] + [p['corners'][0]])
        held = t in holdout
        ax.plot(cs[:, 0], cs[:, 1], color=HOLDOUT_COLOR if held else 'k', ls='--' if held else '-',
                lw=1.4 if held else 0.8, zorder=4)
        ax.text(p['end'][0], min(p['end'][1] + 4, Y_MAX + 2), f'T{t}' + (' (withheld)' if held else ''),
                ha='center', va='bottom', fontsize=8, color=HOLDOUT_COLOR if held else 'k', zorder=5)
    # the surveyed envelopes (dashed seaward fringe / dotted landward ends, interpolated between the
    # transects) were drawn here; removed 2026-09-21 at the user's request ('remove the two lines that
    # are connecting all transects, they serve no purpose'). The seaward edge of the sand band still
    # shows the fringe line, and the envelopes remain the clip used in postprocess.py.
    ax.plot([X_MIN, X_MAX, X_MAX, X_MIN, X_MIN], [Y_MIN, Y_MIN, Y_MAX, Y_MAX, Y_MIN], 'k-', lw=1, zorder=4)
    for (x0, y0), (x1, y1) in zone_boundary_lines():
        # thin faint dash-dot, under the trees: emphasising them (0.25 / lw 1.3 / zorder 3.5) was
        # tried on 2026-09-21 and reverted at the user's request ('I think it looked better before').
        ax.plot([x0, x1], [y0, y1], color='0.5', lw=0.6, ls='-.', zorder=1)
    ax.set_xlim(X_MIN - 5, X_MAX + 5)
    ax.set_ylim(Y_MIN - 5, Y_MAX + 12)     # model orientation (user's choice): sea at the bottom
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    n = len(stand)
    n_meas = int((stand['source'] == 'anchor').sum())
    if holdout:
        kept = [t for t in tg.PARALLEL_IDS if t not in holdout]
        ax.set_title(f'{n} trees: {n_meas} measured on T{", T".join(map(str, kept))}, {n - n_meas} reconstructed; '
                     f'T{" and T".join(map(str, sorted(holdout)))} withheld', fontsize=9)
    else:
        ax.set_title(f'{n} trees: {n_meas} measured on the transects, {n - n_meas} reconstructed', fontsize=9)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(plt.Line2D([], [], marker='o', ls='none', mfc='none', mec='0.4', ms=5))
    labels.append('below maturity threshold (no seeds)')
    handles.append(plt.Rectangle((0, 0), 1, 1, fc=SAND_FILL, ec='k', lw=0.4))
    labels.append('sand bank')
    if holdout:
        handles.append(plt.Line2D([], [], color='k', ls='-', lw=0.8))
        kept_ids = [t for t in tg.PARALLEL_IDS if t not in holdout]
        labels.append('transect used for the reconstruction (' + ', '.join(f'T{t}' for t in kept_ids) + ')')
        handles.append(plt.Line2D([], [], color=HOLDOUT_COLOR, ls='--', lw=1.4))
        labels.append('transect withheld from the reconstruction (' + ', '.join(f'T{t}' for t in sorted(holdout)) + ': adult evaluation)')
    ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=2, fontsize=7.0,
              markerscale=1.0, title='filled = parent tree, black edge = measured', title_fontsize=7.5,
              frameon=False, columnspacing=1.2, handletextpad=0.6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--page', default=DEFAULT_PAGE)
    ap.add_argument('--stand', default=DEFAULT_STAND)
    ap.add_argument('--geom', default=DEFAULT_GEOM)
    ap.add_argument('--out', default=C.FIG_DIR)
    ap.add_argument('--name', default='fig_layout_stand.png')
    ap.add_argument('--holdout', type=int, nargs='*', default=DEFAULT_HOLDOUT, help='transect ids withheld from the reconstruction (drawn red dashed)')
    ap.add_argument('--paper-orient', choices=['flipud', 'rot180', 'none'], default='flipud',
                    help="orientation of the printed map in panel A: flipud (default) = mirrored "
                         "top-bottom, i.e. the model frame's own transform (sea at the bottom, same "
                         "left-to-right transect order as the model panel), rot180 = turned upside "
                         "down (sea at the bottom but alongshore order reversed), none = as printed "
                         "(sea at the top)")
    args = ap.parse_args()
    if args.stand is None or args.geom is None:
        import json
        mfile = os.path.join(STAND3T_DIR, 'stand_meta.json')
        member = 1
        if os.path.exists(mfile):
            with open(mfile) as f:
                member = int(json.load(f).get('representative_member', 1))
        args.stand = args.stand or os.path.join(STAND3T_DIR, f'reconstructed_stand_m{member}.csv')
        args.geom = args.geom or os.path.join(STAND3T_DIR, f'geom_m{member}')
    stand = pd.read_csv(args.stand)
    meta = load_geom(args.geom)['meta']
    print('stand', args.stand, '| geom', args.geom)

    fig, axes = plt.subplots(1, 2, figsize=(14.6, 6.8), gridspec_kw=dict(width_ratios=[1.0, 1.15], wspace=0.30))
    ax_a, ax = axes
    panel_a = load_panel_a(args.page)
    if args.paper_orient == 'rot180':
        panel_a = panel_a[::-1, ::-1]      # proper 180 deg rotation: sea at the bottom, T6 on the left
    elif args.paper_orient == 'flipud':
        panel_a = panel_a[::-1]            # mirror top-bottom: the transform the model frame itself uses
    ax_a.imshow(panel_a)
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    for sp in ax_a.spines.values():
        sp.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PAPER_COLORS[k], ec='k', lw=0.4) for k, _ in LEGEND_ENTRIES]
    ax_a.legend(handles, [lab for _, lab in LEGEND_ENTRIES], loc='upper center', bbox_to_anchor=(0.5, -0.02),
                ncol=3, fontsize=7.5, frameon=False, handlelength=1.2, handleheight=1.2, columnspacing=1.2,
                handletextpad=0.5)

    draw_stand(ax, stand, meta, holdout=set(args.holdout or []), geom_dir=args.geom if args.holdout else None)
    fig.tight_layout()
    out = os.path.join(args.out, args.name)
    fig.savefig(out, dpi=300, bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main()
