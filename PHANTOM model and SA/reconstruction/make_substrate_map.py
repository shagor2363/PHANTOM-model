"""Substrate map of the model frame from the vegetation map of Dahdouh-Guebas et al. (2004, Fig. 6a).

Every pixel of the 300-dpi page render inside the map frame is classified by the nearest
colour of the paper's legend, transformed into the model frame with the same digitisation
(transect_geometry._px_to_model) and binned into 1 m cells by majority vote. White pixels
(quadrat symbols, seaward margin) carry no class and are filled from the nearest classified
cell; patches smaller than MIN_PATCH_CELLS are speckle and take the nearest larger patch's class.
Classes are grouped into physiographic substrate types used by the establishment
submodel:

    0 mangrove         Avicennia, Ceriops, Rhizophora, Sonneratia, Bruguiera dominated, mixed mangrove
    1 sandy stretch    sand within the forest
    2 sand bank        intertidal sand bank seaward of the fringe
    3 water            water and seagrass beds
    4 terrestrial      terrestrial vegetation

Outputs: results/reconstruction/substrate_map.npz: 'group' (MODEL INPUT: 0 forest floor, 1 sand band =
         'treeless', everything seaward of the forest edge; no rooting and no parent trees there),
         'paper_group' / 'cls' (the paper's classes, figures only), 'in_frame', origin, cell size, names
         results/reconstruction/figures/fig_substrate_map.png

Usage: python reconstruction/make_substrate_map.py [--page PAGE] [--cell 1.0]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402
from scipy import ndimage  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recon_config as C  # noqa: E402
import transect_geometry as tg  # noqa: E402
from tidal_dispersal_model import X_MIN, X_MAX, Y_MIN, Y_MAX, zone_boundary_lines  # noqa: E402

# Colours sampled from the legend of Dahdouh-Guebas et al. (2004), Fig. 6
PAPER_COLORS = {
    'empty quadrat': '#ffffff',
    'Avicennia marina': '#a7a4b7',
    'Bruguiera gymnorrhiza': '#017898',
    'Ceriops tagal': '#d08c34',
    'Rhizophora mucronata': '#458149',
    'Sonneratia alba': '#d2dc34',
    'Mixed mangrove': '#f46843',
    'Terrestrial vegetation': '#d9a993',
    'Sand (incl. sandy stretch)': '#fdf109',
    'Sand bank': '#fffac3',
    'Water + seagrasses': '#90c7eb',
}
SPECIES_KEY = {'Avicennia': 'Avicennia marina', 'Bruguiera': 'Bruguiera gymnorrhiza',
               'Ceriops': 'Ceriops tagal', 'Rhizophora': 'Rhizophora mucronata',
               'Sonneratia': 'Sonneratia alba'}
# legend entries in the order of the published legend
LEGEND_ENTRIES = [
    ('empty quadrat', 'empty quadrat'),
    ('Avicennia marina', r'$\it{Avicennia\ marina}$ dominated'),
    ('Bruguiera gymnorrhiza', r'$\it{Bruguiera\ gymnorrhiza}$ dominated'),
    ('Ceriops tagal', r'$\it{Ceriops\ tagal}$ dominated'),
    ('Rhizophora mucronata', r'$\it{Rhizophora\ mucronata}$ dominated'),
    ('Sonneratia alba', r'$\it{Sonneratia\ alba}$ dominated'),
    ('Mixed mangrove', 'Mixed mangrove'),
    ('Terrestrial vegetation', 'Terrestrial vegetation'),
    ('Sand (incl. sandy stretch)', 'Sand (incl. sandy stretch)'),
    ('Sand bank', 'Sand bank'),
    ('Water + seagrasses', 'Water + seagrasses'),
]
CROP = (236, 349, 1458, 1526)   # black frame of the published panel in the 300 dpi render
LABEL_BOX = (16, 58, 75, 106)   # the printed '(a)' inside that crop


def load_panel_a(page):
    img = np.asarray(Image.open(page).convert('RGB').crop(CROP)).astype(np.uint8).copy()
    x0, y0, x1, y1 = LABEL_BOX
    ring = np.concatenate([img[max(0, y0 - 12):y0, x0:x1].reshape(-1, 3),
                           img[y1:y1 + 12, x0:x1].reshape(-1, 3),
                           img[y0:y1, max(0, x0 - 12):x0].reshape(-1, 3),
                           img[y0:y1, x1:x1 + 12].reshape(-1, 3)])
    img[y0:y1, x0:x1] = np.median(ring, axis=0).astype(np.uint8)
    return img


GROUP_NAMES = ['mangrove', 'sandy stretch', 'sand bank', 'water + seagrass', 'terrestrial']   # paper classes grouped (figure only)
MODEL_NAMES = ['forest', 'sand']          # what the model reads: forest floor / sand band (no rooting)
MODEL_COLORS = ['#458149', '#fdf6a8']
MIN_PATCH_CELLS = 25   # patches smaller than this (m2) are classification speckle, not map polygons
GROUP_COLORS = ['#458149', '#fdf109', '#fffac3', '#90c7eb', '#d9a993']
CLASS_TO_GROUP = {
    'Avicennia marina': 0, 'Bruguiera gymnorrhiza': 0, 'Ceriops tagal': 0, 'Rhizophora mucronata': 0,
    'Sonneratia alba': 0, 'Mixed mangrove': 0, 'Sand (incl. sandy stretch)': 1, 'Sand bank': 2,
    'Water + seagrasses': 3, 'Terrestrial vegetation': 4, 'empty quadrat': -1,
}
DEFAULT_PAGE = os.path.join(C.FIG_DIR, 'dahdouh_guebas_2004_page09_300dpi.png')
OUT_NPZ = os.path.join(C.OUT_ROOT, 'substrate_map.npz')


USER_EDGE_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'forest_edge_user.csv')


def transect_ends_line(xs):
    """Piecewise-linear line through the seaward ends of the parallel transects, T4 at the start of
    its surveyed part (tape 120 m), constant beyond T4 and T6."""
    pts = []
    for t in tg.ALL_TRANSECT_IDS:
        fr = tg.transect_frame(t)
        if not fr.get('parallel', True):
            continue
        pts.append((fr['start'][0] + fr['tape_offset'] * fr['along'][0],
                    fr['start'][1] + fr['tape_offset'] * fr['along'][1]))
    pts.sort()
    return np.interp(np.asarray(xs, dtype=float), [q[0] for q in pts], [q[1] for q in pts])


def forest_fringe_line(xs):
    """Seaward forest edge y(x) in the model frame. Seaward of it no adult of the modelled species
    stands (the 'treeless' band, clipped from the stand by postprocess.py --clip-treeless and shaded
    as sand in the schematic maps). West of T1 the edge is the line drawn by hand on the stand map
    (reconstruction/forest_edge_user.csv, 2026-09-06): it leaves T4's surveyed start hugging the
    left edge of the box and bends seaward to run 10-15 m above the seaward boundary, i.e. the
    sand band narrows towards T4 as on the printed map; east of its last knot, and everywhere when
    the file is missing, the line through the transect seaward ends applies."""
    xs = np.asarray(xs, dtype=float)
    base = transect_ends_line(xs)
    if not os.path.exists(USER_EDGE_CSV):
        return base
    k = np.loadtxt(USER_EDGE_CSV, delimiter=',')
    user = np.interp(xs, k[:, 0], k[:, 1])
    return np.where(xs <= k[-1, 0], user, base)


def classify_page(page):
    img = np.asarray(Image.open(page).convert('RGB')).astype(float)
    x0, y0, x1, y1 = CROP
    sub = img[y0:y1, x0:x1]
    names = list(PAPER_COLORS)
    cols = np.array([tuple(int(PAPER_COLORS[k].lstrip('#')[i:i + 2], 16) for i in (0, 2, 4)) for k in names], float)
    flat = sub.reshape(-1, 3)
    cls = np.empty(len(flat), dtype=np.int16)
    step = 200000
    for i in range(0, len(flat), step):
        d = ((flat[i:i + step, None, :] - cols[None, :, :]) ** 2).sum(axis=2)
        cls[i:i + step] = d.argmin(axis=1)
    py, px = np.mgrid[y0:y1, x0:x1]
    return cls, px.ravel().astype(float), py.ravel().astype(float), names


def to_model(px, py):
    """Vectorised transect_geometry._px_to_model (same affine map, checked against it)."""
    xm = (np.asarray(px, float) - tg.MAP_ORIGIN_PX[0]) / tg.MAP_SCALE_PX_PER_M
    ym = (np.asarray(py, float) - tg.MAP_ORIGIN_PX[1]) / tg.MAP_SCALE_PX_PER_M
    c, s_ = np.cos(tg._tilt_rad), np.sin(tg._tilt_rad)
    xs = xm * c - ym * s_ + tg.MAP_OFFSET[0]
    ys = xm * s_ + ym * c + tg.MAP_OFFSET[1]
    chk = tg._px_to_model((float(px[0]), float(py[0])))
    assert abs(chk[0] - xs[0]) < 1e-9 and abs(chk[1] - ys[0]) < 1e-9
    return xs, ys


def despeckle(grid, min_cells=MIN_PATCH_CELLS):
    """Reassign connected patches smaller than min_cells (colour-classification speckle from the
    map's stipple textures and polygon outlines) to the class of the nearest cell of a larger
    patch, so that the raster carries only the map's polygons."""
    keep = np.zeros(grid.shape, dtype=bool)
    for g in np.unique(grid):
        m = grid == g
        lab, k = ndimage.label(m, structure=np.ones((3, 3)))
        sizes = ndimage.sum(m, lab, range(1, k + 1))
        big = np.isin(lab, np.where(sizes >= min_cells)[0] + 1)
        keep |= big
    _, (ii, jj) = ndimage.distance_transform_edt(~keep, return_indices=True)
    out = grid[ii, jj]
    return out, int((out != grid).sum())


def build(page, cell):
    cls, px, py, names = classify_page(page)
    group = np.array([CLASS_TO_GROUP[n] for n in names], dtype=np.int16)[cls]
    mx, my = to_model(px, py)
    nx = int(np.ceil((X_MAX - X_MIN) / cell))
    ny = int(np.ceil((Y_MAX - Y_MIN) / cell))
    ix = np.floor((mx - X_MIN) / cell).astype(int)
    iy = np.floor((my - Y_MIN) / cell).astype(int)
    inside = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    ok = inside & (group >= 0)
    counts = np.zeros((ny, nx, len(GROUP_NAMES)), dtype=np.int32)
    np.add.at(counts, (iy[ok], ix[ok], group[ok]), 1)
    grid = counts.argmax(axis=2).astype(np.int16)
    empty = counts.sum(axis=2) == 0
    # the paper's own classes (majority vote) and the cells covered by the printed map frame,
    # kept for the digitisation figure (fig_layout_substrate.py); the model reads 'group' only
    ccounts = np.zeros((ny, nx, len(names)), dtype=np.int32)
    np.add.at(ccounts, (iy[ok], ix[ok], cls[ok]), 1)
    cls_grid = ccounts.argmax(axis=2).astype(np.int16)
    in_frame = np.zeros((ny, nx), dtype=bool)
    in_frame[iy[inside], ix[inside]] = True
    if empty.any():
        _, (ii, jj) = ndimage.distance_transform_edt(empty, return_indices=True)
        grid = grid[ii, jj]
        cls_grid = cls_grid[ii, jj]
    grid, n_spk = despeckle(grid)
    cls_grid, _ = despeckle(cls_grid)
    print(f'speckle removed: {n_spk} cells reassigned (patches < {MIN_PATCH_CELLS} m2)')
    # treeless band (no adult of the modelled species stands there; drawn as sand in the schematic
    # maps): everything seaward of the forest edge (forest_fringe_line). The substrate classes are
    # not changed by it (outside the printed map: nearest mapped class).
    xc = X_MIN + (np.arange(nx) + 0.5) * cell
    yc = Y_MIN + (np.arange(ny) + 0.5) * cell
    treeless = yc[:, None] < forest_fringe_line(xc)[None, :]
    print(f'treeless band: {treeless.sum() * cell ** 2 / 1e4:.2f} ha')
    return grid, empty, names, cls_grid, in_frame, treeless


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--page', default=DEFAULT_PAGE)
    ap.add_argument('--cell', type=float, default=1.0)
    ap.add_argument('--out', default=OUT_NPZ)
    ap.add_argument('--fig', default=os.path.join(C.FIG_DIR, 'fig_substrate_map.png'))
    ap.add_argument('--no-user-edge', action='store_true',
                    help='sensitivity (2026-10-11): forest edge = the line through the transect seaward ends only, no hand-drawn segment west of T1')
    ap.add_argument('--add-mapped-sand', action='store_true',
                    help="sensitivity (2026-10-11): also refuse rooting on the paper's 'Sand bank' and 'Water + seagrasses' cells")
    args = ap.parse_args()
    if args.no_user_edge:
        global USER_EDGE_CSV
        USER_EDGE_CSV = ''
    grid, empty, names, cls_grid, in_frame, treeless = build(args.page, args.cell)
    model_grid = treeless.astype(np.int16)          # 0 forest floor, 1 sand band (no rooting, no parent trees)
    if args.add_mapped_sand:
        mapped = np.isin(cls_grid, [names.index('Sand bank'), names.index('Water + seagrasses')])
        model_grid = (treeless | mapped).astype(np.int16)
        print(f'mapped sand bank / water added to the no-rooting class: {int(mapped.sum())} cells')
    np.savez(args.out, group=model_grid, x_min=X_MIN, y_min=Y_MIN, cell=args.cell, names=np.array(MODEL_NAMES),
             paper_group=grid, paper_group_names=np.array(GROUP_NAMES), filled_from_neighbour=empty,
             cls=cls_grid, cls_names=np.array(names), in_frame=in_frame, treeless=treeless)
    print('substrate map', grid.shape, 'cells of', args.cell, 'm; paper digitisation filled from neighbour:',
          int(empty.sum()), f'; outside the printed map frame: {100 * (~in_frame).mean():.1f}% of the box')
    from tidal_dispersal_model import zone_at
    ys = Y_MIN + (np.arange(grid.shape[0]) + 0.5) * args.cell
    xs = X_MIN + (np.arange(grid.shape[1]) + 0.5) * args.cell
    zones = np.array([[zone_at(x, y) for x in xs] for y in ys])
    print('model substrate (share of each zone):')
    print(f"{'class':<18}" + ''.join(f'{z:>9}' for z in 'CBA') + '    total')
    for g, name in enumerate(MODEL_NAMES):
        row = [100 * (model_grid[zones == z] == g).mean() for z in 'CBA']
        print(f'{name:<18}' + ''.join(f'{v:8.1f}%' for v in row) + f'{100 * (model_grid == g).mean():8.1f}%')
    band_ha = model_grid.sum() * args.cell ** 2 / 1e4
    print('sand band: %.2f ha, of which zone C %.0f%%, zone B %.0f%%' % (
        band_ha, 100 * (model_grid[zones == 'C'].sum() / model_grid.sum()),
        100 * (model_grid[zones == 'B'].sum() / model_grid.sum())))
    # figure: the model substrate
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    from matplotlib.colors import ListedColormap
    ax.imshow(model_grid, origin='lower', extent=(X_MIN, X_MAX, Y_MIN, Y_MAX), cmap=ListedColormap(MODEL_COLORS),
              vmin=-0.5, vmax=len(MODEL_NAMES) - 0.5, interpolation='nearest')
    for (x0, y0), (x1, y1) in zone_boundary_lines():
        ax.plot([x0, x1], [y0, y1], color='0.3', lw=0.8, ls='-.')
    for t, p in tg.build_transect_polygons(tg.ALL_TRANSECT_IDS).items():
        cs = np.array(p['corners'] + [p['corners'][0]])
        ax.plot(cs[:, 0], cs[:, 1], 'k-', lw=0.8)
        ax.text(0.5 * (p['corners'][0][0] + p['corners'][1][0]), Y_MAX + 4, f'T{t}', ha='center', va='bottom', fontsize=8)
    ax.contour(X_MIN + (np.arange(grid.shape[1]) + 0.5) * args.cell, Y_MIN + (np.arange(grid.shape[0]) + 0.5) * args.cell,
               in_frame.astype(float), levels=[0.5], colors='k', linewidths=0.8, linestyles='--')
    ax.contour(xs, ys, treeless.astype(float), levels=[0.5], colors='#8a7a00', linewidths=1.0)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX + 14)     # model orientation: sea at the bottom
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title('Substrate of the model frame: forest floor and the sand band seaward of the forest edge', fontsize=9)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec='k', lw=0.4) for c in MODEL_COLORS]
    handles.append(plt.Line2D([], [], color='k', ls='--', lw=0.8))
    ax.legend(handles, ['forest floor: propagules root, parent trees',
                        'sand band: no rooting (Eq. 10), no parent trees',
                        'edge of the printed map of Dahdouh-Guebas et al. (2004)'],
              loc='lower center', bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    out_fig = args.fig
    fig.savefig(out_fig, dpi=300, bbox_inches='tight')
    print('wrote', args.out, 'and', out_fig)


if __name__ == '__main__':
    main()
