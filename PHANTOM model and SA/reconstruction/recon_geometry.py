"""Geometry for the reconstruction: rotated frame, belts, corridors, bands, tiles.

Transect frames come from `transect_geometry.py` (map-digitised layout).  The
reconstruction works in a frame rotated so that the five parallel transects are
vertical: xp = cross-transect (west -> east), yp = along transect (seaward -> landward).
T5 is oblique in this frame; every belt segment is therefore a rotated rectangle
(p0, unit along u, along-range [s0, s1), half width hw), and all membership tests use
that representation.  Tiles (corridor x band of the parallel transects) stay
axis-aligned, as MTPPR requires for its reconstruction window.
"""
import math

import numpy as np
import pandas as pd

import recon_config as C
import transect_geometry as tg
from tidal_dispersal_model import X_MIN, X_MAX, Y_MIN, Y_MAX


# --- Frames ---------------------------------------------------------------------------

def to_rot(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return x * tg.PERP_DX + y * tg.PERP_DY, x * tg.ALONG_DX + y * tg.ALONG_DY


def to_model(xp, yp):
    xp = np.asarray(xp, dtype=float)
    yp = np.asarray(yp, dtype=float)
    return xp * tg.PERP_DX + yp * tg.ALONG_DX, xp * tg.PERP_DY + yp * tg.ALONG_DY


def domain_bbox_rot():
    cx = np.array([X_MIN, X_MAX, X_MAX, X_MIN])
    cy = np.array([Y_MIN, Y_MIN, Y_MAX, Y_MAX])
    xp, yp = to_rot(cx, cy)
    return float(xp.min()), float(xp.max()), float(yp.min()), float(yp.max())


def in_domain(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return (x >= X_MIN) & (x <= X_MAX) & (y >= Y_MIN) & (y <= Y_MAX)


def transect_frames(ids=None):
    """Per transect (rotated frame): u=0 point (x0, y0), along unit (ux, uy), right-hand
    perp (nx, ny), map-per-tape factor f, plot count, tape length, tape offset."""
    frames = {}
    for t in (ids or C.STAND_TRANSECTS):
        fr = tg.transect_frame(t)
        x0, y0 = to_rot(*fr['start'])
        ux, uy = to_rot(*fr['along'])
        nx, ny = uy, -ux
        frames[t] = dict(transect=t, x0=float(x0), y0=float(y0), ux=float(ux), uy=float(uy),
                         nx=float(nx), ny=float(ny), f=float(fr['f']), n_plots=int(fr['n_plots']),
                         tape_length=float(fr['tape_length']), tape_offset=float(fr['tape_offset']),
                         parallel=bool(t in tg.PARALLEL_IDS), L=float(fr['tape_length']))
    return frames


def parallel_ids(frames):
    return sorted([t for t in frames if frames[t]['parallel']], key=lambda t: frames[t]['x0'])


def corridor_bounds(frames):
    """Corridor xp bounds for the parallel transects: midlines between neighbours, bbox outside."""
    xp_min, xp_max, _, _ = domain_bbox_rot()
    order = parallel_ids(frames)
    bounds = {}
    for i, t in enumerate(order):
        lo = xp_min if i == 0 else 0.5 * (frames[order[i - 1]]['x0'] + frames[t]['x0'])
        hi = xp_max if i == len(order) - 1 else 0.5 * (frames[t]['x0'] + frames[order[i + 1]]['x0'])
        bounds[t] = (float(lo), float(hi))
    return bounds


# --- Bands -----------------------------------------------------------------------------

def band_segments(L, u_start=0.0, u_min=None, u_max=None, breaks=None, min_remnant=None):
    """Band segments [(k, u0, u1)] in tape metres for a belt sampled over [u_start, L].

    Breaks (fringe distance u, last band open to L) default to C.BAND_BREAKS; a remnant
    shorter than min_remnant (default C.MIN_REMNANT) is merged into the previous band.
    u_min/u_max extend the first/last segment (anchors overshooting the belt ends).
    2026-09-20 B1: `breaks` and `min_remnant` are explicit arguments so that the
    three-transect design can pass its own layout (persisted in prep_meta).
    """
    breaks = [float(b) for b in (C.BAND_BREAKS if breaks is None else breaks)]
    min_remnant = C.MIN_REMNANT if min_remnant is None else float(min_remnant)
    segs = []
    for k, b0 in enumerate(breaks):
        b1 = breaks[k + 1] if k + 1 < len(breaks) else L
        u0, u1 = max(b0, u_start), min(b1, L)
        if u1 - u0 <= 1e-9:
            continue
        segs.append([k, float(u0), float(u1)])
    if len(segs) > 1 and (segs[-1][2] - segs[-1][1]) < min_remnant:
        segs[-2][2] = segs[-1][2]
        segs.pop()
    if len(segs) > 1 and (segs[0][2] - segs[0][1]) < min_remnant:
        segs[1][1] = segs[0][1]
        segs.pop(0)
    if u_min is not None:
        segs[0][1] = min(segs[0][1], float(u_min))
    if u_max is not None:
        segs[-1][2] = max(segs[-1][2], float(u_max))
    return [tuple(s) for s in segs]


def band_of_u(segments, u):
    u = np.asarray(u, dtype=float)
    out = np.full(u.shape, -1, dtype=int)
    for k, u0, u1 in segments:
        out[(u >= u0) & (u < u1)] = k
    return out


# --- Positions along transects ------------------------------------------------------------

def plot_u(t, plot):
    return tg.plot_u(t, plot)


def position(frame, u):
    """Rotated-frame position of fringe distance u (tape m) on a transect."""
    u = np.asarray(u, dtype=float)
    return frame['x0'] + u * frame['f'] * frame['ux'], frame['y0'] + u * frame['f'] * frame['uy']


def local_coords(frame, xp, yp):
    """(along map m from u=0, signed perp m) of points in a transect's frame."""
    dx = np.asarray(xp, dtype=float) - frame['x0']
    dy = np.asarray(yp, dtype=float) - frame['y0']
    return dx * frame['ux'] + dy * frame['uy'], dx * frame['nx'] + dy * frame['ny']


def pcqm_points(frames):
    rows = []
    for t, fr in frames.items():
        for k in range(1, fr['n_plots'] + 1):
            u = float(plot_u(t, k))
            px, py = position(fr, u)
            rows.append(dict(transect=t, plot=k, u=u, px=float(px), py=float(py),
                             ux=fr['ux'], uy=fr['uy'], nx=fr['nx'], ny=fr['ny']))
    return pd.DataFrame(rows)


# --- Sectors ---------------------------------------------------------------------------------

def in_sector(xp, yp, sectors, exclude=None):
    """True for points inside any quarter sector.

    sectors: DataFrame(px, py, ux, uy, sx, sy, radius): px/py = PCQM point, (ux, uy) =
    transect along unit, sx/sy = quarter signs on the perp/along axes.
    """
    xp = np.atleast_1d(np.asarray(xp, dtype=float))
    yp = np.atleast_1d(np.asarray(yp, dtype=float))
    if len(sectors) == 0:
        return np.zeros(xp.shape, dtype=bool)
    keep = np.ones(len(sectors), dtype=bool) if exclude is None else ~np.asarray(exclude, dtype=bool)
    s = sectors[keep]
    dx = xp[:, None] - s['px'].values[None, :]
    dy = yp[:, None] - s['py'].values[None, :]
    ux, uy = s['ux'].values[None, :], s['uy'].values[None, :]
    along = dx * ux + dy * uy
    perp = dx * uy - dy * ux
    hit = (s['sx'].values[None, :] * perp > 0) & (s['sy'].values[None, :] * along > 0) & \
          (dx * dx + dy * dy < (s['radius'].values ** 2)[None, :])
    return hit.any(axis=1)


# --- Rotated rectangles (belt segments) ------------------------------------------------------

def rect_local(rects, xp, yp):
    """Along/perp coordinates of points relative to every rect: arrays (n_points, n_rects)."""
    xp = np.atleast_1d(np.asarray(xp, dtype=float))
    yp = np.atleast_1d(np.asarray(yp, dtype=float))
    dx = xp[:, None] - rects['x0'].values[None, :]
    dy = yp[:, None] - rects['y0'].values[None, :]
    ux, uy = rects['ux'].values[None, :], rects['uy'].values[None, :]
    return dx * ux + dy * uy, dx * uy - dy * ux


def _clip_cols(rects):
    n = len(rects)
    out = []
    for c, d in (('cx0', -np.inf), ('cx1', np.inf), ('cy0', -np.inf), ('cy1', np.inf)):
        out.append(rects[c].values.astype(float) if c in rects.columns else np.full(n, d))
    return out


def _hits(rects, xp, yp):
    xp = np.atleast_1d(np.asarray(xp, dtype=float))
    yp = np.atleast_1d(np.asarray(yp, dtype=float))
    along, perp = rect_local(rects, xp, yp)
    cx0, cx1, cy0, cy1 = _clip_cols(rects)
    return ((along >= rects['s0'].values[None, :]) & (along < rects['s1'].values[None, :]) &
            (np.abs(perp) < rects['hw'].values[None, :]) &
            (xp[:, None] >= cx0[None, :]) & (xp[:, None] < cx1[None, :]) &
            (yp[:, None] >= cy0[None, :]) & (yp[:, None] < cy1[None, :]))


def in_rects(xp, yp, rects):
    """Membership in any (optionally box-clipped) rotated rectangle."""
    if len(rects) == 0:
        return np.zeros(np.atleast_1d(xp).shape, dtype=bool)
    return _hits(rects, xp, yp).any(axis=1)


def rect_index(xp, yp, rects):
    """Index of the first rect containing each point, -1 if none."""
    if len(rects) == 0:
        return np.full(np.atleast_1d(xp).shape, -1, dtype=int)
    hit = _hits(rects, xp, yp)
    out = np.full(hit.shape[0], -1, dtype=int)
    any_hit = hit.any(axis=1)
    out[any_hit] = hit[any_hit].argmax(axis=1)
    return out


def rect_corners(r):
    ux, uy = r.ux, r.uy
    nx, ny = uy, -ux
    pts = []
    for s, w in ((r.s0, -r.hw), (r.s0, r.hw), (r.s1, r.hw), (r.s1, -r.hw)):
        pts.append((r.x0 + s * ux + w * nx, r.y0 + s * uy + w * ny))
    return pts


def sample_in_rect(r, rng, n):
    """Uniform samples inside a (clipped) rect; returns xp, yp and the acceptance fraction."""
    ms = rng.uniform(r.s0, r.s1, n)
    mw = rng.uniform(-r.hw, r.hw, n)
    xp = r.x0 + ms * r.ux + mw * r.uy
    yp = r.y0 + ms * r.uy - mw * r.ux
    cx0 = getattr(r, 'cx0', -np.inf)
    cx1 = getattr(r, 'cx1', np.inf)
    cy0 = getattr(r, 'cy0', -np.inf)
    cy1 = getattr(r, 'cy1', np.inf)
    keep = (xp >= cx0) & (xp < cx1) & (yp >= cy0) & (yp < cy1)
    return xp[keep], yp[keep], float(keep.mean())


def rect_area(r, rng, n=20000):
    _, _, frac = sample_in_rect(r, rng, n)
    return float((r.s1 - r.s0) * 2 * r.hw * frac)


def tile_of(xp, yp, tiles):
    xp = np.atleast_1d(np.asarray(xp, dtype=float))
    yp = np.atleast_1d(np.asarray(yp, dtype=float))
    out = np.full(len(xp), -1, dtype=int)
    for tl in tiles.itertuples(index=False):
        m = in_box(xp, yp, tl) & (out < 0)
        out[m] = tl.tile_id
    return out


def band_of_tile(tile_ids, tiles):
    lut = dict(zip(tiles['tile_id'], tiles['band']))
    return np.array([lut.get(int(t), -1) for t in np.atleast_1d(tile_ids)], dtype=int)


def build_core_rects(frames, segments, tiles, rng):
    """Belt segments as rotated rectangles.

    Parallel transects: one rect per band segment (its own corridor's tile rows).
    Oblique transects (T5): the whole belt clipped to every tile box it intersects, one
    piece per tile, banded by that tile (rect_id = 5000 + tile_id).
    """
    rows = []
    for t, fr in frames.items():
        if fr['parallel']:
            for k, u0, u1 in segments[t]:
                rows.append(dict(rect_id=100 * t + k, tile_id=100 * t + k, transect=t, band=k,
                                 x0=fr['x0'], y0=fr['y0'], ux=fr['ux'], uy=fr['uy'],
                                 s0=u0 * fr['f'], s1=u1 * fr['f'], hw=C.BELT_HALF, u0=u0, u1=u1,
                                 cx0=-np.inf, cx1=np.inf, cy0=-np.inf, cy1=np.inf,
                                 area=(u1 - u0) * fr['f'] * 2 * C.BELT_HALF, f=fr['f'], clipped=0))
        else:
            u0, u1 = segments[t][0][1], segments[t][-1][2]
            full = dict(x0=fr['x0'], y0=fr['y0'], ux=fr['ux'], uy=fr['uy'], s0=u0 * fr['f'], s1=u1 * fr['f'], hw=C.BELT_HALF)
            for tl in tiles.itertuples(index=False):
                r = pd.Series(dict(full, cx0=tl.xp_min, cx1=tl.xp_max, cy0=tl.yp_min, cy1=tl.yp_max))
                area = rect_area(r, rng, 40000)
                if area < 25.0:
                    continue
                rows.append(dict(rect_id=5000 + int(tl.tile_id), tile_id=int(tl.tile_id), transect=t, band=int(tl.band),
                                 x0=fr['x0'], y0=fr['y0'], ux=fr['ux'], uy=fr['uy'], s0=full['s0'], s1=full['s1'],
                                 hw=C.BELT_HALF, u0=u0, u1=u1, cx0=tl.xp_min, cx1=tl.xp_max, cy0=tl.yp_min, cy1=tl.yp_max,
                                 area=area, f=fr['f'], clipped=1))
    return pd.DataFrame(rows)


def build_tiles(frames, segments, corridors, fill_domain=True):
    """Corridor x band tiles.  With fill_domain the first band of every corridor extends
    to the seaward domain edge and the last band to the landward edge, so the tiles cover
    the whole box (the extensions are extrapolations of the fringe / landward bands)."""
    _, _, yb0, yb1 = domain_bbox_rot()
    rows = []
    for t in parallel_ids(frames):
        fr = frames[t]
        lo, hi = corridors[t]
        segs = segments[t]
        for i, (k, u0, u1) in enumerate(segs):
            y0 = fr['y0'] + u0
            y1 = fr['y0'] + u1
            if fill_domain and i == 0:
                y0 = min(y0, yb0)
            if fill_domain and i == len(segs) - 1:
                y1 = max(y1, yb1)
            rows.append(dict(tile_id=100 * t + k, transect=t, band=k, xp_min=lo, xp_max=hi,
                             yp_min=y0, yp_max=y1, u0=y0 - fr['y0'], u1=y1 - fr['y0'],
                             area=(hi - lo) * (y1 - y0)))
    return pd.DataFrame(rows)


def in_box(xp, yp, tile):
    xp = np.atleast_1d(np.asarray(xp, dtype=float))
    yp = np.atleast_1d(np.asarray(yp, dtype=float))
    return (xp >= tile.xp_min) & (xp < tile.xp_max) & (yp >= tile.yp_min) & (yp < tile.yp_max)


def tile_eligible_areas(tiles, npz_path=None):
    """2026-09-20 B1: eligible area of every tile = tile area minus the treeless cells of the
    habitat mask (grid 'treeless' of substrate_map.npz, model frame, 1 m cells).  Cell centres
    are converted to the rotated frame and counted inside each tile box.

    Returns a DataFrame (tile_id, n_cells, n_treeless, treeless_area, eligible_area)."""
    d = np.load(npz_path or C.SUBSTRATE_NPZ)
    mask = d['treeless'].astype(bool)
    cell = float(d['cell'])
    x0, y0 = float(d['x_min']), float(d['y_min'])
    ny, nx = mask.shape
    cx = x0 + (np.arange(nx) + 0.5) * cell
    cy = y0 + (np.arange(ny) + 0.5) * cell
    X, Y = np.meshgrid(cx, cy)
    xp, yp = to_rot(X.ravel(), Y.ravel())
    tl = mask.ravel()
    rows = []
    for t in tiles.itertuples(index=False):
        inb = in_box(xp, yp, t)
        n_cells = int(inb.sum())
        n_tl = int((inb & tl).sum())
        rows.append(dict(tile_id=int(t.tile_id), n_cells=n_cells, n_treeless=n_tl,
                         treeless_area=n_tl * cell * cell, eligible_area=(n_cells - n_tl) * cell * cell))
    return pd.DataFrame(rows)


def rect_touches_box(r, tile):
    if getattr(r, 'clipped', 0):
        return int(r.tile_id) == int(tile.tile_id)
    cs = np.array(rect_corners(r))
    if in_box(cs[:, 0], cs[:, 1], tile).any():
        return True
    s = np.linspace(r.s0, r.s1, 60)
    return in_box(r.x0 + s * r.ux, r.y0 + s * r.uy, tile).any()


def union_area(rects, rng, n=60000):
    """Monte Carlo area of the union of rects (handles overlaps at crossings)."""
    if len(rects) == 0:
        return 0.0
    cs = np.vstack([rect_corners(r) for r in rects.itertuples(index=False)])
    x0, x1 = cs[:, 0].min(), cs[:, 0].max()
    y0, y1 = cs[:, 1].min(), cs[:, 1].max()
    mx = rng.uniform(x0, x1, n)
    my = rng.uniform(y0, y1, n)
    return float(in_rects(mx, my, rects).mean() * (x1 - x0) * (y1 - y0))


# --- Envelope (fringe and landward limit interpolated between parallel belts) ---------------

def envelope(frames, xp, segments=None):
    """(yp_fringe, yp_land) at cross positions xp from the parallel transects.

    Constant across each belt (x0 +/- BELT_HALF), piecewise-linear across the gaps,
    constant beyond the outer belts.  With `segments` the belt ends follow the
    (possibly extended) first/last segment.
    """
    order = parallel_ids(frames)
    xs, yf, yl = [], [], []
    for t in order:
        fr = frames[t]
        if segments is None:
            f, l = fr['y0'] + fr['tape_offset'], fr['y0'] + fr['tape_length']
        else:
            f, l = fr['y0'] + segments[t][0][1], fr['y0'] + segments[t][-1][2]
        for dx in (-C.BELT_HALF, C.BELT_HALF):
            xs.append(fr['x0'] + dx)
            yf.append(f)
            yl.append(l)
    xp = np.asarray(xp, dtype=float)
    return np.interp(xp, np.array(xs), np.array(yf)), np.interp(xp, np.array(xs), np.array(yl))


def corridor_of(xp, corridors):
    xp = np.atleast_1d(np.asarray(xp, dtype=float))
    out = np.zeros(xp.shape, dtype=int)
    for t, (lo, hi) in corridors.items():
        out[(xp >= lo) & (xp < hi)] = t
    return out


def fringe_distance(xp, yp, frames, corridors):
    """u relative to the corridor transect's u=0 point (parallel transects)."""
    t = corridor_of(xp, corridors)
    y0 = np.array([frames[tt]['y0'] if tt in frames else np.nan for tt in t])
    return np.asarray(yp, dtype=float) - y0, t


def geometry_check():
    frames = transect_frames()
    for t, fr in frames.items():
        x, y = to_model(fr['x0'], fr['y0'])
        sx, sy = tg.transect_frame(t)['start']
        assert abs(x - sx) < 1e-9 and abs(y - sy) < 1e-9
    ax, ay = to_rot(tg.ALONG_DX, tg.ALONG_DY)
    px, py = to_rot(tg.PERP_DX, tg.PERP_DY)
    assert abs(ax) < 1e-12 and abs(ay - 1) < 1e-12 and abs(px - 1) < 1e-12 and abs(py) < 1e-12
    return frames


if __name__ == '__main__':
    fr = geometry_check()
    cb = corridor_bounds(fr)
    print('rotated domain bbox:', [round(v, 1) for v in domain_bbox_rot()])
    for t in parallel_ids(fr) + [5]:
        f = fr[t]
        segs = band_segments(f['tape_length'], u_start=f['tape_offset'])
        print(f"T{t}: x0={f['x0']:.1f} y0={f['y0']:.1f} u=({f['ux']:.3f},{f['uy']:.3f}) f={f['f']:.3f} "
              f"plots={f['n_plots']} tape={f['tape_length']:.0f} offset={f['tape_offset']:.0f} "
              f"corridor={tuple(round(v, 1) for v in cb[t]) if t in cb else '-'} bands={segs}")
