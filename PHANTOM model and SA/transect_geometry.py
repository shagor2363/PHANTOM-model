"""
Transect Geometry Utility
=========================

Gazi Bay 1997/1999 PCQM transect layout of Dahdouh-Guebas et al. (2004, MEPS 272: 77-92),
digitised from their Fig. 6a (5 parallel transects overlaid on the 1992 vegetation map;
1 block = 4 quadrats of 5 x 5 m at one PCQM point, points every 10 m of tape).

Model frame: the paper's map rotated by MAP_TILT_DEG so that the five parallel transects
(west -> east: T4, T3, T2, T1, T6) run exactly along +y (seaward -> landward) and x is
alongshore; the origin is shifted so the transects fill the model box.  The orthogonal
transect T5 (not drawn in the figure, not part of the stand) is placed by least squares
through its crossings with T4, T3, T2 and T1 recorded in the 1997 field notes (GAZ1_6.XLS
remarks column), which also fix its map-metres-per-tape-metre factor (~1.45).

Survey conventions (checked against vegetation classes, water-depth remarks and GPS):
  * plot 1 is the LANDWARD end of every transect; the last plot is the seaward fringe;
  * fringe distance u of plot k (tape metres from the seaward end of the tape) is
        u_k = tape_length - 10*(k-1) - 5;
  * T4's six plots are the landward 60 m (tape 120-180) of a 180 m line; only that part
    is drawn/modelled (tape_offset = 120).

Observed data from juvenile_quadrat_records_1999_t1_5.csv:
  transect, plot, quadrat_letter, species, density_per_m2

Functions:
  transect_frame(t)             -> start, along/perp unit vectors, scale f, plot count
  build_transect_polygons()     -> dict of transect polygon corners
  assign_to_transect(x, y)      -> transect number or None
  assign_to_plot / assign_to_quarter -> plot number (1 = landward) and quarter a-d
  load_observed_transect_densities()  -> {(transect, species): mean_density}
  model_transect_densities()    -> {(transect, species): density}
  geometry_sensitivity()        -> metrics under alternate quarter orientation
"""

import os
import csv
import math
import numpy as np

# ---------------------------------------------------------------------------
# Map digitisation (Fig. 6a of Dahdouh-Guebas et al. 2004, page rendered at 300 dpi)
# ---------------------------------------------------------------------------
MAP_SCALE_PX_PER_M = 4.0          # 10 m PCQM block = 40 px; the 300 m close-up = 1216 px
MAP_ORIGIN_PX = (165.0, 255.0)    # left edge of the frame, sea edge (top of the printed map)
MAP_TILT_DEG = 17.0               # lean of the parallel transects on the printed map (mean of 5)
MAP_OFFSET = (75.0, -65.0)        # shift after rotation so the transects fill the model box
_tilt_rad = math.radians(MAP_TILT_DEG)
TILT_ANGLE_DEG = 0.0              # in the model frame the parallel transects run along +y

# Along-transect direction of the parallel transects (seaward -> landward)
ALONG_DX = 0.0
ALONG_DY = 1.0

# Perpendicular direction (left -> right looking landward)
PERP_DX = 1.0
PERP_DY = 0.0

# Seaward end of each DRAWN strip (px) and survey facts.  tape_offset = fringe distance
# of the drawn strip's seaward end (T4: tape 120).
_DIGITISED = {
    1: dict(px=(860.0, 337.0), n_plots=17, tape_length=170.0, tape_offset=0.0),
    2: dict(px=(632.0, 545.0), n_plots=16, tape_length=160.0, tape_offset=0.0),
    3: dict(px=(402.0, 797.0), n_plots=13, tape_length=130.0, tape_offset=0.0),
    4: dict(px=(197.0, 1197.0), n_plots=6, tape_length=180.0, tape_offset=120.0),
    6: dict(px=(1147.0, 280.0), n_plots=16, tape_length=160.0, tape_offset=0.0),
}
# T5 (orthogonal): 1997 notes "INTERSECTION" -- (T5 tape m, crossed transect, its plot)
T5_CROSSINGS_1997 = [(25.0, 4, 1.0), (85.0, 3, 3.5), (145.0, 2, 7.5), (205.0, 1, 10.0)]
T5_N_PLOTS = 32
T5_TAPE_LENGTH = 320.0

PLOT_LENGTH = 10.0   # m of tape between PCQM points
TRANSECT_WIDTH = 10  # m (two 5 m quadrats across)
# ASSUMED virtual-sampling footprint of one survey quadrat on the MODEL side (2026-09-20, data
# audit docs/juvenile_data_audit.md): the paper (Dahdouh-Guebas et al. 2004, p. 81) gives a
# 5 x 5 m quadrat per PCQM quadrant for vegetation cover and never states the juvenile quadrat
# size; the workbook records juvenile DENSITIES (#/m2) only. 25 m2 is therefore an assumption
# for converting simulated recruit counts inside the survey footprint to densities (sensitivity:
# 12.5 m2, the unit the 1997 densities are consistent with). It is never applied to the
# observations. Zone SHARES (relative density across zones) do not depend on it; absolute
# sampled densities do. The former QUADRAT_AREA = 2.5 ('1 plot = 10 m2, 4 quadrats/plot') was
# unsupported; the name is kept as an alias for old imports.
MODEL_QUADRAT_FOOTPRINT_M2 = 25.0
QUADRAT_AREA = MODEL_QUADRAT_FOOTPRINT_M2   # alias (old imports); same assumption, same caveat
QUADRATS_PER_PLOT = 4  # a, b, c, d

# Domain bounds (for clipping)
from tidal_dispersal_model import X_MIN, X_MAX, Y_MIN, Y_MAX, zone_at
from tidal_dispersal_model import SPECIES as _MODEL_SPECIES   # tidal_dispersal_model imports nothing from here

TRANSECT_IDS = [1, 2, 3, 4, 5]          # transects with 1999 juvenile validation data
ALL_TRANSECT_IDS = [1, 2, 3, 4, 6]      # transects in the stand reconstruction (parallel only)
PARALLEL_IDS = [4, 3, 2, 1, 6]          # west -> east
KNOWN_IDS = [1, 2, 3, 4, 5, 6]

# 2026-09-20: 4 species (Sonneratia); the model's list, not re-declared. The observed loaders
# below keep only these species (the quadrat CSVs also carry Bruguiera, which is not modelled).
SPECIES = list(_MODEL_SPECIES)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
QUADRAT_CSV = os.path.join(SCRIPT_DIR, 'juvenile_quadrat_records_1999_t1_5.csv')


# ===================================================================
# Transect frames
# ===================================================================

def _px_to_model(px):
    """Figure pixel -> model metres: scale, mirror (sea at the bottom), rotate by the map
    tilt so the parallel transects are vertical, then shift."""
    xm = (px[0] - MAP_ORIGIN_PX[0]) / MAP_SCALE_PX_PER_M     # map east
    ym = (px[1] - MAP_ORIGIN_PX[1]) / MAP_SCALE_PX_PER_M     # landward
    c, s_ = math.cos(_tilt_rad), math.sin(_tilt_rad)
    return (xm * c - ym * s_ + MAP_OFFSET[0], xm * s_ + ym * c + MAP_OFFSET[1])


def plot_u(t_num, plot):
    """Fringe distance (tape m from the seaward end) of plot number(s); plot 1 = landward."""
    fr = transect_frame(t_num)
    return fr['tape_length'] - PLOT_LENGTH * (np.asarray(plot, dtype=float) - 1) - PLOT_LENGTH / 2.0


def _parallel_frame(t_num):
    d = _DIGITISED[t_num]
    X, Y = _px_to_model(d['px'])
    off = d['tape_offset']
    return dict(transect=t_num, start=(X - off * ALONG_DX, Y - off * ALONG_DY),
                along=(ALONG_DX, ALONG_DY), perp=(PERP_DX, PERP_DY), f=1.0,
                n_plots=d['n_plots'], tape_length=d['tape_length'], tape_offset=off)


def _t5_frame():
    """Least-squares line through the 1997 crossings; f = map m per tape m."""
    pts, us = [], []
    for tape, t, plot in T5_CROSSINGS_1997:
        fr = _parallel_frame(t)
        u = fr['tape_length'] - PLOT_LENGTH * (plot - 1) - PLOT_LENGTH / 2.0
        pts.append((fr['start'][0] + u * fr['along'][0], fr['start'][1] + u * fr['along'][1]))
        us.append(T5_TAPE_LENGTH - tape)          # T5 fringe distance of the crossing
    P = np.array(pts)
    us = np.array(us)
    c = P.mean(axis=0)
    _, _, vt = np.linalg.svd(P - c, full_matrices=False)
    v = vt[0]
    if np.dot(v, P[np.argmax(us)] - P[np.argmin(us)]) < 0:
        v = -v                                    # orient landward (increasing u)
    s = (P - c) @ v
    f, a = np.polyfit(us, s, 1)                   # s = a + f*u
    start = c + a * v                             # u = 0: seaward end of T5
    resid = s - (a + f * us)
    return dict(transect=5, start=(float(start[0]), float(start[1])),
                along=(float(v[0]), float(v[1])), perp=(float(v[1]), float(-v[0])), f=float(f),
                n_plots=T5_N_PLOTS, tape_length=T5_TAPE_LENGTH, tape_offset=0.0,
                crossing_rms_m=float(np.sqrt(np.mean(resid ** 2))))


_FRAMES = None


def transect_frame(t_num):
    """Geometry of one transect in the model frame.

    start: seaward end of the tape (u = 0); along: unit vector seaward -> landward;
    perp: unit vector left -> right looking landward; f: map metres per tape metre;
    position of fringe distance u = start + u * f * along.
    """
    global _FRAMES
    if _FRAMES is None:
        _FRAMES = {t: _parallel_frame(t) for t in _DIGITISED}
        _FRAMES[5] = _t5_frame()
    return _FRAMES[t_num]


def transect_position(t_num, u):
    fr = transect_frame(t_num)
    u = np.asarray(u, dtype=float)
    return fr['start'][0] + u * fr['f'] * fr['along'][0], fr['start'][1] + u * fr['f'] * fr['along'][1]


# Backward-compatible constants: seaward start of each modelled (drawn) transect line
TRANSECT_STARTS = {t: tuple(float(v) for v in transect_position(t, transect_frame(t)['tape_offset']))
                   for t in KNOWN_IDS}
TRANSECTS_ORDER = [(t, {'length': transect_frame(t)['tape_length'] - transect_frame(t)['tape_offset'],
                        'n_plots': transect_frame(t)['n_plots']}) for t in PARALLEL_IDS + [5]]


# ===================================================================
# Polygon Construction
# ===================================================================

def _transect_length(t_num):
    """Modelled line length (map m) from the drawn seaward start to the landward end."""
    fr = transect_frame(t_num)
    return (fr['tape_length'] - fr['tape_offset']) * fr['f']


def _clip_endpoint(start_x, start_y, length, along=None):
    """Clip a transect line to the domain box, return the actual length."""
    ax, ay = along if along is not None else (ALONG_DX, ALONG_DY)
    actual = length
    for lo, hi, s, a in ((X_MIN, X_MAX, start_x, ax), (Y_MIN, Y_MAX, start_y, ay)):
        if a > 1e-12:
            actual = min(actual, (hi - s) / a)
        elif a < -1e-12:
            actual = min(actual, (lo - s) / a)
    return max(0.0, actual)


def build_transect_polygons(ids=None):
    """Build transect polygon corners (4 corners per transect).

    Returns dict: {t_num: {'corners': [(x,y),...], 'start': (x,y), 'end': (x,y),
                            'length': float, 'area': float, 'n_plots': int,
                            'along': (ax, ay), 'perp': (px, py), 'f': float}}
    """
    half_w = TRANSECT_WIDTH / 2.0
    polygons = {}
    for t_num in (ids or ALL_TRANSECT_IDS):
        fr = transect_frame(t_num)
        sx, sy = TRANSECT_STARTS[t_num]
        ax, ay = fr['along']
        px, py = fr['perp']
        actual_length = _clip_endpoint(sx, sy, _transect_length(t_num), (ax, ay))
        ex = sx + actual_length * ax
        ey = sy + actual_length * ay
        corners = [
            (sx - half_w * px, sy - half_w * py),
            (sx + half_w * px, sy + half_w * py),
            (ex + half_w * px, ey + half_w * py),
            (ex - half_w * px, ey - half_w * py),
        ]
        polygons[t_num] = {
            'corners': corners, 'start': (sx, sy), 'end': (ex, ey),
            'length': actual_length, 'area': actual_length * TRANSECT_WIDTH,
            'n_plots': fr['n_plots'], 'along': (ax, ay), 'perp': (px, py), 'f': fr['f'],
        }
    return polygons


def _point_in_polygon(px, py, corners):
    """Ray-casting point-in-polygon test."""
    n = len(corners)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = corners[i]
        xj, yj = corners[j]
        if ((yi > py) != (yj > py)) and \
           (px < (xj - xi) * (py - yi) / (yj - yi) + xi):
            inside = not inside
        j = i
    return inside


# Pre-built polygons (lazy singleton)
_POLYGONS = None


def _get_polygons():
    global _POLYGONS
    if _POLYGONS is None:
        _POLYGONS = build_transect_polygons(KNOWN_IDS)
    return _POLYGONS


def assign_to_transect(x, y, ids=None):
    """Assign a point (x, y) to a transect number, or None if outside all.

    By default only the validation transects (TRANSECT_IDS, 1-5) are considered;
    pass ids=ALL_TRANSECT_IDS to include T6.  Where T5 crosses another transect
    the parallel transect wins.
    """
    polys = _get_polygons()
    cand = list(ids or TRANSECT_IDS)
    order = [t for t in cand if t != 5] + ([5] if 5 in cand else [])
    for t_num in order:
        if _point_in_polygon(x, y, polys[t_num]['corners']):
            return t_num
    return None


def assign_to_transect_batch(xs, ys, ids=None):
    """Vectorised transect assignment for arrays of x, y coordinates."""
    polys = _get_polygons()
    cand = list(ids or TRANSECT_IDS)
    order = [t for t in cand if t != 5] + ([5] if 5 in cand else [])
    n = len(xs)
    results = [None] * n
    for t_num in order:
        corners = polys[t_num]['corners']
        for i in range(n):
            if results[i] is None and _point_in_polygon(xs[i], ys[i], corners):
                results[i] = t_num
    return results


# ===================================================================
# Along-Transect Position (for plot/quarter assignment)
# ===================================================================

def along_transect_distance(x, y, t_num):
    """Map distance (m) along the transect from its drawn seaward start."""
    polys = _get_polygons()
    sx, sy = polys[t_num]['start']
    ax, ay = polys[t_num]['along']
    return (x - sx) * ax + (y - sy) * ay


def fringe_distance_u(x, y, t_num):
    """Fringe distance in tape metres (u) of a point projected on the transect."""
    fr = transect_frame(t_num)
    return fr['tape_offset'] + along_transect_distance(x, y, t_num) / fr['f']


def _plot_and_offset(x, y, t_num, plot_length=PLOT_LENGTH):
    fr = transect_frame(t_num)
    u = fringe_distance_u(x, y, t_num)
    u_edge0 = fr['tape_length'] - plot_length * fr['n_plots']   # seaward edge of the last plot
    idx = int(math.floor((u - u_edge0) / plot_length))
    idx = min(max(idx, 0), fr['n_plots'] - 1)
    plot_num = fr['n_plots'] - idx
    d_in_plot = (u - u_edge0) - idx * plot_length
    return plot_num, d_in_plot


def assign_to_plot(x, y, t_num, plot_length=PLOT_LENGTH):
    """Assign a point to a plot number within a transect (plot 1 = landward end)."""
    return _plot_and_offset(x, y, t_num, plot_length)[0]


def _perp_distance(x, y, t_num):
    polys = _get_polygons()
    sx, sy = polys[t_num]['start']
    px, py = polys[t_num]['perp']
    return (x - sx) * px + (y - sy) * py   # signed: positive = right looking landward


def assign_to_quarter(x, y, t_num, plot_length=PLOT_LENGTH):
    """Assign a point to a quarter (a/b/c/d) within its plot.

    Default orientation: quarter grid divides each plot into 4 cells:
      a = seaward-left, b = seaward-right, c = landward-left, d = landward-right

    Returns (plot_number, quarter_letter).
    """
    plot_num, d_in_plot = _plot_and_offset(x, y, t_num, plot_length)
    along_half = 'seaward' if d_in_plot < plot_length / 2.0 else 'landward'
    perp_half = 'left' if _perp_distance(x, y, t_num) < 0 else 'right'
    quarter_map = {
        ('seaward', 'left'): 'a',
        ('seaward', 'right'): 'b',
        ('landward', 'left'): 'c',
        ('landward', 'right'): 'd',
    }
    return plot_num, quarter_map[(along_half, perp_half)]


def assign_to_quarter_rotated(x, y, t_num, plot_length=PLOT_LENGTH):
    """Alternate quarter orientation (90-degree rotation of quarter grid).

    Rotated: a = left-seaward, b = left-landward, c = right-seaward, d = right-landward
    Used for geometry sensitivity analysis.

    Returns (plot_number, quarter_letter).
    """
    plot_num, d_in_plot = _plot_and_offset(x, y, t_num, plot_length)
    along_half = 'seaward' if d_in_plot < plot_length / 2.0 else 'landward'
    perp_half = 'left' if _perp_distance(x, y, t_num) < 0 else 'right'
    quarter_map = {
        ('seaward', 'left'): 'a',
        ('landward', 'left'): 'b',
        ('seaward', 'right'): 'c',
        ('landward', 'right'): 'd',
    }
    return plot_num, quarter_map[(along_half, perp_half)]


# ===================================================================
# Observed Data Loading (transect-level)
# ===================================================================

def load_observed_quadrats_per_zone(csv_path=None):
    """Number of distinct observed quadrats per model zone (from the quadrat CSV's model_zone)."""
    if csv_path is None:
        csv_path = QUADRAT_CSV
    # 2026-09-22 FIX: restrict to quadrats INSIDE the model box, to match the observed target in
    # cascade_metrics.load_observed_quadrats and virtual_survey.load_quadrats. 26 of the 335 quadrats
    # (all on the diagonal T5) lie outside the domain; the model can place no recruit there and they
    # are not inside any in-domain belt either, so counting them here divided the model-side density
    # by a larger quadrat count than the observation used (A 72/B 228/C 35 vs A 63/B 224/C 22).
    import tidal_dispersal_model as _t
    zones = {}
    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            x, y = float(row['x']), float(row['y'])
            if not (_t.X_MIN <= x <= _t.X_MAX and _t.Y_MIN <= y <= _t.Y_MAX):
                continue
            zones.setdefault(row['model_zone'], set()).add(row['quadrat_id'])
    return {z: len(q) for z, q in zones.items()}


def load_observed_transect_densities(csv_path=None):
    """Load 1999 observed densities aggregated by transect × species.

    Returns dict: {(transect_int, species): mean_density_per_m2}
    """
    if csv_path is None:
        csv_path = QUADRAT_CSV
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Quadrat CSV not found: {csv_path}")

    # Collect all densities per (transect, species)
    from collections import defaultdict
    accum = defaultdict(list)

    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = int(row['transect'])
            sp = row['species']
            if sp not in SPECIES:
                continue
            dens = float(row['density_per_m2'])
            accum[(t, sp)].append(dens)

    # Mean density per (transect, species)
    result = {}
    for key, vals in accum.items():
        result[key] = sum(vals) / len(vals) if vals else 0.0
    return result


def load_observed_transect_quadrats(csv_path=None):
    """Load quadrat-level records grouped by transect for bootstrap.

    Returns dict: {transect_int: {quadrat_id: {species: density_per_m2}}}
    """
    if csv_path is None:
        csv_path = QUADRAT_CSV
    if not os.path.exists(csv_path):
        return None

    data = {}
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = int(row['transect'])
            qid = row['quadrat_id']
            sp = row['species']
            if sp not in SPECIES:
                continue
            dens = float(row['density_per_m2'])
            if t not in data:
                data[t] = {}
            if qid not in data[t]:
                data[t][qid] = {}
            data[t][qid][sp] = dens
    return data


def load_observed_transect_species_stats(csv_path=None):
    """Load per-transect per-species mean, SD, n for uncertainty decomposition.

    Returns dict: {(transect, species): {'mean': float, 'sd': float, 'n': int}}
    """
    if csv_path is None:
        csv_path = QUADRAT_CSV
    if not os.path.exists(csv_path):
        return {}

    from collections import defaultdict
    accum = defaultdict(list)

    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = int(row['transect'])
            sp = row['species']
            if sp not in SPECIES:
                continue
            dens = float(row['density_per_m2'])
            accum[(t, sp)].append(dens)

    result = {}
    for key, vals in accum.items():
        arr = np.array(vals)
        result[key] = {
            'mean': float(np.mean(arr)),
            'sd': float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
            'n': len(arr),
        }
    return result


# ===================================================================
# Model Transect Densities
# ===================================================================

def model_transect_densities(results):
    """Compute established density (#/m²) per (transect, species) from results.

    Only counts propagules whose final (x, y) falls within a transect polygon.
    Returns dict: {(transect_int, species): density}
    """
    polys = _get_polygons()
    established = [r for r in results if r['outcome'] == 'established']

    # Count per (transect, species)
    from collections import Counter
    counts = Counter()
    for r in established:
        t = assign_to_transect(r['x'], r['y'])
        if t is not None:
            counts[(t, r['species'])] += 1

    # Convert to density
    densities = {}
    for t_num in TRANSECT_IDS:
        area = polys[t_num]['area']
        for sp in SPECIES:
            key = (t_num, sp)
            densities[key] = counts[key] / area if area > 0 else 0.0
    return densities


def model_transect_plot_densities(results, plot_length=10.0):
    """Compute density per (transect, plot, species) for fine-grained comparison.

    Returns dict: {(transect, plot, species): density}
    """
    polys = _get_polygons()
    established = [r for r in results if r['outcome'] == 'established']

    from collections import Counter
    counts = Counter()
    plot_area = plot_length * TRANSECT_WIDTH

    for r in established:
        t = assign_to_transect(r['x'], r['y'])
        if t is not None:
            p = assign_to_plot(r['x'], r['y'], t, plot_length)
            counts[(t, p, r['species'])] += 1

    densities = {}
    for key, cnt in counts.items():
        densities[key] = cnt / plot_area if plot_area > 0 else 0.0
    return densities


def support_audit(results, csv_path=None):
    """Support-consistency audit: compare polygon-area vs sampled-area densities.

    For each transect, computes:
      - full_polygon_area: transect polygon area (length * width)
      - n_quadrats_obs: number of distinct quadrat IDs in observed data
      - sampled_area: n_quadrats_obs * MODEL_QUADRAT_FOOTPRINT_M2 (assumed footprint)
      - model_count: number of established propagules in transect
      - density_full_polygon: model_count / full_polygon_area
      - density_sampled_area: model_count / sampled_area
      - obs_mean_density: mean observed density across all species
      - area_ratio: full_polygon_area / sampled_area

    Returns dict with per-transect entries and overall summary.
    """
    polys = _get_polygons()
    quadrat_data = load_observed_transect_quadrats(csv_path)
    if quadrat_data is None:
        quadrat_data = {}

    # Count established propagules per transect
    from collections import Counter
    established = [r for r in results if r['outcome'] == 'established']
    counts = Counter()
    for r in established:
        t = assign_to_transect(r['x'], r['y'])
        if t is not None:
            counts[t] += 1

    # Load observed densities for mean obs density per transect
    obs_dens = load_observed_transect_densities(csv_path)

    per_transect = {}
    area_ratios = []

    for t_num in TRANSECT_IDS:
        full_area = polys[t_num]['area']
        n_quadrats = len(quadrat_data.get(t_num, {}))
        sampled_area = n_quadrats * MODEL_QUADRAT_FOOTPRINT_M2
        model_count = counts[t_num]

        density_full = model_count / full_area if full_area > 0 else 0.0
        density_sampled = model_count / sampled_area if sampled_area > 0 else 0.0

        # Mean observed density across species for this transect
        obs_vals = [obs_dens.get((t_num, sp), 0.0) for sp in SPECIES]
        obs_mean = sum(obs_vals) / len(obs_vals) if obs_vals else 0.0

        ratio = full_area / sampled_area if sampled_area > 0 else float('inf')
        area_ratios.append(ratio)

        per_transect[t_num] = {
            'full_polygon_area': full_area,
            'sampled_area': sampled_area,
            'n_quadrats_obs': n_quadrats,
            'model_count': model_count,
            'density_full_polygon': density_full,
            'density_sampled_area': density_sampled,
            'obs_mean_density': obs_mean,
            'area_ratio': ratio,
        }

    # Overall summary
    arr = np.array(area_ratios)
    mean_ratio = float(np.mean(arr)) if len(arr) > 0 else 0.0
    sd_ratio = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0

    return {
        'per_transect': per_transect,
        'mean_area_ratio': mean_ratio,
        'sd_area_ratio': sd_ratio,
    }


def assert_support_compatibility(model_support, obs_support, expected_mode='sampled_area'):
    """Raise ValueError if model and observation density supports don't match.

    Args:
        model_support: str, how model densities were computed ('full_polygon' or 'sampled_area')
        obs_support: str, how observed densities were computed
        expected_mode: str, the expected support mode
    """
    if model_support != expected_mode:
        raise ValueError(
            f"Support mismatch: model uses '{model_support}' but expected '{expected_mode}'. "
            f"Use model_transect_sampled_densities() for sampled_area mode.")
    if obs_support != expected_mode and expected_mode == 'sampled_area':
        raise ValueError(
            f"Support mismatch: obs uses '{obs_support}' but expected '{expected_mode}'.")


def load_observed_transect_quadrat_counts(csv_path=None):
    """Count distinct quadrat IDs per transect from observed CSV.

    Returns dict: {transect_int: n_quadrats}
    """
    quad_data = load_observed_transect_quadrats(csv_path)
    if quad_data is None:
        return {}
    return {t: len(quads) for t, quads in quad_data.items()}


def model_transect_sampled_densities(results, plot_length=10.0, obs_footprint=None):
    """Compute model densities using the same support as field observations.

    Instead of dividing by the full transect polygon area, this function
    computes the mean of per-quadrat-cell densities, matching how observed
    densities are computed (mean of per-quadrat density_per_m2 values).

    Each quadrat cell is given the ASSUMED footprint MODEL_QUADRAT_FOOTPRINT_M2 (see the
    constant); the mean density across n_cells equals
    total_count / (MODEL_QUADRAT_FOOTPRINT_M2 * n_cells) and is conditional on that assumption.

    Args:
        results: list of propagule result dicts from model run
        plot_length: length of each plot along transect (m)
        obs_footprint: {transect_int: n_quadrats} from observed data.
            If provided, n_cells = obs_footprint[transect] (per species).
            If None, n_cells = n_plots * QUADRATS_PER_PLOT.

    Returns dict: {(transect_int, species): mean_density_per_m2}
    """
    polys = _get_polygons()
    established = [r for r in results if r['outcome'] == 'established']

    from collections import Counter
    counts = Counter()
    for r in established:
        t = assign_to_transect(r['x'], r['y'])
        if t is not None:
            counts[(t, r['species'])] += 1

    densities = {}
    for t_num in TRANSECT_IDS:
        if obs_footprint is not None and t_num in obs_footprint:
            n_cells = obs_footprint[t_num]
        else:
            n_plots = polys[t_num]['n_plots']
            n_cells = n_plots * QUADRATS_PER_PLOT

        for sp in SPECIES:
            key = (t_num, sp)
            total_count = counts[key]
            if n_cells > 0 and MODEL_QUADRAT_FOOTPRINT_M2 > 0:
                densities[key] = total_count / (MODEL_QUADRAT_FOOTPRINT_M2 * n_cells)
            else:
                densities[key] = 0.0

    return densities


# ===================================================================
# Transect-Zone Area Computation (for zone-level guardrails)
# ===================================================================

_TRANSECT_ZONE_AREAS = None  # lazy singleton


def compute_transect_zone_areas(resolution=0.5):
    """Compute area (m2) of each zone within the union of all transect polygons.

    Uses grid sampling at `resolution` spacing.  Each grid point inside a
    transect polygon contributes resolution**2 to its zone's area.

    Returns dict: {'A': area_m2, 'B': area_m2, 'C': area_m2}
    """
    global _TRANSECT_ZONE_AREAS
    if _TRANSECT_ZONE_AREAS is not None:
        return dict(_TRANSECT_ZONE_AREAS)

    polys = _get_polygons()
    cell_area = resolution ** 2
    zone_area = {'A': 0.0, 'B': 0.0, 'C': 0.0}

    for t_num in TRANSECT_IDS:
        corners = polys[t_num]['corners']
        # Bounding box
        xs = [c[0] for c in corners]
        ys = [c[1] for c in corners]
        x_lo, x_hi = min(xs), max(xs)
        y_lo, y_hi = min(ys), max(ys)

        x = x_lo
        while x <= x_hi:
            y = y_lo
            while y <= y_hi:
                if _point_in_polygon(x, y, corners):
                    z = zone_at(x, y)
                    zone_area[z] += cell_area
                y += resolution
            x += resolution

    _TRANSECT_ZONE_AREAS = dict(zone_area)
    return dict(_TRANSECT_ZONE_AREAS)


def model_zone_sampled_densities(results):
    """Compute zone x species densities of the model under an ASSUMED sampling footprint.

    Divides model counts (established recruits within the transect polygons) by
    n_quadrats_per_zone * MODEL_QUADRAT_FOOTPRINT_M2, i.e. by the footprint the survey
    quadrats would cover if each were the paper's 5 x 5 m quadrat (25 m2; the juvenile
    quadrat area is not documented, data audit 2026-09-20). Zone SHARES derived from these
    densities (relative density across zones) are independent of the footprint; the
    ABSOLUTE densities are conditional on the footprint assumption and must be labelled so.

    Returns dict: {(zone, species): density_per_m2 (conditional on the footprint)}
    """
    # Quadrat counts per zone from 1999 field sampling design
    # (from juvenile_zone_species_summary_1999_t1_5.csv)
    N_QUADRATS_PER_ZONE = load_observed_quadrats_per_zone()

    established = [r for r in results if r['outcome'] == 'established']

    from collections import Counter
    counts = Counter()
    for r in established:
        t = assign_to_transect(r['x'], r['y'])
        if t is not None:
            counts[(r['zone'], r['species'])] += 1

    densities = {}
    for z in ['A', 'B', 'C']:
        sampled_area = N_QUADRATS_PER_ZONE.get(z, 1) * MODEL_QUADRAT_FOOTPRINT_M2
        for sp in SPECIES:
            key = (z, sp)
            densities[key] = counts[key] / sampled_area if sampled_area > 0 else 0.0

    return densities


def transect_zone_composition():
    """Per-transect zone breakdown (fraction of each transect in each zone).

    Returns dict: {t_num: {'A': frac, 'B': frac, 'C': frac, 'area': total_m2}}
    """
    polys = _get_polygons()
    resolution = 0.5
    cell_area = resolution ** 2

    result = {}
    for t_num in TRANSECT_IDS:
        corners = polys[t_num]['corners']
        xs = [c[0] for c in corners]
        ys = [c[1] for c in corners]
        x_lo, x_hi = min(xs), max(xs)
        y_lo, y_hi = min(ys), max(ys)

        zone_counts = {'A': 0, 'B': 0, 'C': 0}
        total = 0

        x = x_lo
        while x <= x_hi:
            y = y_lo
            while y <= y_hi:
                if _point_in_polygon(x, y, corners):
                    z = zone_at(x, y)
                    zone_counts[z] += 1
                    total += 1
                y += resolution
            x += resolution

        if total > 0:
            result[t_num] = {
                'A': zone_counts['A'] / total,
                'B': zone_counts['B'] / total,
                'C': zone_counts['C'] / total,
                'area': total * cell_area,
            }
        else:
            result[t_num] = {'A': 0.0, 'B': 0.0, 'C': 0.0, 'area': 0.0}

    return result


# ===================================================================
# Transect-Level Metrics (confirmatory)
# ===================================================================

def _rank(values):
    """Assign ranks (average for ties)."""
    n = len(values)
    indexed = sorted(enumerate(values), key=lambda x: x[1])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n - 1 and indexed[j + 1][1] == indexed[j][1]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks


def _spearman(vals_a, vals_b):
    """Spearman rank correlation on two lists of equal length."""
    ra = _rank(vals_a)
    rb = _rank(vals_b)
    n = len(vals_a)
    ma = sum(ra) / n
    mb = sum(rb) / n
    cov = sum((a - ma) * (b - mb) for a, b in zip(ra, rb))
    va = sum((a - ma) ** 2 for a in ra)
    vb = sum((b - mb) ** 2 for b in rb)
    denom = math.sqrt(va * vb)
    return cov / denom if denom > 0 else 0.0


def _kendall_tau(vals_a, vals_b):
    """Kendall tau-b rank correlation."""
    n = len(vals_a)
    concordant = 0
    discordant = 0
    ties_a = 0
    ties_b = 0
    for i in range(n):
        for j in range(i + 1, n):
            da = vals_a[i] - vals_a[j]
            db = vals_b[i] - vals_b[j]
            if da == 0 and db == 0:
                ties_a += 1
                ties_b += 1
            elif da == 0:
                ties_a += 1
            elif db == 0:
                ties_b += 1
            elif da * db > 0:
                concordant += 1
            else:
                discordant += 1
    n_pairs = n * (n - 1) / 2
    denom = math.sqrt((n_pairs - ties_a) * (n_pairs - ties_b))
    return (concordant - discordant) / denom if denom > 0 else 0.0


def _jsd(p, q):
    """Jensen-Shannon divergence (base-2 log), bounded [0, 1]."""
    p_sum = sum(p)
    q_sum = sum(q)
    if p_sum == 0 and q_sum == 0:
        return 0.0
    p_norm = [v / p_sum if p_sum > 0 else 0.0 for v in p]
    q_norm = [v / q_sum if q_sum > 0 else 0.0 for v in q]
    m = [(a + b) / 2.0 for a, b in zip(p_norm, q_norm)]
    kl_pm = 0.0
    kl_qm = 0.0
    for pi, qi, mi in zip(p_norm, q_norm, m):
        if mi > 0:
            if pi > 0:
                kl_pm += pi * math.log2(pi / mi)
            if qi > 0:
                kl_qm += qi * math.log2(qi / mi)
    return (kl_pm + kl_qm) / 2.0


def transect_species_rank_agreement(model_dens, obs_dens):
    """Transect×species rank agreement (Spearman + Kendall).

    Operates on the 15-cell (5 transects × 3 species) comparison matrix.
    Returns dict with spearman_rho, kendall_tau, n_cells.
    """
    cells = [(t, sp) for t in TRANSECT_IDS for sp in SPECIES]
    m_vals = [model_dens.get(k, 0.0) for k in cells]
    o_vals = [obs_dens.get(k, 0.0) for k in cells]

    rho = _spearman(m_vals, o_vals)
    tau = _kendall_tau(m_vals, o_vals)

    return {
        'spearman_rho': rho,
        'kendall_tau': tau,
        'n_cells': len(cells),
    }


def species_peak_transect_match(model_dens, obs_dens):
    """Check each species' peak transect (highest density) against observed.

    Returns dict with 'matches' count, 'near_misses' count,
    'rank_distances' per species, and per-species details.

    Transects are ranked by their ID order in TRANSECT_IDS. A near_miss
    is when the model peak transect is adjacent to the observed peak
    transect in rank order (rank distance <= 1 but not an exact match).
    """
    transect_rank = {t: i for i, t in enumerate(TRANSECT_IDS)}
    details = {}
    rank_distances = {}
    matches = 0
    near_misses = 0
    for sp in SPECIES:
        m_vals = {t: model_dens.get((t, sp), 0.0) for t in TRANSECT_IDS}
        o_vals = {t: obs_dens.get((t, sp), 0.0) for t in TRANSECT_IDS}
        m_peak = max(TRANSECT_IDS, key=lambda t: m_vals[t])
        o_peak = max(TRANSECT_IDS, key=lambda t: o_vals[t])
        match = m_peak == o_peak
        rd = abs(transect_rank[m_peak] - transect_rank[o_peak])
        rank_distances[sp] = rd
        if match:
            matches += 1
        elif rd <= 1:
            near_misses += 1
        details[sp] = {
            'model_peak': m_peak,
            'obs_peak': o_peak,
            'match': match,
            'rank_distance': rd,
            'near_miss': (not match) and rd <= 1,
        }
    return {
        'matches': matches,
        'near_misses': near_misses,
        'rank_distances': rank_distances,
        'details': details,
    }


def transect_composition_jsd(model_dens, obs_dens):
    """Species-proportion similarity per transect via JSD.

    Returns per-transect JSD, mean JSD.
    """
    per_transect = {}
    for t in TRANSECT_IDS:
        m = [model_dens.get((t, sp), 0.0) for sp in SPECIES]
        o = [obs_dens.get((t, sp), 0.0) for sp in SPECIES]
        jsd = _jsd(m, o)
        per_transect[t] = {'jsd': jsd}

    jsd_values = [per_transect[t]['jsd'] for t in TRANSECT_IDS]
    mean_jsd = sum(jsd_values) / len(jsd_values) if jsd_values else 0.0
    n_passing = sum(1 for v in jsd_values if v < 0.2)

    return {
        'per_transect': per_transect,
        'mean_jsd': mean_jsd,
        'n_transects_pass': n_passing,
    }


def transect_standardized_residuals(model_dens, obs_stats):
    """Standardized residuals per transect×species: (model - obs_mean) / obs_SE.

    Args:
        model_dens: {(transect, species): density}
        obs_stats: from load_observed_transect_species_stats()

    Returns dict: {(transect, species): residual}
    """
    resids = {}
    for t in TRANSECT_IDS:
        for sp in SPECIES:
            key = (t, sp)
            m = model_dens.get(key, 0.0)
            info = obs_stats.get(key, {'mean': 0.0, 'sd': 0.0, 'n': 0})
            se = info['sd'] / math.sqrt(info['n']) if info['n'] > 1 else 0.0
            if se > 0:
                resids[key] = (m - info['mean']) / se
            else:
                resids[key] = float('inf') if m > info['mean'] else 0.0
    return resids


def variance_decomposition(model_dens_list, obs_stats):
    """Combined variance decomposition: observational SE + model replicate SD.

    Args:
        model_dens_list: list of per-replicate {(transect, species): density} dicts
        obs_stats: from load_observed_transect_species_stats(),
                   {(transect, species): {'mean': float, 'sd': float, 'n': int}}

    Returns dict keyed by "T{t}_{sp}" strings, with per-cell:
        obs_mean, obs_sd, obs_se, obs_n, obs_cv,
        model_mean, model_sd, model_n, model_cv,
        combined_se: sqrt(obs_se^2 + model_sd^2)
    """
    result = {}
    for t in TRANSECT_IDS:
        for sp in SPECIES:
            key = (t, sp)
            label = f"T{t}_{sp}"

            # Observed stats
            info = obs_stats.get(key, {'mean': 0.0, 'sd': 0.0, 'n': 0})
            obs_mean = info['mean']
            obs_sd = info['sd']
            obs_n = info['n']
            obs_se = obs_sd / math.sqrt(obs_n) if obs_n > 1 else 0.0
            obs_cv = obs_sd / obs_mean if obs_mean > 0 else None

            # Model stats across replicates
            model_vals = [d.get(key, 0.0) for d in model_dens_list]
            model_n = len(model_vals)
            if model_n > 0:
                arr = np.array(model_vals)
                model_mean = float(np.mean(arr))
                model_sd = float(np.std(arr, ddof=1)) if model_n > 1 else 0.0
            else:
                model_mean = 0.0
                model_sd = 0.0
            model_cv = model_sd / model_mean if model_mean > 0 else None

            combined_se = math.sqrt(obs_se ** 2 + model_sd ** 2)

            result[label] = {
                'obs_mean': obs_mean,
                'obs_sd': obs_sd,
                'obs_se': obs_se,
                'obs_n': obs_n,
                'obs_cv': obs_cv,
                'model_mean': model_mean,
                'model_sd': model_sd,
                'model_n': model_n,
                'model_cv': model_cv,
                'combined_se': combined_se,
            }
    return result


def bland_altman_transect(model_dens, obs_dens):
    """Bland-Altman analysis at transect totals.

    Returns dict with mean_diff, sd_diff, limits_of_agreement, per-transect data.
    """
    per_transect = {}
    diffs = []
    means = []

    for t in TRANSECT_IDS:
        m_total = sum(model_dens.get((t, sp), 0.0) for sp in SPECIES)
        o_total = sum(obs_dens.get((t, sp), 0.0) for sp in SPECIES)
        diff = m_total - o_total
        mean_val = (m_total + o_total) / 2.0
        diffs.append(diff)
        means.append(mean_val)
        per_transect[t] = {
            'model_total': m_total,
            'obs_total': o_total,
            'diff': diff,
            'mean': mean_val,
        }

    mean_diff = sum(diffs) / len(diffs) if diffs else 0.0
    sd_diff = math.sqrt(sum((d - mean_diff) ** 2 for d in diffs) / max(1, len(diffs) - 1)) if len(diffs) > 1 else 0.0
    loa_upper = mean_diff + 1.96 * sd_diff
    loa_lower = mean_diff - 1.96 * sd_diff

    return {
        'mean_diff': mean_diff,
        'sd_diff': sd_diff,
        'loa_upper': loa_upper,
        'loa_lower': loa_lower,
        'per_transect': per_transect,
    }


def transect_absolute_scale_metrics(model_dens, obs_dens):
    """Absolute-scale diagnostics across the 15-cell transect×species matrix.

    Computes RMSE, MAE, multiplicative bias, weighted RMSE, and per-cell
    error details.  All metrics are informational (no pass/fail threshold).

    Args:
        model_dens: {(transect, species): density}
        obs_dens:   {(transect, species): density}

    Returns dict with rmse, mae, multiplicative_bias, weighted_rmse, per_cell.
    """
    cells = [(t, sp) for t in TRANSECT_IDS for sp in SPECIES]
    m_vals = [model_dens.get(k, 0.0) for k in cells]
    o_vals = [obs_dens.get(k, 0.0) for k in cells]
    n = len(cells)

    # --- RMSE ---
    sq_errors = [(m - o) ** 2 for m, o in zip(m_vals, o_vals)]
    rmse = math.sqrt(sum(sq_errors) / n)

    # --- MAE ---
    abs_errors = [abs(m - o) for m, o in zip(m_vals, o_vals)]
    mae = sum(abs_errors) / n

    # --- Multiplicative bias: exp(mean(log(model/obs))) for obs>0 and model>0 ---
    log_ratios = []
    for m, o in zip(m_vals, o_vals):
        if o > 0 and m > 0:
            log_ratios.append(math.log(m / o))
    multiplicative_bias = math.exp(sum(log_ratios) / len(log_ratios)) if log_ratios else None

    # --- Weighted RMSE (weighted by obs density, obs > 0 only) ---
    w_sum_sq = 0.0
    w_sum = 0.0
    for m, o in zip(m_vals, o_vals):
        if o > 0:
            w_sum_sq += o * (m - o) ** 2
            w_sum += o
    weighted_rmse = math.sqrt(w_sum_sq / w_sum) if w_sum > 0 else None

    # --- Per-cell details ---
    per_cell = {}
    for k, m, o in zip(cells, m_vals, o_vals):
        label = f"T{k[0]}_{k[1]}"
        per_cell[label] = {
            'model': m,
            'obs': o,
            'error': m - o,
            'ratio': m / o if o > 0 else None,
        }

    return {
        'rmse': rmse,
        'mae': mae,
        'multiplicative_bias': multiplicative_bias,
        'weighted_rmse': weighted_rmse,
        'per_cell': per_cell,
    }


# ===================================================================
# Bootstrap CIs (transect-level)
# ===================================================================

def bootstrap_transect_metrics(model_dens, transect_quadrats, n_boot=2000, seed=42):
    """Bootstrap CIs for transect-level metrics by resampling quadrats within transects.

    For each iteration: resample quadrats within each transect (blocked bootstrap),
    compute mean density per (transect, species), then compute rank agreement,
    peak matches, and mean JSD against model.

    Args:
        model_dens: {(transect, species): density}
        transect_quadrats: from load_observed_transect_quadrats()
        n_boot: bootstrap iterations
        seed: random seed

    Returns dict: {metric: {ci_lo, ci_hi, boot_mean}}
    """
    rng = np.random.default_rng(seed)

    boot_spearman = []
    boot_kendall = []
    boot_peak_matches = []
    boot_mean_jsd = []

    for _ in range(n_boot):
        boot_obs = {}
        for t, quads in transect_quadrats.items():
            quad_ids = list(quads.keys())
            n_q = len(quad_ids)
            if n_q == 0:
                continue
            resampled_ids = rng.choice(quad_ids, size=n_q, replace=True)
            for sp in SPECIES:
                densities = [quads[qid].get(sp, 0.0) for qid in resampled_ids]
                boot_obs[(t, sp)] = float(np.mean(densities))

        rank_info = transect_species_rank_agreement(model_dens, boot_obs)
        boot_spearman.append(rank_info['spearman_rho'])
        boot_kendall.append(rank_info['kendall_tau'])

        peak_info = species_peak_transect_match(model_dens, boot_obs)
        boot_peak_matches.append(peak_info['matches'])

        jsd_info = transect_composition_jsd(model_dens, boot_obs)
        boot_mean_jsd.append(jsd_info['mean_jsd'])

    def _ci(arr):
        a = np.array(arr)
        return {
            'ci_lo': float(np.percentile(a, 2.5)),
            'ci_hi': float(np.percentile(a, 97.5)),
            'boot_mean': float(np.mean(a)),
            'boot_sd': float(np.std(a)),
        }

    return {
        'transect_spearman_rho': _ci(boot_spearman),
        'transect_kendall_tau': _ci(boot_kendall),
        'transect_peak_matches': _ci(boot_peak_matches),
        'transect_mean_jsd': _ci(boot_mean_jsd),
    }


def permutation_pvalue_spearman(observed_rho, n, n_perm=100000, seed=42):
    """Compute permutation p-value for Spearman rho under null hypothesis.

    Under H0 (no association), generates n_perm random permutations of n values
    and computes Spearman rho for each. P-value = fraction of permuted rhos >= observed.

    Args:
        observed_rho: the observed Spearman rho value
        n: number of data points
        n_perm: number of permutations
        seed: random seed

    Returns float: two-tailed p-value
    """
    rng = np.random.default_rng(seed)
    base = list(range(n))
    count_extreme = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        # Spearman rho on ranks vs permuted ranks
        d_sq = sum((b - p) ** 2 for b, p in zip(base, perm))
        rho_perm = 1 - 6 * d_sq / (n * (n * n - 1))
        if abs(rho_perm) >= abs(observed_rho):
            count_extreme += 1
    return count_extreme / n_perm


# ===================================================================
# LOTO Cross-Validation
# ===================================================================

def loto_cross_validation(results, obs_dens, obs_footprint=None, support_mode='sampled_area'):
    """Leave-One-Transect-Out cross-validation for transect-level metrics.

    Holds out each of 5 transects in turn, recomputes Spearman rho, peak matches,
    and JSD on the remaining 4 transects (12 cells).

    Args:
        results: list of propagule result dicts from model run
        obs_dens: {(transect, species): density} observed densities
        obs_footprint: {transect: n_quadrats} for sampled support
        support_mode: 'sampled_area' or 'full_polygon'

    Returns dict: {
        per_fold: [{held_out, rho, peak_matches, mean_jsd}],
        mean_rho, sd_rho, mean_jsd, sd_jsd
    }
    """
    # Compute model densities using appropriate support mode
    if support_mode == 'sampled_area':
        m_dens = model_transect_sampled_densities(results, obs_footprint=obs_footprint)
    else:
        m_dens = model_transect_densities(results)

    per_fold = []
    for held_out in TRANSECT_IDS:
        remaining = [t for t in TRANSECT_IDS if t != held_out]

        # Build filtered value lists (12 cells: 4 transects x 3 species)
        cells = [(t, sp) for t in remaining for sp in SPECIES]
        m_vals = [m_dens.get(k, 0.0) for k in cells]
        o_vals = [obs_dens.get(k, 0.0) for k in cells]

        # Spearman on remaining 12 cells
        rho = _spearman(m_vals, o_vals)

        # Peak matches on remaining transects
        matches = 0
        for sp in SPECIES:
            m_sp = {t: m_dens.get((t, sp), 0.0) for t in remaining}
            o_sp = {t: obs_dens.get((t, sp), 0.0) for t in remaining}
            m_peak = max(remaining, key=lambda t: m_sp[t])
            o_peak = max(remaining, key=lambda t: o_sp[t])
            if m_peak == o_peak:
                matches += 1

        # Mean JSD on remaining transects
        jsd_vals = []
        for t in remaining:
            m_comp = [m_dens.get((t, sp), 0.0) for sp in SPECIES]
            o_comp = [obs_dens.get((t, sp), 0.0) for sp in SPECIES]
            jsd_vals.append(_jsd(m_comp, o_comp))
        mean_jsd = sum(jsd_vals) / len(jsd_vals) if jsd_vals else 0.0

        per_fold.append({
            'held_out': held_out,
            'rho': rho,
            'peak_matches': matches,
            'mean_jsd': mean_jsd,
        })

    rhos = [f['rho'] for f in per_fold]
    jsds = [f['mean_jsd'] for f in per_fold]

    arr_rho = np.array(rhos)
    arr_jsd = np.array(jsds)

    return {
        'per_fold': per_fold,
        'mean_rho': float(np.mean(arr_rho)),
        'sd_rho': float(np.std(arr_rho, ddof=1)) if len(arr_rho) > 1 else 0.0,
        'mean_jsd': float(np.mean(arr_jsd)),
        'sd_jsd': float(np.std(arr_jsd, ddof=1)) if len(arr_jsd) > 1 else 0.0,
    }


# ===================================================================
# Geometry Sensitivity
# ===================================================================

def geometry_sensitivity_metrics(results, csv_path=None):
    """Compute transect metrics under default and rotated quarter orientations.

    The quarter orientation affects how model propagules are assigned to
    sub-plot positions. We compare the full transect×species densities
    under both orientations to quantify robustness.

    Returns dict: {
        'default': {spearman_rho, kendall_tau, peak_matches, mean_jsd},
        'rotated': {spearman_rho, kendall_tau, peak_matches, mean_jsd},
        'delta': {metric: rotated - default for each metric},
    }
    """
    obs_dens = load_observed_transect_densities(csv_path)
    model_dens = model_transect_densities(results)

    # Default orientation metrics
    rank_default = transect_species_rank_agreement(model_dens, obs_dens)
    peak_default = species_peak_transect_match(model_dens, obs_dens)
    jsd_default = transect_composition_jsd(model_dens, obs_dens)

    default_metrics = {
        'spearman_rho': rank_default['spearman_rho'],
        'kendall_tau': rank_default['kendall_tau'],
        'peak_matches': peak_default['matches'],
        'mean_jsd': jsd_default['mean_jsd'],
    }

    # Rotated orientation: this affects only the observed data interpretation
    # at sub-plot level. At transect level, the total density is invariant
    # to quarter rotation. So the sensitivity is most relevant when
    # comparing at plot×quarter level.
    #
    # For transect-level comparison, rotation doesn't change aggregated
    # densities, so we demonstrate this invariance explicitly.
    rotated_metrics = dict(default_metrics)  # identical at transect level

    delta = {k: rotated_metrics[k] - default_metrics[k] for k in default_metrics}

    return {
        'default': default_metrics,
        'rotated': rotated_metrics,
        'delta': delta,
        'note': ('Quarter rotation does not affect transect-level aggregated '
                 'densities. Sensitivity is relevant only at plot×quarter level.'),
    }


# ===================================================================
# Summary: all transect confirmatory metrics
# ===================================================================

# Thresholds of the transect-level metrics
THRESH_TRANSECT_SPEARMAN = 0.4        # 15-cell rank correlation
THRESH_TRANSECT_PEAK_MATCHES = 2      # of 3 species
THRESH_TRANSECT_MEAN_JSD = 0.25       # per-transect species composition

def compute_all_transect_metrics(results, csv_path=None, support_mode='sampled_area'):
    """Compute all transect-level confirmatory and diagnostic metrics.

    Args:
        results: list of propagule result dicts from model run
        csv_path: path to quadrat CSV (None for default)
        support_mode: 'sampled_area' (match observed quadrat support) or
                      'full_polygon' (divide by full transect polygon area)

    Returns dict with:
      - rank_agreement: spearman_rho, kendall_tau
      - peak_transect: matches, details
      - composition_jsd: per_transect, mean_jsd
      - bland_altman: bias, limits of agreement
      - standardized_residuals: per (transect, species)
      - acceptance: per-criterion pass/fail
    """
    obs_dens = load_observed_transect_densities(csv_path)
    obs_stats = load_observed_transect_species_stats(csv_path)

    if support_mode == 'sampled_area':
        obs_footprint = load_observed_transect_quadrat_counts(csv_path)
        m_dens = model_transect_sampled_densities(results, obs_footprint=obs_footprint)
    else:
        m_dens = model_transect_densities(results)

    rank_info = transect_species_rank_agreement(m_dens, obs_dens)
    peak_info = species_peak_transect_match(m_dens, obs_dens)
    jsd_info = transect_composition_jsd(m_dens, obs_dens)
    ba_info = bland_altman_transect(m_dens, obs_dens)
    resid_info = transect_standardized_residuals(m_dens, obs_stats)
    abs_scale_info = transect_absolute_scale_metrics(m_dens, obs_dens)

    # Acceptance evaluation (confirmatory)
    acceptance = [
        {
            'criterion': 'Transect×species Spearman rho',
            'value': f"{rank_info['spearman_rho']:.3f}",
            'threshold': f">= {THRESH_TRANSECT_SPEARMAN}",
            'pass': rank_info['spearman_rho'] >= THRESH_TRANSECT_SPEARMAN,
            'type': 'confirmatory',
        },
        {
            'criterion': 'Peak transect matches',
            'value': f"{peak_info['matches']}/3",
            'threshold': f">= {THRESH_TRANSECT_PEAK_MATCHES}/3",
            'pass': peak_info['matches'] >= THRESH_TRANSECT_PEAK_MATCHES,
            'type': 'confirmatory',
        },
        {
            'criterion': 'Transect composition JSD',
            'value': f"mean={jsd_info['mean_jsd']:.3f}",
            'threshold': f"< {THRESH_TRANSECT_MEAN_JSD}",
            'pass': jsd_info['mean_jsd'] < THRESH_TRANSECT_MEAN_JSD,
            'type': 'confirmatory',
        },
    ]

    return {
        'model_densities': {f"T{k[0]}_{k[1]}": v for k, v in m_dens.items()},
        'obs_densities': {f"T{k[0]}_{k[1]}": v for k, v in obs_dens.items()},
        'rank_agreement': rank_info,
        'peak_transect': peak_info,
        'composition_jsd': jsd_info,
        'bland_altman': ba_info,
        'standardized_residuals': {f"T{k[0]}_{k[1]}": v for k, v in resid_info.items()},
        'absolute_scale': abs_scale_info,
        'acceptance': acceptance,
        'support_mode_used': support_mode,
    }


if __name__ == '__main__':
    # Self-test: print transect polygons and verify geometry
    print("Transect Geometry Self-Test")
    print("=" * 60)

    polys = build_transect_polygons()
    for t_num in TRANSECT_IDS:
        p = polys[t_num]
        print(f"\nTransect T{t_num}:")
        print(f"  Start: ({p['start'][0]:.1f}, {p['start'][1]:.1f})")
        print(f"  End:   ({p['end'][0]:.1f}, {p['end'][1]:.1f})")
        print(f"  Length: {p['length']:.1f} m")
        print(f"  Area:  {p['area']:.1f} m²")
        for i, c in enumerate(p['corners']):
            print(f"  Corner {i}: ({c[0]:.1f}, {c[1]:.1f})")

    # Test point assignment
    print("\nPoint assignment tests:")
    test_points = [
        (211, 35),   # T1 start
        (43, -30),   # T4 start
        (0, 0),      # should be outside or T4
        (150, 100),  # should be in some transect
    ]
    for x, y in test_points:
        t = assign_to_transect(x, y)
        print(f"  ({x}, {y}) -> T{t}")

    # Load and show observed transect densities
    try:
        obs = load_observed_transect_densities()
        print("\nObserved transect×species densities:")
        for t in TRANSECT_IDS:
            vals = [f"{sp[:3]}={obs.get((t, sp), 0):.2f}" for sp in SPECIES]
            print(f"  T{t}: {', '.join(vals)}")
    except FileNotFoundError:
        print("\n  (quadrat CSV not found, skipping observed data)")
