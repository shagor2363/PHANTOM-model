"""
Tide-mediated propagule dispersal and establishment model
=========================================================

Individual-based Lagrangian model of mangrove propagule transport and post-stranding
establishment on a planar intertidal gradient. The ODD description is in CLAUDE.md; the
provenance and calibration history of every constant is in docs/parameter_provenance.md.

Equations
  Eq. 1    z(x, y) = Z_SEA + (Z_LAND - Z_SEA) d / D                        bed elevation
           d = (x - X_MIN) sin(a) + (y - Y_MIN) cos(a), the distance along the unit vector
           GRAD_DIR (angle a = GRAD_ANGLE_DEG from +y toward +x) from the seaward corner;
           D = DIAG, the value of d at the landward corner
  Eq. 4    h(t) = H0 + A(t) sin(2 pi t / T + phi)                          water level
           A(t): cosine envelope between A_NEAP and A_SPRING over the spring-neap cycle
  Eq. T1   per tidal cycle: landward m F tau + N(0, 0.3 m F), lateral N(0, sigma F)
           (m mobility, F zone drag factor, tau tidal scale); the cycle displacement is
           spread over the 2 h substeps and applied only while the propagule floats
  Eq. T2   P_settle = 1 / (1 + exp(-k (z - z_crit))), z_crit = H0 + c_A A - c_d draft,
           evaluated at the end of each cycle from cycle MIN_FLOAT_STEPS + 1 on
  Eq. 2-3  seeds per tree ~ Poisson(F_base (r_stem / r_min)^gamma)         fecundity
  Eq. 8    burial = max(0, R_base (H_wet / 24) T_factor eps - erosion)      daily siltation
           eps ~ lognormal(0, 0.3) clipped to [0.5, 2]
  Eq. 9    P_pred = min(MAX_DAILY_P, 1 - exp(-lambda D V E A H))           daily predation
           H = max(0.1, 1 - h / H_REF), the height term
  Eq. 10   rooting on forest floor only (no rooting on the sand band)       substrate

Eqs. T1 and T2 replace the hydrodynamic ODD Eqs. 5-7 by a biased random walk with a
logistic settlement rule keyed to propagule draft.

Scenarios: (a) full_mechanistic = boundary influx + parent trees, vegetation drag;
(b) internal_only = parent trees only, drag; (c) boundary_only = boundary influx, no drag.

Every constant is a module attribute read at call time, so experiments override it by
assignment (see transport_treatments.run_with_overrides). Species-keyed dictionaries
must contain every species in SPECIES.
"""

import argparse
import csv
import math
import os
from collections import Counter

import numpy as np

# ===================================================================
# 1. Parameters
# ===================================================================

# Species. The first three were calibrated on the 1999 juvenile zonation; Sonneratia alba
# (added 2026-09-20) is parameterised by rule and its placement is a prediction to report.
SPECIES = ('Avicennia', 'Ceriops', 'Rhizophora', 'Sonneratia')
CALIBRATED_SPECIES = SPECIES[:3]

# --- 1a. Domain (ODD 1.2.4): the Gazi transect map (Dahdouh-Guebas et al. 2004, Fig. 6)
# rotated so the parallel transects run along +y; x alongshore, y seaward -> landward.
X_MIN, X_MAX = 0.0, 325.0        # m
Y_MIN, Y_MAX = 0.0, 225.0        # m
Z_SEA, Z_LAND = 1.0, 2.2         # m, bed elevation at the seaward and landward corners
GRAD_ANGLE_DEG = 15.0            # deg, angle a of the gradient from +y toward +x (orientation of the belt relative to the transects)
GRAD_DIR = np.array([math.sin(math.radians(GRAD_ANGLE_DEG)),
                     math.cos(math.radians(GRAD_ANGLE_DEG))])   # unit vector seaward -> landward (flood direction)
DIAG = (X_MAX - X_MIN) * GRAD_DIR[0] + (Y_MAX - Y_MIN) * GRAD_DIR[1]   # m, gradient extent D in Eq. 1
ZONE_B_MIN, ZONE_A_MIN = 1.4, 1.8   # m, zone boundaries C|B and B|A (thirds of the gradient)
RELEASE_DEPTH = 55.0             # m along the gradient: boundary propagules enter where d < RELEASE_DEPTH

# --- 1b. Tidal forcing (Eq. 4)
H0 = 1.6                         # m, mean water level
A_SPRING, A_NEAP = 1.25, 0.5     # m, tidal amplitude at spring and neap
T_PERIOD = 12.42                 # h, M2 period
SN_CYCLE = 15                    # days, spring-neap cycle: days 1-7 spring, 8-15 neap
SN_SPRING_CENTER = 4.0           # day, peak of the cosine amplitude envelope
PHASE_OFFSET = -4.50             # rad, high water at noon on day 0

# --- 1c. Transport (Eqs. T1, T2). Status: [lit] literature, [cal] calibrated on the 1999
# zonation, [asm] assumed (Sonneratia borrows Avicennia's transport values).
MOBILITY = {'Avicennia': 55.0, 'Ceriops': 45.0, 'Rhizophora': 40.0, 'Sonneratia': 55.0}   # m per cycle, net landward displacement [cal; Son asm]
SIGMA_LATERAL = {'Avicennia': 12.0, 'Ceriops': 8.0, 'Rhizophora': 5.0, 'Sonneratia': 12.0}  # m per cycle, lateral spread [lit; Son asm]
DRAFT = {'Avicennia': 0.03, 'Ceriops': 0.12, 'Rhizophora': 0.175, 'Sonneratia': 0.02}      # m, submerged depth afloat [lit mid-range; Son asm]
TIDE_SCALE = {'spring': 1.0, 'neap': 0.4}      # displacement scaling by tidal state (amplitude ratio)
N_CYCLES = {'spring': 8, 'neap': 4}            # tidal cycles in the transport window [cal]
MIN_FLOAT_STEPS = {'spring': 1, 'neap': 1}     # settlement (Eq. T2) first evaluated at the end of cycle MIN_FLOAT_STEPS + 1 [cal]
SUBSTEP_DT = 2.0                               # h, transport substep (six per cycle)
GROUNDING_MARGIN = 0.02                        # m, depth buffer above the draft (ODD Eq. 7, delta_g)
SETTLE_STEEPNESS = 2.0                         # 1/m, logistic steepness k (Eq. T2)
SETTLE_TIDAL_COEFF = 0.6                       # c_A: spring amplitude raises z_crit (longer transport)
SETTLE_DRAFT_COEFF = 6.0                       # c_d: draft lowers z_crit (earlier settlement) [cal]
ZONE_DRAG = {'C': 0.4, 'B': 0.7, 'A': 1.0}     # movement retention factor F (1 = no attenuation), strongest drag seaward [cal]

# --- 1d. Establishment (Eqs. 8-10)
R_BASE = {'C': 0.65, 'B': 0.3, 'A': 0.1}       # mm/day, siltation at full inundation [lit; C cal]
T_FACTOR = {'spring': 1.3, 'neap': 0.7}        # tidal modifier of siltation
EROSION = {'C': 0.005, 'B': 0.003, 'A': 0.002}  # mm/day, ebb erosion
BURIAL_TOL = {'Avicennia': 30, 'Ceriops': 50, 'Rhizophora': 80, 'Sonneratia': 20}         # mm, lethal cumulative burial [lit-inferred]
LAMBDA_PRED = 0.003                            # base encounter rate (Clarke & Kerrigan 2002) [cal]
PRED_D = {'C': 24.0, 'B': 7.97, 'A': 4.24}     # effective predator density per zone [C cal]
VULN = {'Avicennia': 0.85, 'Ceriops': 0.55, 'Rhizophora': 0.25, 'Sonneratia': 0.85}        # species vulnerability [lit ordering; Son asm]
TIDE_E = {'spring': 0.8, 'neap': 1.2}          # tidal modifier of predation
AGE_A = [(14, 1.5), (30, 1.0), (90, 0.7)]      # (last day, modifier): age modifier of predation
H_REF = {'Avicennia': 10.0, 'Ceriops': 28.0, 'Rhizophora': 30.0, 'Sonneratia': 8.0}       # cm, reference height of the height term H = max(0.1, 1 - h / H_REF); growth stops at 2 H_REF [Son asm]
MAX_DAILY_P = 0.15                             # cap on the daily predation probability
INIT_HEIGHT = {'Avicennia': 1.4, 'Ceriops': 8.25, 'Rhizophora': 18.0, 'Sonneratia': 2.0}  # cm, seedling height at rooting
DAILY_GROWTH = {'Avicennia': 0.3, 'Ceriops': 0.22, 'Rhizophora': 0.25, 'Sonneratia': 0.05}  # cm/day
# Substrate (Eq. 10): two classes, forest floor (rooting) and the treeless sand band seaward
# of the forest edge (no rooting), 1 m cells from reconstruction/make_substrate_map.py.
SUBSTRATE_MAP_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'results', 'reconstruction', 'substrate_map.npz')
SUBSTRATE_SUITABILITY = {'forest': 1.0, 'sand': 0.0}   # rooting probability per substrate class

# --- 1e. Propagule supply (Eqs. 2-3 and boundary influx)
F_BASE = {'Avicennia': 12.8, 'Ceriops': 48.0, 'Rhizophora': 8.0, 'Sonneratia': 0.00384}   # seeds per year of a tree at r_min [cal; Son fitted supply scale]
GAMMA_FECUNDITY = {'Avicennia': 2.0, 'Ceriops': 2.5, 'Rhizophora': 2.5, 'Sonneratia': 2.0}  # allometric exponent [Son asm]
MATURITY = {'Avicennia': 0.027, 'Ceriops': 0.027, 'Rhizophora': 0.024, 'Sonneratia': 0.025}  # m, minimum reproductive r_stem
BOUNDARY_SPRING = {'Avicennia': 448, 'Ceriops': 384, 'Rhizophora': 160, 'Sonneratia': 0.1344}  # propagules/day released in the seaward strip (d < RELEASE_DEPTH) at peak season, spring days; x seasonal_scale otherwise [cal; Son fitted]
BOUNDARY_NEAP = {'Avicennia': 224, 'Ceriops': 192, 'Rhizophora': 64, 'Sonneratia': 0.0672}     # same at neap
# Fruiting season (shoulder start, peak start, peak end, shoulder end) as day of year; Gazi Bay
# phenology (Wang'ondu et al. 2013; Ochieng & Erftemeijer 2002; Slim et al. 1996).
FRUITING_WINDOW = {
    'Avicennia':  (60, 91, 152, 182),    # Mar 1 - Jul 1, peak Apr - May
    'Ceriops':    (335, 1, 59, 90),      # Dec 1 - Mar 31, peak Jan - Feb (wraps the year end)
    'Rhizophora': (152, 244, 304, 335),  # Jun 1 - Dec 1, peak Sep - Oct
    'Sonneratia': (213, 244, 319, 349),  # Aug 1 - Dec 15, peak Sep 1 - Nov 15
}

# Transport states of a propagule (outcomes in the results use the strings 'established',
# 'dead', 'lost').
FLOATING, STRANDED, DEAD = 0, 1, 2


# ===================================================================
# 2. Domain and tidal physics (Eqs. 1, 4)
# ===================================================================

def along_gradient(x, y):
    """Distance d (m) along GRAD_DIR from the seaward corner; not clamped."""
    return (x - X_MIN) * GRAD_DIR[0] + (y - Y_MIN) * GRAD_DIR[1]


def elevation(x, y):
    """Eq. 1: planar bed elevation (m); positions are clamped to the domain box."""
    xc = max(X_MIN, min(X_MAX, x))
    yc = max(Y_MIN, min(Y_MAX, y))
    return Z_SEA + (Z_LAND - Z_SEA) * along_gradient(xc, yc) / DIAG


def zone_at(x, y):
    """Elevation zone: 'C' (seaward, z < 1.4 m), 'B' (middle) or 'A' (landward, z >= 1.8 m)."""
    z = elevation(x, y)
    if z < ZONE_B_MIN:
        return 'C'
    if z < ZONE_A_MIN:
        return 'B'
    return 'A'


def zone_boundary_lines(elevations=(ZONE_B_MIN, ZONE_A_MIN)):
    """Iso-elevation segments ((x0, y0), (x1, y1)) clipped to the domain box, for drawing."""
    segs = []
    gx, gy = GRAD_DIR
    for z in elevations:
        d = (z - Z_SEA) / (Z_LAND - Z_SEA) * DIAG
        pts = []
        if abs(gy) > 1e-12:
            for x in (X_MIN, X_MAX):
                y = Y_MIN + (d - (x - X_MIN) * gx) / gy
                if Y_MIN - 1e-9 <= y <= Y_MAX + 1e-9:
                    pts.append((x, y))
        if abs(gx) > 1e-12:
            for y in (Y_MIN, Y_MAX):
                x = X_MIN + (d - (y - Y_MIN) * gy) / gx
                if X_MIN - 1e-9 <= x <= X_MAX + 1e-9:
                    pts.append((x, y))
        pts = sorted(set((round(a, 6), round(b, 6)) for a, b in pts))
        if len(pts) >= 2:
            segs.append((pts[0], pts[-1]))
    return segs


def _is_spring(day):
    """Days 1-7 of each 15-day cycle are spring, 8-15 neap (day <= 0 counts as spring)."""
    return day <= 0 or (day - 1) % SN_CYCLE < 7


def tidal_state(day):
    """'spring' or 'neap' for a 1-indexed simulation day."""
    return 'spring' if _is_spring(day) else 'neap'


def amplitude(t_hours):
    """Tidal amplitude A(t): cosine envelope between A_NEAP and A_SPRING over the spring-neap cycle."""
    day = t_hours / 24.0
    a_mean = 0.5 * (A_SPRING + A_NEAP)
    a_range = 0.5 * (A_SPRING - A_NEAP)
    return a_mean + a_range * math.cos(2.0 * math.pi * (day - SN_SPRING_CENTER) / SN_CYCLE)


def water_level(t_hours):
    """Eq. 4: water surface elevation h(t) in m."""
    return H0 + amplitude(t_hours) * math.sin(2.0 * math.pi * t_hours / T_PERIOD + PHASE_OFFSET)


def water_depth(x, y, t_hours):
    """Local water depth max(0, h(t) - z) in m."""
    return max(0.0, water_level(t_hours) - elevation(x, y))


def hours_wet_on_day(x, y, day):
    """Hours of inundation (> 1 cm) at (x, y) on *day*, sampled hourly."""
    return float(sum(1 for h in range(24) if water_depth(x, y, day * 24.0 + h) > 0.01))


# ===================================================================
# 3. Transport (Eqs. T1, T2)
# ===================================================================

def settlement_probability(species, x, y, tstate):
    """Eq. T2: logistic settlement probability of a floating propagule at (x, y)."""
    a = A_SPRING if tstate == 'spring' else A_NEAP
    z_crit = H0 + SETTLE_TIDAL_COEFF * a - SETTLE_DRAFT_COEFF * DRAFT[species]
    return 1.0 / (1.0 + math.exp(-SETTLE_STEEPNESS * (elevation(x, y) - z_crit)))


def _is_deep_enough(x, y, t_hours, species):
    """True when the water depth exceeds the species draft plus the grounding margin."""
    return water_depth(x, y, t_hours) > DRAFT[species] + GROUNDING_MARGIN


def simulate_propagule(x0, y0, species, release_day, use_drag, rng, release_t_hours=None):
    """Eq. T1: transport one propagule until it strands, is washed out or its window expires.

    Each tidal cycle is split into substeps of SUBSTEP_DT hours; the propagule moves only
    while the water is deeper than its draft (float when wet, pause when dry). Stranding
    occurs by the settlement rule (Eq. T2, evaluated at the end of each cycle from cycle
    MIN_FLOAT_STEPS + 1 on), by reaching the
    landward edge ('boundary_clamp') or when the transport window expires ('forced');
    crossing the seaward edge washes the propagule out ('washed_out').

    Returns (x, y, state, settle_method, t_ground, x_ground, y_ground), the last three
    being the time (h) and position of the first contact with the bed (None if washed out
    before touching it).
    """
    tstate = tidal_state(release_day)
    mob, sigma_lat, tide = MOBILITY[species], SIGMA_LATERAL[species], TIDE_SCALE[tstate]
    n_sub = max(1, round(T_PERIOD / SUBSTEP_DT))
    lat_dir = (-GRAD_DIR[1], GRAD_DIR[0])

    x, y = float(x0), float(y0)
    t_hours = release_t_hours if release_t_hours is not None else release_day * 24.0
    state, settle_method = FLOATING, 'forced'
    grounding = None   # (t, x, y) of the first contact with the bed

    for cycle in range(N_CYCLES[tstate]):
        drag = ZONE_DRAG[zone_at(x, y)] if use_drag else 1.0
        landward = mob * tide * drag + rng.normal(0, mob * 0.3 * drag)
        lateral = rng.normal(0, sigma_lat * drag)
        dx = (landward * GRAD_DIR[0] + lateral * lat_dir[0]) / n_sub
        dy = (landward * GRAD_DIR[1] + lateral * lat_dir[1]) / n_sub

        for _ in range(n_sub):
            t_hours += SUBSTEP_DT
            if not _is_deep_enough(x, y, t_hours, species):        # dry: pause in place
                grounding = grounding or (t_hours, x, y)
                continue
            nx = max(X_MIN, min(X_MAX, x + dx))
            ny = y + dy
            if ny < Y_MIN:                                          # carried past the seaward edge
                state, settle_method = DEAD, 'washed_out'
                break
            if ny > Y_MAX:                                          # reached the landward edge
                x, y = nx, Y_MAX
                state, settle_method = STRANDED, 'boundary_clamp'
                grounding = grounding or (t_hours, x, y)
                break
            if not _is_deep_enough(nx, ny, t_hours, species):      # destination too shallow: wait
                grounding = grounding or (t_hours, nx, ny)
                continue
            x, y = nx, ny

        if state != FLOATING:
            break
        if cycle >= MIN_FLOAT_STEPS[tstate] and rng.random() < settlement_probability(species, x, y, tstate):
            state, settle_method = STRANDED, 'probability'
            grounding = grounding or (t_hours, x, y)
            break

    if state == FLOATING:                                           # transport window expired
        state = STRANDED
        grounding = grounding or (t_hours, x, y)

    t_g, x_g, y_g = grounding if grounding else (None, None, None)
    return x, y, state, settle_method, t_g, x_g, y_g


# ===================================================================
# 4. Propagule supply (Eqs. 2-3, boundary influx, phenology)
# ===================================================================

def find_release_hour(x, y, day, species):
    """First hour of *day* at which (x, y) is deep enough for the species to float; None if never."""
    draft_req = DRAFT[species] + GROUNDING_MARGIN
    for h in range(24):
        t = day * 24.0 + h
        if water_depth(x, y, t) > draft_req:
            return t
    return None


def boundary_release_position(rng):
    """Uniform random position in the seaward strip d < RELEASE_DEPTH (parallel to the zone boundaries)."""
    while True:
        x = rng.uniform(X_MIN, X_MAX)
        y = rng.uniform(Y_MIN, Y_MAX)
        if along_gradient(x, y) < RELEASE_DEPTH:
            return x, y


def seasonal_scale(species, day):
    """Trapezoidal fruiting-season weight: 0 off season, 1 at peak, linear over the shoulders."""
    s_start, p_start, p_end, s_end = FRUITING_WINDOW[species]
    doy = (day - 1) % 365 + 1
    if s_start > s_end:                 # window wraps the year end: shift it to start on day 1
        shift = s_start - 1
        doy = (doy - shift - 1) % 365 + 1
        s_end = (s_end - shift - 1) % 365 + 1
        p_start = (p_start - shift - 1) % 365 + 1
        p_end = (p_end - shift - 1) % 365 + 1
        s_start = 1
    if p_start <= doy <= p_end:
        return 1.0
    if s_start <= doy < p_start:
        return (doy - s_start) / max(1, p_start - s_start)
    if p_end < doy <= s_end:
        return 1.0 - (doy - p_end) / max(1, s_end - p_end)
    return 0.0


def _fruiting_days(species):
    """Days of year inside the fruiting window, shoulder to shoulder."""
    s_start, _, _, s_end = FRUITING_WINDOW[species]
    if s_start <= s_end:
        return list(range(s_start, s_end + 1))
    return list(range(s_start, 366)) + list(range(1, s_end + 1))


def schedule_annual_seeds(trees, species, rng, fbase_override=None):
    """Eqs. 2-3: one year of propagules of *species* from the parent trees.

    Each tree with r_stem >= MATURITY produces Poisson(F_base (r_stem / r_min)^gamma)
    propagules at its position; they are spread over the fruiting window with weights
    seasonal_scale(). Returns {day_of_year: [(x, y), ...]}.
    """
    r_min = MATURITY[species]
    f_base = (fbase_override or F_BASE)[species]
    gamma = GAMMA_FECUNDITY[species]

    seeds = []
    for t in trees:
        if t['species'] != species or t['r_stem'] < r_min:
            continue
        n = rng.poisson(f_base * (t['r_stem'] / r_min) ** gamma)
        sx = max(X_MIN, min(X_MAX, t['x']))
        sy = max(Y_MIN, min(Y_MAX, t['y']))
        seeds.extend([(sx, sy)] * n)
    if not seeds:
        return {}

    days = _fruiting_days(species)
    weights = np.array([seasonal_scale(species, d) for d in days])
    weights /= weights.sum()
    assigned = rng.choice(days, size=len(seeds), p=weights)
    schedule = {}
    for pos, day in zip(seeds, assigned):
        schedule.setdefault(day, []).append(pos)
    return schedule


def load_trees(csv_path):
    """Tree inventory (columns plant, species, x, y, r_stem); species outside SPECIES are skipped."""
    trees = []
    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            if row['species'] in SPECIES:
                trees.append({'plant': row['plant'], 'species': row['species'],
                              'x': float(row['x']), 'y': float(row['y']),
                              'r_stem': float(row['r_stem'])})
    return trees


# ===================================================================
# 5. Establishment (Eqs. 8-10)
# ===================================================================

_SUBSTRATE = None   # lazily loaded substrate map


def _load_substrate():
    global _SUBSTRATE
    if _SUBSTRATE is None:
        if os.path.exists(SUBSTRATE_MAP_FILE):
            d = np.load(SUBSTRATE_MAP_FILE)
            _SUBSTRATE = {'group': d['group'], 'names': [str(n) for n in d['names']],
                          'x_min': float(d['x_min']), 'y_min': float(d['y_min']),
                          'cell': float(d['cell'])}
        else:
            _SUBSTRATE = {}
            print(f"  Warning: substrate map not found ({SUBSTRATE_MAP_FILE}); substrate rule inactive")
    return _SUBSTRATE


def substrate_at(x, y):
    """Substrate class at (x, y): 'forest' or 'sand' ('forest' everywhere without a map)."""
    sub = _load_substrate()
    if not sub:
        return 'forest'
    g = sub['group']
    i = min(max(int((y - sub['y_min']) / sub['cell']), 0), g.shape[0] - 1)
    j = min(max(int((x - sub['x_min']) / sub['cell']), 0), g.shape[1] - 1)
    return sub['names'][int(g[i, j])]


def daily_siltation(zone, hours_wet, tstate, rng):
    """Eq. 8: net sediment deposition (mm) on one day; eps ~ lognormal(0, 0.3) clipped to [0.5, 2]."""
    wet_frac = max(0.0, min(1.0, hours_wet / 24.0))
    noise = float(np.clip(rng.lognormal(0, 0.3), 0.5, 2.0))
    gross = R_BASE[zone] * wet_frac * T_FACTOR[tstate] * noise
    return max(0.0, gross - EROSION[zone])


def daily_predation_prob(zone, species, height_cm, day_idx, tstate):
    """Eq. 9: daily predation probability, capped at MAX_DAILY_P."""
    age_mod = next((mod for last_day, mod in AGE_A if day_idx <= last_day), AGE_A[-1][1])
    protection = max(0.1, min(1.0, 1.0 - height_cm / H_REF[species]))
    raw = 1.0 - math.exp(-LAMBDA_PRED * PRED_D[zone] * VULN[species] * TIDE_E[tstate]
                         * age_mod * protection)
    return max(0.0, min(MAX_DAILY_P, raw))


def evaluate_establishment(species, x, y, strand_day, rng):
    """90-day establishment filter at the stranding position.

    *strand_day* is the release day as passed by run_scenario (transport lasts at most
    N_CYCLES cycles, i.e. up to 4 days, so the window is counted from the release).
    Day 0: rooting according to the substrate (Eq. 10). Then daily: growth, siltation
    (Eq. 8, lethal when the cumulative burial reaches BURIAL_TOL) and predation (Eq. 9).
    Returns {outcome, cause, survival_days, burial_mm}.
    """
    zone = zone_at(x, y)
    suit = SUBSTRATE_SUITABILITY[substrate_at(x, y)]
    if suit < 1.0 and rng.random() >= suit:
        return {'outcome': 'dead', 'cause': 'substrate', 'survival_days': 0, 'burial_mm': 0.0}

    height = INIT_HEIGHT[species]
    burial = 0.0
    for d in range(1, 91):
        abs_day = strand_day + d
        tstate = tidal_state(int(abs_day))
        height = min(height + DAILY_GROWTH[species], 2.0 * H_REF[species])
        burial += daily_siltation(zone, hours_wet_on_day(x, y, abs_day), tstate, rng)
        if burial >= BURIAL_TOL[species]:
            return {'outcome': 'dead', 'cause': 'siltation', 'survival_days': d, 'burial_mm': burial}
        if rng.random() < daily_predation_prob(zone, species, height, d, tstate):
            return {'outcome': 'dead', 'cause': 'predation', 'survival_days': d, 'burial_mm': burial}
    return {'outcome': 'established', 'cause': None, 'survival_days': 90, 'burial_mm': burial}


# ===================================================================
# 6. Scenario runner and output
# ===================================================================

# key -> (name, boundary influx, parent trees, vegetation drag, default tree inventory)
SCENARIOS = {
    'a': ('full_mechanistic', True, True, True, 'demo.csv'),
    'b': ('internal_only', False, True, True, 'demo.csv'),
    'c': ('boundary_only', True, False, False, None),
}

# seed-source mode -> (boundary influx, parent trees)
SEED_SOURCE_MODES = {
    'combined': (True, True),
    'boundary_only': (True, False),
    'internal_only': (False, True),
}


def resolve_seed_source_mode(scenario_key, seed_source_mode=None):
    """(use_boundary, use_internal, mode name) from the scenario defaults or an explicit mode."""
    if seed_source_mode is not None:
        if seed_source_mode not in SEED_SOURCE_MODES:
            raise ValueError(f"Unknown seed_source_mode '{seed_source_mode}'. "
                             f"Valid modes: {list(SEED_SOURCE_MODES)}")
        use_boundary, use_internal = SEED_SOURCE_MODES[seed_source_mode]
        return use_boundary, use_internal, seed_source_mode
    _, use_boundary, use_internal, _, _ = SCENARIOS[scenario_key]
    for mode, flags in SEED_SOURCE_MODES.items():
        if flags == (use_boundary, use_internal):
            return use_boundary, use_internal, mode
    return use_boundary, use_internal, 'combined'


def stochastic_round(lam, rng):
    """Unbiased integer draw for an expected count: floor(lam) + Bernoulli(lam - floor(lam))."""
    if lam <= 0:
        return 0
    base = int(lam)
    return base + (1 if rng.random() < (lam - base) else 0)


def _release_generators(seed):
    """One release-count generator per species, spawned from the run seed.

    Keeping these draws out of the main generator means a change in a release rate cannot
    shift the transport and establishment draws of the other species; spawning from a
    SeedSequence makes them reproducible across processes.
    """
    children = np.random.SeedSequence(seed).spawn(len(SPECIES))
    return {sp: np.random.default_rng(child) for sp, child in zip(SPECIES, children)}


def _track_propagule(sp, x, y, day, source, use_drag, rng, ledger, results):
    """Release one propagule at (x, y) on *day*, transport it and apply the establishment
    filter; append one result row (none when the site never floods that day)."""
    led = ledger[sp]
    led['produced'] += 1
    rel_t = find_release_hour(x, y, day, sp)
    if rel_t is None:
        led['not_released'] += 1
        return
    led['released'] += 1
    fx, fy, fstate, settle_method, t_g, x_g, y_g = simulate_propagule(
        x, y, sp, day, use_drag, rng, release_t_hours=rel_t)
    led['settle_' + settle_method] += 1
    settle_day = int(t_g / 24.0) + 1 if t_g is not None else day
    row = {
        'species': sp, 'x': fx, 'y': fy, 'zone': zone_at(fx, fy), 'source': source, 'day': day,
        'settle_method': settle_method,
        'first_grounding_t': t_g, 'first_grounding_x': x_g, 'first_grounding_y': y_g,
        'first_grounding_z': elevation(x_g, y_g) if x_g is not None else None,
        'terminal_z': elevation(fx, fy),
        'tidal_state_at_settle': tidal_state(settle_day),
    }
    if fstate == DEAD:
        led['lost'] += 1
        row.update(outcome='lost', cause='washed_out', survival_days=0, burial_mm=0.0)
    else:
        est = evaluate_establishment(sp, fx, fy, day, rng)
        led['established' if est['outcome'] == 'established' else 'dead_' + est['cause']] += 1
        row.update(est)
    results.append(row)


def run_scenario(scenario_key, n_days=365, seed=42, csv_dir=None,
                 boundary_spring_override=None, boundary_neap_override=None,
                 fbase_override=None, seed_source_mode=None, return_ledger=False,
                 tree_csv=None):
    """Run a scenario for *n_days* and return one result dict per released propagule.

    Result keys: species, x, y, zone, source ('boundary' | 'internal'), day, outcome
    ('established' | 'dead' | 'lost'), cause, survival_days, burial_mm, settle_method,
    first_grounding_{t,x,y,z}, terminal_z, tidal_state_at_settle.

    csv_dir / tree_csv: directory of the scenario's default inventory, or an explicit
    inventory path. The override dicts replace BOUNDARY_SPRING, BOUNDARY_NEAP and F_BASE
    for this run; seed_source_mode overrides the scenario's source flags (drag stays per
    scenario). With return_ledger the per-species stage counts are returned as well.
    """
    name, _, _, use_drag, csv_file = SCENARIOS[scenario_key]
    if tree_csv is not None and csv_file is not None:
        csv_file, csv_dir = os.path.abspath(tree_csv), None
    use_boundary, use_internal, mode = resolve_seed_source_mode(scenario_key, seed_source_mode)
    rng = np.random.default_rng(seed)
    release_rng = _release_generators(seed)
    print(f"\n--- Scenario ({scenario_key}): {name} ---")
    print(f"    seed_source_mode={mode}, boundary={use_boundary}, internal={use_internal}, "
          f"drag={use_drag}, days={n_days}")

    trees = []
    if csv_file is not None:
        path = os.path.join(csv_dir, csv_file) if csv_dir else csv_file
        if os.path.exists(path):
            trees = load_trees(path)
            print(f"    Loaded {len(trees)} trees from {path}")
        elif use_internal:
            raise FileNotFoundError(f"Tree file '{path}' not found but seed_source_mode='{mode}' "
                                    f"requires internal seeds.")
        else:
            print(f"    Warning: {path} not found, running without internal trees")
    elif use_internal:
        raise FileNotFoundError(f"No tree inventory for scenario '{scenario_key}' but "
                                f"seed_source_mode='{mode}' requires internal seeds.")

    results = []
    ledger = {sp: {'produced': 0, 'released': 0, 'not_released': 0,
                   'settle_probability': 0, 'settle_forced': 0,
                   'settle_boundary_clamp': 0, 'settle_washed_out': 0,
                   'established': 0, 'dead_siltation': 0, 'dead_predation': 0,
                   'dead_substrate': 0, 'lost': 0} for sp in SPECIES}

    internal_schedule = {}
    if use_internal and trees:
        internal_schedule = {sp: schedule_annual_seeds(trees, sp, rng, fbase_override=fbase_override)
                             for sp in SPECIES}

    for day in range(1, n_days + 1):
        if tidal_state(day) == 'spring':
            rates = boundary_spring_override or BOUNDARY_SPRING
        else:
            rates = boundary_neap_override or BOUNDARY_NEAP
        doy = (day - 1) % 365 + 1
        for sp in SPECIES:
            if use_boundary:
                lam = rates[sp] * seasonal_scale(sp, day)
                for _ in range(stochastic_round(lam, release_rng[sp])):
                    x, y = boundary_release_position(rng)
                    _track_propagule(sp, x, y, day, 'boundary', use_drag, rng, ledger, results)
            if use_internal:
                for x, y in internal_schedule.get(sp, {}).get(doy, []):
                    _track_propagule(sp, x, y, day, 'internal', use_drag, rng, ledger, results)

        if day % 30 == 0:
            n_est = sum(1 for r in results if r['outcome'] == 'established')
            print(f"    Day {day:>4d}/{n_days}: {len(results)} propagules, {n_est} established")

    n_est = sum(1 for r in results if r['outcome'] == 'established')
    print(f"    Done: {len(results)} total, {n_est} established")
    return (results, ledger) if return_ledger else results


def summarize(results):
    """Print the species x zone table of established recruits and the mortality by zone."""
    zones = ('A', 'B', 'C')
    est = [r for r in results if r['outcome'] == 'established']
    dead = [r for r in results if r['outcome'] == 'dead']
    lost = [r for r in results if r['outcome'] == 'lost']

    print(f"\n  Total propagules: {len(results)}")
    print(f"  Established: {len(est)}, Dead: {len(dead)}, Lost: {len(lost)}")
    causes = Counter(r['cause'] for r in dead)
    if causes:
        print(f"  Mortality causes: {dict(causes)}")

    print(f"\n  {'Species':<14s}" + ''.join(f"  Zone {z}" for z in zones) + "   Total")
    print("  " + "-" * 48)
    for sp in SPECIES:
        counts = [sum(1 for r in est if r['species'] == sp and r['zone'] == z) for z in zones]
        print(f"  {sp:<14s}" + ''.join(f"  {n:>5d}" for n in counts) + f"  {sum(counts):>6d}")
    print("  " + "-" * 48)
    totals = [sum(1 for r in est if r['zone'] == z) for z in zones]
    print(f"  {'Total':<14s}" + ''.join(f"  {n:>5d}" for n in totals) + f"  {sum(totals):>6d}")

    print("\n  Mortality by zone (%):")
    for z in zones:
        n_zone = sum(1 for r in results if r['zone'] == z and r['outcome'] in ('established', 'dead'))
        n_dead = sum(1 for r in dead if r['zone'] == z)
        pct = 100.0 * n_dead / n_zone if n_zone else 0
        print(f"    Zone {z}: {pct:.1f}% ({n_dead}/{n_zone})")


CSV_FIELDS = ['species', 'x', 'y', 'zone', 'source', 'day', 'outcome', 'cause',
              'survival_days', 'burial_mm', 'settle_method',
              'first_grounding_t', 'first_grounding_x', 'first_grounding_y',
              'first_grounding_z', 'terminal_z', 'tidal_state_at_settle']


def write_csv(results, path):
    """Write one row per propagule."""
    if not results:
        return
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(results)
    print(f"  Results written to {path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tide-mediated propagule dispersal model')
    parser.add_argument('--scenario', '-s', default='c', choices=['a', 'b', 'c', 'all'],
                        help='a = full, b = parent trees only, c = boundary only, all')
    parser.add_argument('--days', '-d', type=int, default=365, help='simulation days')
    parser.add_argument('--seed', type=int, default=42, help='random seed')
    parser.add_argument('--output', '-o', default='results', help='output directory for the CSVs')
    parser.add_argument('--seed-source-mode', default=None, choices=list(SEED_SOURCE_MODES),
                        help='override the scenario source flags')
    parser.add_argument('--trees', default=None,
                        help='tree inventory CSV (default: demo.csv next to this script)')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    for key in (['a', 'b', 'c'] if args.scenario == 'all' else [args.scenario]):
        results = run_scenario(key, n_days=args.days, seed=args.seed, csv_dir=here,
                               seed_source_mode=args.seed_source_mode, tree_csv=args.trees)
        write_csv(results, os.path.join(args.output, f'{SCENARIOS[key][0]}.csv'))
        summarize(results)
