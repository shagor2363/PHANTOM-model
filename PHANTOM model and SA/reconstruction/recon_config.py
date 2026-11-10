"""Configuration for the Gazi Bay whole-domain adult stand reconstruction.

All stages (PCQM preparation, MTPPR tiles, assembly, validation, LOTO) import
this module.  Paths are resolved relative to the project root so the scripts can
be run from anywhere.
"""
import os
import sys

RECON_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(RECON_DIR)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# --- Inputs -------------------------------------------------------------------
XLS_PATH = os.path.join(ROOT, 'GAZ1_6 (1).XLS')
XLS_SHEET = 'GAZ1 env.summary'
XLS_HEADER_ROWS = 4
ENGINE_R = os.path.join(
    ROOT, 'ChrisWudel-Multi-trait-point-pattern-reconstruction-6982c39',
    'func', 'MTPPR-v.2.0.0.R')
RSCRIPT = os.environ.get('RSCRIPT', r'C:\Program Files\R\R-4.5.3\bin\Rscript.exe')
# 2026-09-20 B0/B1: sha256 of the vendored, unmodified engine (upstream commit 6982c39);
# run_mtppr.check_engine() asserts it before every MTPPR job.
ENGINE_SHA256 = '68fedcb2cd41e687280e2abea09c7def3387c8e857028fb336301ab728e323d8'
# 2026-09-20 B1: habitat mask read for the eligible tile area (grid 'treeless', model frame)
SUBSTRATE_NPZ = os.path.join(ROOT, 'results', 'reconstruction', 'substrate_map.npz')

# --- Outputs ------------------------------------------------------------------
OUT_ROOT = os.path.join(ROOT, 'results', 'reconstruction')
GEOM_DIR = os.path.join(OUT_ROOT, 'geom')
MTPPR_DIR = os.path.join(OUT_ROOT, 'mtppr')
LOG_DIR = os.path.join(OUT_ROOT, 'logs')
FIG_DIR = os.path.join(OUT_ROOT, 'figures')
TABLE_DIR = os.path.join(OUT_ROOT, 'tables')

# --- Survey design (PCQM+, Dahdouh-Guebas & Koedam 2006) ------------------------
# Transect layout, plot counts, tape lengths and T4's tape offset live in transect_geometry.py
# (map-digitised from Dahdouh-Guebas et al. 2004, Fig. 6).  Plot k sits at fringe distance
# u = tape_length - PLOT_LEN*(k-1) - PLOT_LEN/2 (plot 1 = landward end).
STAND_TRANSECTS = [1, 2, 3, 4, 6]           # the five parallel transects (T5, orthogonal, is not used)
VALID_TRANSECTS = STAND_TRANSECTS
PLOT_LEN = 10.0                             # m of tape per PCQM point
# 2026-09-20 B0: PCQM search limit.  Dahdouh-Guebas et al. (2004, MEPS 272:77-92, p. 80-81) state
# only "Sampling at 10 m intervals was undertaken using the point-centred quarter method, PCQM
# (Cottam & Curtis 1956) as described by Cintron & Schaeffer Novelli (1984). The individual tree
# closest to the sample point was identified and recorded in each of 4 quadrants"; no maximum
# search distance is given there, and the workbook remarks (GAZ1 env.summary, cols 42-43) carry
# only intersections, GPS fixes, water depths and sky conditions.  R = 20 m is therefore an
# ANALYTICAL ASSUMPTION: empty quarters are right-censored at R ("no tree recorded within R"),
# recorded distances > R are censored at R in the main variant.  Sensitivity variants in
# prep_pcqm.band_density_variants: R = 10 m, and 'recorded distances uncapped' (records > R enter
# at their recorded distance; EMPTY quarters stay censored at R, never an infinite search area).
DMAX = 20.0                                 # assumed search radius R (m): censoring distance of empty quarters
SEARCH_LIMIT_ASSUMED = True                 # True: R is not documented by the field protocol (see above)
SEARCH_LIMIT_VARIANTS = [20.0, 10.0]        # R values tabulated by band_density_variants
BELT_HALF = 10.0                            # half width of the reference belt (m); holds 97 % of PCQM distances
# Reference density mode (lambda_used of each band):
#   'pcqm_mle':  2026-09-20 B0, DEFAULT: censored PCQM maximum-likelihood estimator
#                lambda = 4 n_uncensored / (pi sum_q min(d_q, R)^2) over the band's quarters
#                (uncensored = recorded d <= R; empty quarters and d > R contribute R^2; missing d
#                excluded), plot-level bootstrap CI, times the modelled-species share of the
#                recorded adults (Bruguiera removed from the total).  The reference belts keep the
#                recorded trees only (belts as wide as the trees they hold); the belt area is NOT a
#                density denominator any more; tile counts are pinned (COUNT_CONTROL).
#   'recorded':  earlier stand's convention (Pambala recipe): recorded trees per belt area, no
#                infill (~235-500/ha by band; a count-per-assigned-area reading, abandoned as a
#                density on 2026-09-20, kept selectable for traceability)
#   'pcqm':      belts filled by dart-throwing infill to the Cottam-Curtis density (not used: a null
#                model would enter the reference)
DENSITY_MODE = 'pcqm_mle'
DENSITY_MODES = ('pcqm_mle', 'recorded', 'pcqm')
# 2026-09-20 B0 ADAPTATION "count control": the MTPPR tile count is set to lambda_used x tile area
# through the wrapper's --pin-n / --lambda-target path (add/delete proposals off) instead of the
# engine's emergent count; verified on one tile first (b0_check.py).
# 2026-09-20 B0, RETIRED (user decision after the checks): count control fails criterion (iv) of the
# one-tile check (b0_check/report.md); the code path is kept as a record.  False = emergent count.
COUNT_CONTROL = False
# 2026-09-20 B0, RETIRED: "reference moment rescaling" (reference per-point second-order moments x
# lambda_target / lambda_belt, independent-thinning APPROXIMATION; mtppr_tile.R
# `reference_moment_rescaling`).  NOT ACCEPTED 2026-09-20: fails the synthetic criterion at the target
# density (b0_check/report.md, 'Synthetic-stand test').  Code path kept as a record; never applied.
REF_MOMENT_RESCALING = False
# 2026-09-20 B0 CANDIDATE (not the method yet): SUPERPOSITION.  An explicit ADDITIONAL reconstruction
# assumption beyond MTPPR (leaving the engine unchanged does not make it the original method): per tile,
# k = lambda_target / lambda_belt realisations of MTPPR with the EMERGENT count (no count control, no
# rescaling), independent seeds, superposed; ceil(k) copies with the last copy randomly thinned so that
# the expected count is lambda_target x tile area.  Copy 1 = the standard run (measured trees fixed,
# belts and exclusion sectors forbidden); copies 2..ceil(k) run WITHOUT the measured trees as fixed
# points (they are absent from those copies' output; exclusion sectors stay forbidden), so every
# measured tree appears exactly once.  See superpose.py.  Accepted only if the synthetic criterion of
# b0_check/report.md ('Superposition candidate') is met.
SUPERPOSITION = True
RECORDED_MIN_HALF = 2.5                     # minimum half width of a recorded-tree belt (m)
QUARTER_SIGNS = {                           # (sign on xp, sign on yp); left = negative xp
    'a': (-1, -1),   # seaward-left
    'b': (+1, -1),   # seaward-right
    'c': (-1, +1),   # landward-left
    'd': (+1, +1),   # landward-right
}

# --- Stratification -------------------------------------------------------------
# 2026-09-20 B1: BAND_BREAKS is only the legacy default of the published five-transect stand;
# recon_geometry.band_segments(..., breaks=) and prep_pcqm.run(..., band_breaks=) take the breaks
# explicitly and persist them in prep_meta['settings']['band_breaks'] (validate.py labels from
# there).  The three-transect design (run_stand3T.py) uses BAND_BREAKS_3T.
BAND_BREAKS = [0.0, 40.0, 80.0, 120.0]          # fringe distance u (m); last band open to L_t (legacy default)
BAND_BREAKS_3T = [0.0, 60.0, 120.0]              # three-transect design: composition changes at u ~55-65 and ~115-125 m
BAND_BREAKS_3T_2BAND = [0.0, 80.0]               # two-band sensitivity layout
RETAINED_TRANSECTS_3T = [4, 2, 6]                # T3, T1 withheld; T5 never used (a priori choice, plan 2026-09-20)
MIN_REMNANT = 20.0                               # remnants shorter than this merge into previous band
CORRIDOR_ORDER = [4, 3, 2, 1, 6]                 # west -> east
HARDCORE = 0.5                                   # minimum inter-stem distance for infill (m)

# --- Species ----------------------------------------------------------------------
# Species reconstructed: the four modelled species (2026-09-20: Sonneratia added as the 4th
# modelled species; the earlier statement that S. alba has no recruits was wrong: the workbook
# holds 11 S. alba juvenile records in 1997 and 3 in 1999, <1 % of the juvenile density).
# B. gymnorhiza (7 adults, 5 juvenile records per year) is not modelled; its recorded trees
# still define exclusion sectors (nearest-tree distances of any species).
STAND_SPECIES = ['Avicennia', 'Ceriops', 'Rhizophora', 'Sonneratia']
SPECIES_LEVELS = STAND_SPECIES
ALL_SPECIES = ['Avicennia', 'Bruguiera', 'Ceriops', 'Rhizophora', 'Sonneratia']
SPECIES_MAP = {
    'r.muc': 'Rhizophora', 's.alb': 'Sonneratia', 'c.tag': 'Ceriops',
    'a.mar': 'Avicennia', 'b.gym': 'Bruguiera',
}
PLANT_PREFIX = {
    'Rhizophora': 'Rhizpop', 'Ceriops': 'Ceripop', 'Avicennia': 'Avipop',
    'Sonneratia': 'Sonnpop', 'Bruguiera': 'Brugpop',
}

# --- Marks -------------------------------------------------------------------------
R_STEM_SCALE = 0.45     # r_stem (m) / 0.45 -> ~[0, 1]; raised from 0.4 on 2026-09-20 (largest S. alba r_stem 0.438 m)
H_STEM_SCALE = 17.0     # h_stem (m) / 17 -> ~[0, 1]
# r_root = coeff * DBH_cm ** exponent (m); Pambala methods note; proxies for Ceriops/Sonneratia
ROOT_ALLOM = {
    'Rhizophora': (0.21, 0.68), 'Ceriops': (0.21, 0.68),
    'Avicennia': (0.25, 0.62), 'Sonneratia': (0.25, 0.62),
    'Bruguiera': (0.19, 0.67),
}

# --- MTPPR settings -----------------------------------------------------------------
N_SEEDS = 5
R_MAX = 8.0
DR = 0.25
BW = 0.5
MIN_STEPS = 40000
STEPS_PER_POINT = 200
MAX_STEPS = 400000
CONV_START = 5000
CONV_WINDOW = 5000
CONV_TOL = 2e-6         # relative energy improvement per step over CONV_WINDOW (runs cost seconds)
TRACE_EVERY = 500
N_WORKERS = 12
TILE_MARGIN = 8.0       # m: each tile is reconstructed with this buffer and clipped afterwards (edge artefacts)
ANCHOR_SEED = 42


def ensure_dirs():
    for d in (OUT_ROOT, GEOM_DIR, MTPPR_DIR, LOG_DIR, FIG_DIR, TABLE_DIR):
        os.makedirs(d, exist_ok=True)
