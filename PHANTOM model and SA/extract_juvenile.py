"""
Extract the observed juvenile (seedling) records of one survey year from the Gazi Bay field
workbook, transects 1-5. The survey year is the command-line argument:
    python extract_juvenile.py 1997
    python extract_juvenile.py 1999

WHAT A ROW IS (data audit 2026-09-20, docs/juvenile_data_audit.md)
==================================================================
Source:      GAZ1_6 (1).XLS, sheet 'GAZ1 env.summary' (Dahdouh-Guebas et al. 2004, MEPS 272)
Row:         one quadrant (a-d) of one PCQM sample point ('plot'); the juvenile densities are
             recorded per quadrant (they are NOT copied across the four quadrants of a point)
Year scope:  juvenile density columns of YEAR only (1997: cols 11-18, 1999: cols 19-26; four
             species/density pairs, header '#/m2')
Transects:   1, 2, 3, 4, 5 written out; transect 6 is read for the duplicate check only
Species:     all five recorded species: Avicennia marina ('A. mar'), Rhizophora mucronata
             ('R. muc' / 'R.muc'), Ceriops tagal ('C. tag'), Sonneratia alba ('S. alb'),
             Bruguiera gymnorhiza ('B. gym'). Which of them the model uses is decided
             downstream by tidal_dispersal_model.SPECIES.
Row filter:  Tree type = 'A' (adult PCQM records, which carry the juvenile densities)
Quadrat ID:  Parsed as {transect}.{plot}{letter}; non-numeric IDs and summary rows excluded
Zone map:    Veg class -> model zone per ODD Section 1.2.4 (the position-based re-zoning of
             reconstruction/zone_quadrats_by_position.py overwrites model_zone afterwards and
             keeps this class-based zone in veg_zone):
               'R. muc'                -> A (landward, 1.8-2.2 m)
               'A. mar', 'Mixed'       -> B (middle, 1.4-1.8 m)
               'Open space', 'S. alb'  -> C (seaward, 1.0-1.4 m)

THE OBSERVATION IS THE RECORDED DENSITY; AREAS AND COUNTS ARE HYPOTHESES
=======================================================================
`density_per_m2` is the recorded value ('#/m2'), written unchanged. The workbook and the paper
state the quadrat area of the juvenile survey nowhere: the paper (p. 81) gives a 5 m x 5 m
quadrat per quadrant for vegetation cover and 1 m2 quadrats for crabs and snails, and says
only that 'relative density and relative frequency were also calculated for young and
juvenile mangrove trees'. Rounded densities do not determine the quadrat area or the number of
individuals. The earlier column `count_per_quadrat` (= density x 2.5 m2, from the unsupported
reading '1 plot = 10 m2, 4 quadrats per plot') produced fractional pseudo-counts (0.1, 0.275)
and was REMOVED (2026-09-20). Two clearly labelled hypothesis columns replace it; they are used
nowhere in calibration or testing:

  unit_area_m2_hypothesis   1997: 12.5 for every quadrat row (all 17 distinct 1997 densities are
                            multiples of 0.08 = 1/12.5, consistent with, not proof of, a 12.5 m2
                            unit and single individuals; the uniform-unit hypothesis is applied
                            to the whole year).
                            1999: 25, 9 or 4 when EVERY non-zero density of the quadrat row
                            (all species slots) is an integer multiple of exactly one of 1/25,
                            1/9, 1/4 (compared after rounding to 2 decimals, as stored,
                            tolerance 1e-6); NaN when the values fit none, more than one
                            family (integers fit all three), or the row has no seedlings.
  count_hypothesis          density x unit_area_m2_hypothesis, rounded to the nearest integer
                            when the product lies within the rounding tolerance of the stored
                            2-decimal density (0.005 x area + 1e-6), else NaN; 0 for a zero
                            density (no seedlings under any area); NaN when the area is NaN.

DUPLICATE PLOT BLOCKS (flagged, never deleted)
==============================================
Whole plots (all four quadrants, all species slots of YEAR) are compared across transects 1-6.
Two plots with identical juvenile blocks (at least one non-zero density) are 'suspected'
transcription duplicates; if their ADULT measurement blocks (cols 5-10: species, d, GBH, DBH,
height, CD of all quarters) are identical as well, 'strongly_suspected'. Both members carry
the flag, each other's plot id (`duplicate_pair`, e.g. 'T2.7') and `duplicate_order` (1 =
first occurrence in the sheet order, 2 = second). Nothing is removed from the main data;
sensitivity variants are built downstream.

ZONE SUMMARIES
==============
Per zone x species: n_quadrats, mean_density (mean of the recorded per-quadrat densities), its
SD, SE and t-based 95 % CI, and `sum_of_densities_not_a_total` (the plain sum of the recorded
densities, kept for reconciliation only: it is NOT a number of seedlings, the units have
unknown and possibly unequal areas). Zone shares used downstream are RELATIVE DENSITY ACROSS
ZONES, p_z = mean density in zone z / sum over zones of the mean densities (not the fraction of
individuals per zone, which would need the zone areas).

OUTPUT FILES (suffix = YEAR)
============================
  juvenile_quadrat_records_{YEAR}_t1_5.csv        Row-level recorded density + hypothesis and flag columns
  juvenile_zone_species_summary_{YEAR}_t1_5.csv   Zone x species statistics
  juvenile_zone_totals_{YEAR}_t1_5.csv            Zonewise all-species statistics
  juvenile_qc_report_{YEAR}_t1_5.csv              QC / reconciliation checks, unit families, duplicates
  results/extraction_audit_{YEAR}.csv             rows excluded by the tree-type filter

Script version: 2.0 (2026-09-20: recorded densities only, hypothesis columns, duplicate flags,
Sonneratia and Bruguiera kept)
"""

import os
import sys
import csv
import math
import re
from datetime import datetime
from collections import defaultdict

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
YEAR = sys.argv[1] if len(sys.argv) > 1 else '1999'   # survey year: python extract_juvenile.py 1997 | 1999

EXCEL_FILE = 'GAZ1_6 (1).XLS'
SHEET_NAME = 'GAZ1 env.summary'
HEADER_ROWS = 4  # rows 0-3 are headers/units; data starts at row 4

# Per-year column layout and unit-area hypothesis rule
YEAR_CONFIG = {
    '1997': dict(pairs=[(11, 12), (13, 14), (15, 16), (17, 18)], header_col=11,
                 header='Juveniles 1997', unit_rule='uniform', unit_area=12.5),
    '1999': dict(pairs=[(19, 20), (21, 22), (23, 24), (25, 26)], header_col=19,
                 header='Juveniles 1999', unit_rule='family', families=(25.0, 9.0, 4.0)),
}
CFG = YEAR_CONFIG[YEAR]
JUV_PAIRS = CFG['pairs']                  # (species_col, density_col) of the year
EXPECTED_HEADER_0 = CFG['header']
EXPECTED_HEADER_1_DENSITY = '#/m2'
ADULT_COLS = (5, 6, 7, 8, 9, 10)          # species, d, GBH, DBH, height, CD (adult PCQM record)

VALID_TRANSECTS = {1, 2, 3, 4, 5}         # written out
DUPLICATE_CHECK_TRANSECTS = {1, 2, 3, 4, 5, 6}   # compared for identical plot blocks

# Vegetation class -> model zone (ODD Section 1.2.4)
VEG_TO_ZONE = {
    'R. muc':     'A',
    'A. mar':     'B',
    'Mixed':      'B',
    'Open space': 'C',
    'S. alb':     'C',
}

# Juvenile species code -> model species name
SPECIES_MAP = {
    'A. mar': 'Avicennia',
    'R. muc': 'Rhizophora',
    'R.muc':  'Rhizophora',
    'C. tag': 'Ceriops',
    'S. alb': 'Sonneratia',   # added 2026-09-20 (11 records 1997, 3 records 1999)
    'B. gym': 'Bruguiera',    # added 2026-09-20 (5 records each year)
}

# Species written to the CSVs (one row per quadrat x species). The model's own species list
# (tidal_dispersal_model.SPECIES) selects the subset that is simulated and evaluated.
MODEL_SPECIES = ['Avicennia', 'Rhizophora', 'Ceriops', 'Sonneratia', 'Bruguiera']

# Quadrat ID regex: transect.plotLetter (e.g. 1.1b, 5.26c)
QUADRAT_RE = re.compile(r'^(\d+)\.(\d+)([a-z])$')

NAN = float('nan')

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val):
    """Convert a cell value to float; return None if not numeric."""
    if val is None:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _cell_str(val):
    """Normalize a cell to stripped string; NaN/None -> ''."""
    if val is None:
        return ''
    s = str(val).strip()
    if s.lower() == 'nan':
        return ''
    return s


def _norm_cell(val):
    """Cell value normalised for block comparison: floats rounded, strings stripped."""
    f = _safe_float(val)
    if f is not None:
        return round(f, 6)
    return _cell_str(val)


def _t_critical_95(n):
    """Approximate two-tailed t critical value for 95% CI (lookup; 1.96 for large n)."""
    table = {
        2: 12.706, 3: 4.303, 4: 3.182, 5: 2.776, 6: 2.571,
        7: 2.447, 8: 2.365, 9: 2.306, 10: 2.262, 11: 2.228,
        12: 2.201, 13: 2.179, 14: 2.160, 15: 2.145, 16: 2.131,
        17: 2.120, 18: 2.110, 19: 2.101, 20: 2.093, 25: 2.064,
        30: 2.045, 40: 2.021, 50: 2.009, 60: 2.000, 80: 1.990,
        100: 1.984, 120: 1.980,
    }
    if n <= 1:
        return NAN
    df = n - 1
    if df in table:
        return table[df]
    for k in sorted(table.keys()):
        if k - 1 >= df:
            return table[k]
    return 1.96


def _fmt(v):
    """CSV cell for a float that may be NaN ('' for NaN)."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ''
    return v


# ---------------------------------------------------------------------------
# Unit-area and count hypotheses
# ---------------------------------------------------------------------------

def _fits_unit(d, unit):
    """True when the stored 2-decimal density d is an integer multiple (>= 1) of `unit`,
    compared after rounding both to 2 decimals (tolerance 1e-6)."""
    k = round(d / unit)
    return k >= 1 and abs(round(d, 2) - round(k * unit, 2)) < 1e-6


def unit_area_hypothesis(nonzero_densities):
    """Unit area (m2) hypothesised for one quadrat row from its non-zero densities, and a label
    of the family test: ('uniform' | '/25' | '/9' | '/4' | 'all' | 'two' | 'none' | 'zero').
    1997 (uniform rule): 12.5 for every row. 1999 (family rule): exactly one of 25/9/4."""
    if CFG['unit_rule'] == 'uniform':
        return CFG['unit_area'], ('zero' if not nonzero_densities else 'uniform')
    if not nonzero_densities:
        return NAN, 'zero'
    fits = [a for a in CFG['families'] if all(_fits_unit(d, 1.0 / a) for d in nonzero_densities)]
    if len(fits) == 1:
        return fits[0], f'/{int(fits[0])}'
    if len(fits) == 0:
        return NAN, 'none'
    return NAN, ('all' if len(fits) == len(CFG['families']) else 'two')


def count_hypothesis(density, area):
    """Integer count implied by density x area under the rounding tolerance of the stored
    2-decimal density; 0 for a zero density; NaN when the area is NaN or the product is not
    an integer within tolerance."""
    if density == 0:
        return 0
    if area is None or (isinstance(area, float) and math.isnan(area)):
        return NAN
    prod = density * area
    k = round(prod)
    tol = 0.005 * area + 1e-6
    return int(k) if abs(prod - k) <= tol else NAN


# ---------------------------------------------------------------------------
# Schema validation
# ---------------------------------------------------------------------------

def validate_schema(df_head):
    """Check that the expected juvenile headers of YEAR are present."""
    msgs = []
    ok = True
    hc = CFG['header_col']
    h0 = _cell_str(df_head[0][hc]) if len(df_head[0]) > hc else ''
    if h0 != EXPECTED_HEADER_0:
        msgs.append(f"FAIL: Row 0, Col {hc} expected '{EXPECTED_HEADER_0}', got '{h0}'")
        ok = False
    else:
        msgs.append(f"OK: Row 0, Col {hc} = '{h0}'")
    for sp_col, dens_col in JUV_PAIRS:
        h1_sp = _cell_str(df_head[1][sp_col]) if len(df_head[1]) > sp_col else ''
        h1_d = _cell_str(df_head[1][dens_col]) if len(df_head[1]) > dens_col else ''
        if h1_sp != 'Species':
            msgs.append(f"FAIL: Row 1, Col {sp_col} expected 'Species', got '{h1_sp}'")
            ok = False
        if h1_d != EXPECTED_HEADER_1_DENSITY:
            msgs.append(f"FAIL: Row 1, Col {dens_col} expected '{EXPECTED_HEADER_1_DENSITY}', got '{h1_d}'")
            ok = False
    if ok:
        msgs.append(f"OK: All {YEAR} juvenile column headers validated")
    return ok, msgs


# ---------------------------------------------------------------------------
# Quadrat parsing
# ---------------------------------------------------------------------------

def parse_quadrat(quad_str):
    """Parse quadrat ID into (transect: int, plot: int, letter: str) or None."""
    m = QUADRAT_RE.match(quad_str)
    if m:
        return int(m.group(1)), int(m.group(2)), m.group(3)
    return None


# ---------------------------------------------------------------------------
# Duplicate plot blocks
# ---------------------------------------------------------------------------

def find_duplicate_plots(plot_blocks):
    """plot_blocks: {(transect, plot): {'first_row': int, 'juv': {letter: block}, 'adult': {letter: block}}}
    Returns {(transect, plot): {'status', 'pair', 'order', 'members'}} for every plot whose
    juvenile block (all quadrants, all slots) equals that of another plot and holds at least one
    non-zero density; 'strongly_suspected' when the adult block is identical to that of at least
    one partner as well."""
    groups = defaultdict(list)
    for key, b in plot_blocks.items():
        juv_sig = tuple(sorted((letter, blk) for letter, blk in b['juv'].items()))
        has_data = any(d > 0 for blk in b['juv'].values() for _sp, d in blk)
        if has_data:
            groups[juv_sig].append(key)
    flags = {}
    for sig, members in groups.items():
        if len(members) < 2:
            continue
        members = sorted(members, key=lambda k: plot_blocks[k]['first_row'])
        for i, key in enumerate(members):
            others = [m for m in members if m != key]
            adult_same = any(plot_blocks[key]['adult'] == plot_blocks[m]['adult'] for m in others)
            flags[key] = {
                'status': 'strongly_suspected' if adult_same else 'suspected',
                'pair': ';'.join(f'T{t}.{p}' for t, p in others),
                'order': i + 1,
                'members': members,
            }
    return flags


# ---------------------------------------------------------------------------
# Main extraction
# ---------------------------------------------------------------------------

def extract(excel_path):
    """Run full extraction pipeline.

    Returns (records, qc) where records is the list of row-level dicts (transects 1-5, one
    per quadrat x species) and qc the audit counters, unit-family counts and duplicate flags.
    """
    try:
        import xlrd
    except ImportError:
        sys.exit("ERROR: xlrd is required.  pip install xlrd")

    if not os.path.exists(excel_path):
        sys.exit(f"ERROR: File not found: {excel_path}")

    wb = xlrd.open_workbook(excel_path)
    if SHEET_NAME not in wb.sheet_names():
        sys.exit(f"ERROR: Sheet '{SHEET_NAME}' not found.  Available: {wb.sheet_names()}")

    sh = wb.sheet_by_name(SHEET_NAME)
    print(f"Opened sheet '{SHEET_NAME}': {sh.nrows} rows x {sh.ncols} cols")

    header_rows_data = [[sh.cell_value(r, c) for c in range(sh.ncols)] for r in range(min(HEADER_ROWS, sh.nrows))]
    schema_ok, schema_msgs = validate_schema(header_rows_data)
    for m in schema_msgs:
        print(f"  Schema: {m}")
    if not schema_ok:
        sys.exit(f"ABORT: Schema validation failed: column layout does not match the expected {YEAR} structure.")

    qc = {
        'total_data_rows': 0,
        'rows_not_adult': 0,
        'rows_no_quadrat': 0,
        'rows_malformed_quadrat': 0,
        'rows_transect_6_excluded': 0,
        'rows_transect_other_excluded': 0,
        'rows_unknown_veg_class': 0,
        'rows_summary_excluded': 0,
        'rows_included': 0,
        'species_codes_unknown': set(),
        'species_codes_nonmodel': set(),
        'negative_densities': 0,
        'non_numeric_densities': 0,
        'missing_densities': 0,
        'unit_family_counts': defaultdict(int),
        'unit_family_counts_t6': defaultdict(int),
        'duplicates': {},
        'messages': list(schema_msgs),
    }
    summary_veg_labels = {'Mean', 'Median', 'St.dev.', 'St. dev.'}

    plot_blocks = {}     # (transect, plot) -> blocks of every quadrant, transects 1-6
    quadrat_rows = []    # per included quadrat: dict with the parsed row (transects 1-5)

    for r in range(HEADER_ROWS, sh.nrows):
        qc['total_data_rows'] += 1
        tree_type = _cell_str(sh.cell_value(r, 0))
        if tree_type != 'A':
            qc['rows_not_adult'] += 1
            continue
        quad_raw = _cell_str(sh.cell_value(r, 1))
        if not quad_raw:
            qc['rows_no_quadrat'] += 1
            continue
        veg_class = _cell_str(sh.cell_value(r, 2))
        if veg_class in summary_veg_labels:
            qc['rows_summary_excluded'] += 1
            continue
        parsed = parse_quadrat(quad_raw)
        if parsed is None:
            qc['rows_malformed_quadrat'] += 1
            continue
        transect, plot, letter = parsed

        # raw juvenile slots of YEAR and the adult block, for the duplicate check (T1-T6)
        juv_block = tuple((_cell_str(sh.cell_value(r, sc)), round(_safe_float(sh.cell_value(r, dc)) or 0.0, 6))
                          for sc, dc in JUV_PAIRS)
        adult_block = tuple(_norm_cell(sh.cell_value(r, c)) for c in ADULT_COLS)
        if transect in DUPLICATE_CHECK_TRANSECTS:
            pb = plot_blocks.setdefault((transect, plot), {'first_row': r, 'juv': {}, 'adult': {}})
            pb['juv'][letter] = juv_block
            pb['adult'][letter] = adult_block

        if transect == 6:
            qc['rows_transect_6_excluded'] += 1
            continue
        if transect not in VALID_TRANSECTS:
            qc['rows_transect_other_excluded'] += 1
            continue
        zone = VEG_TO_ZONE.get(veg_class)
        if zone is None:
            qc['rows_unknown_veg_class'] += 1
            continue
        qc['rows_included'] += 1

        species_densities = {}
        for sp_col, dens_col in JUV_PAIRS:
            sp_raw = _cell_str(sh.cell_value(r, sp_col))
            dens_raw = sh.cell_value(r, dens_col)
            if not sp_raw or sp_raw in ('0', '0.0'):     # empty slot (Excel stores 0 as float 0.0)
                continue
            model_sp = SPECIES_MAP.get(sp_raw)
            if model_sp is None:
                qc['species_codes_unknown'].add(sp_raw)
                continue
            d = _safe_float(dens_raw)
            if d is None:
                qc['non_numeric_densities'] += 1
                d = 0.0
            elif d < 0:
                qc['negative_densities'] += 1
            species_densities[model_sp] = species_densities.get(model_sp, 0.0) + d

        nonzero = [d for d in species_densities.values() if d > 0]
        area, family = unit_area_hypothesis(nonzero)
        qc['unit_family_counts'][family] += 1
        quadrat_rows.append(dict(transect=transect, plot=plot, letter=letter, quad_raw=quad_raw,
                                 veg_class=veg_class, zone=zone, species_densities=species_densities,
                                 area=area, family=family))

    # unit families with transect 6 included (for the audit document only)
    for (t, p), b in plot_blocks.items():
        if t == 6:
            for letter, blk in b['juv'].items():
                nz = [d for sp, d in blk if sp and sp not in ('0', '0.0') and d > 0]
                qc['unit_family_counts_t6'][unit_area_hypothesis(nz)[1]] += 1
    for k, v in qc['unit_family_counts'].items():
        qc['unit_family_counts_t6'][k] += v

    dup = find_duplicate_plots(plot_blocks)
    qc['duplicates'] = dup

    records = []
    for q in quadrat_rows:
        flag = dup.get((q['transect'], q['plot']))
        for sp in MODEL_SPECIES:
            dens = q['species_densities'].get(sp, 0.0)
            is_missing = sp not in q['species_densities']
            records.append({
                'transect': q['transect'],
                'plot': q['plot'],
                'quadrat_letter': q['letter'],
                'quadrat_id': q['quad_raw'],
                'veg_class': q['veg_class'],
                'model_zone': q['zone'],
                'species': sp,
                'density_per_m2': round(dens, 6),
                'is_missing': is_missing,
                'unit_area_m2_hypothesis': q['area'],
                'count_hypothesis': count_hypothesis(round(dens, 6), q['area']),
                'duplicate_status': flag['status'] if flag else '',
                'duplicate_pair': flag['pair'] if flag else '',
                'duplicate_order': flag['order'] if flag else '',
            })
            if is_missing:
                qc['missing_densities'] += 1
    return records, qc


# ---------------------------------------------------------------------------
# Summary statistics (recorded densities only)
# ---------------------------------------------------------------------------

def _density_stats(densities):
    n = len(densities)
    mean_d = sum(densities) / n
    var_d = sum((d - mean_d) ** 2 for d in densities) / (n - 1) if n > 1 else 0.0
    sd_d = math.sqrt(var_d)
    se_d = sd_d / math.sqrt(n) if n > 0 else 0.0
    t_crit = _t_critical_95(n)
    ci_lo = mean_d - t_crit * se_d if not math.isnan(t_crit) else NAN
    ci_hi = mean_d + t_crit * se_d if not math.isnan(t_crit) else NAN
    return {
        'n_quadrats': n,
        'mean_density': round(mean_d, 6),
        'sd_density': round(sd_d, 6),
        'se_density': round(se_d, 6),
        'ci95_lo': round(ci_lo, 6) if not math.isnan(ci_lo) else '',
        'ci95_hi': round(ci_hi, 6) if not math.isnan(ci_hi) else '',
        'sum_of_densities_not_a_total': round(sum(densities), 6),
    }


def zone_species_summary(records):
    """Zone x species statistics of the recorded per-quadrat densities: n_quadrats,
    mean_density (with SD, SE, t-based 95 % CI) and sum_of_densities_not_a_total (kept for
    reconciliation only; not a number of seedlings). Zone shares downstream are relative
    densities across zones: p_z = mean density in zone z / sum over zones of the means."""
    groups = defaultdict(list)
    for rec in records:
        groups[(rec['model_zone'], rec['species'])].append(float(rec['density_per_m2']))
    summary = []
    for (zone, species) in sorted(groups.keys()):
        row = {'zone': zone, 'species': species}
        row.update(_density_stats(groups[(zone, species)]))
        summary.append(row)
    return summary


def zone_totals(records):
    """Zonewise statistics of the all-species density per quadrat (sum over the species of a
    quadrat's recorded densities, then the same statistics as zone_species_summary)."""
    quad_totals = defaultdict(float)
    zone_quads = defaultdict(set)
    for rec in records:
        key = (rec['model_zone'], rec['quadrat_id'])
        quad_totals[key] += float(rec['density_per_m2'])
        zone_quads[rec['model_zone']].add(rec['quadrat_id'])
    results = []
    for zone in ['A', 'B', 'C']:
        quads = sorted(zone_quads.get(zone, set()))
        if not quads:
            continue
        row = {'zone': zone}
        row.update(_density_stats([quad_totals[(zone, q)] for q in quads]))
        results.append(row)
    return results


# ---------------------------------------------------------------------------
# QC / audit report
# ---------------------------------------------------------------------------

def build_qc_rows(records, qc, zone_sp_summary, zone_tot):
    """Build QC report rows as list of dicts {check, result, detail}."""
    rows = []

    def add(check, result, detail=''):
        rows.append({'check': check, 'result': result, 'detail': detail})

    # --- Data scope ---
    add('survey_year', YEAR)
    add('total_data_rows', qc['total_data_rows'])
    add('rows_not_adult', qc['rows_not_adult'])
    add('rows_no_quadrat', qc['rows_no_quadrat'])
    add('rows_malformed_quadrat', qc['rows_malformed_quadrat'])
    add('rows_summary_excluded', qc['rows_summary_excluded'])
    add('rows_transect_6_excluded', qc['rows_transect_6_excluded'],
        'transect 6 rows are read for the duplicate-plot check only')
    add('rows_transect_other_excluded', qc['rows_transect_other_excluded'])
    add('rows_unknown_veg_class', qc['rows_unknown_veg_class'])
    add('rows_included', qc['rows_included'],
        f'= {qc["rows_included"]} quadrats x {len(MODEL_SPECIES)} species = '
        f'{qc["rows_included"] * len(MODEL_SPECIES)} records')

    t6_in_records = sum(1 for r in records if r['transect'] == 6)
    add('VERIFY_no_transect_6_in_output', 'PASS' if t6_in_records == 0 else 'FAIL',
        f'{t6_in_records} records from transect 6')

    add('nonmodel_species_codes', ', '.join(sorted(qc['species_codes_nonmodel'])) or 'none')
    add('unknown_species_codes', ', '.join(sorted(qc['species_codes_unknown'])) or 'none')
    add('non_numeric_densities', qc['non_numeric_densities'])
    add('negative_densities', qc['negative_densities'])
    add('missing_densities_set_to_zero', qc['missing_densities'],
        'species not listed in quadrat row; density set to 0')

    # Reconciliation: the zone total mean equals the sum of the species means (same quadrats)
    for zt in zone_tot:
        zone = zt['zone']
        sp_sum = sum(s['mean_density'] for s in zone_sp_summary if s['zone'] == zone)
        match = abs(zt['mean_density'] - sp_sum) < 1e-4
        add(f'RECONCILE_zone_{zone}_mean_total_eq_sum_of_species_means',
            'PASS' if match else 'FAIL',
            f"zone_mean_total={zt['mean_density']:.6f}, sum_species_means={sp_sum:.6f}")

    expected_records = qc['rows_included'] * len(MODEL_SPECIES)
    add('VERIFY_record_count', 'PASS' if len(records) == expected_records else 'FAIL',
        f'expected={expected_records}, actual={len(records)}')

    # --- The observation and the hypotheses ---
    add('observation', 'density_per_m2 as recorded (#/m2)',
        'the recorded density is the observation; no count is derived for analysis')
    add('quadrat_area', 'not documented; hypotheses only',
        'paper p. 81: 5 x 5 m quadrat per quadrant for vegetation cover, 1 m2 for crabs; juvenile '
        'quadrat size stated nowhere; unit_area_m2_hypothesis / count_hypothesis are hypotheses')
    if CFG['unit_rule'] == 'uniform':
        add('unit_area_hypothesis_rule', f'uniform {CFG["unit_area"]} m2 for every quadrat row',
            'all non-zero densities are multiples of 0.08 = 1/12.5 (consistent with, not proof of)')
    else:
        add('unit_area_hypothesis_rule',
            'exactly one of ' + ', '.join(f'1/{int(a)}' for a in CFG['families']) + ' fits every non-zero density of the row',
            'NaN when none, more than one (integers fit all) or no seedlings')
    fam = qc['unit_family_counts']
    order = ['uniform', '/25', '/9', '/4', 'all', 'two', 'none', 'zero']
    for k in order:
        if k in fam:
            label = {'uniform': 'unit_family_uniform_rows', 'all': 'unit_family_all_three_rows',
                     'two': 'unit_family_two_families_rows', 'none': 'unit_family_none_rows',
                     'zero': 'unit_family_no_seedlings_rows'}.get(k, f'unit_family_{k}_rows')
            add(label, fam[k], 'quadrat rows of transects 1-5 (all species slots considered)')
    fam6 = qc['unit_family_counts_t6']
    add('unit_family_counts_incl_transect_6',
        '; '.join(f'{k} {fam6[k]}' for k in order if k in fam6), 'transects 1-6, for the audit document')
    n_count_nan = sum(1 for r in records if r['density_per_m2'] > 0 and
                      isinstance(r['count_hypothesis'], float) and math.isnan(r['count_hypothesis']))
    n_count_int = sum(1 for r in records if r['density_per_m2'] > 0 and not
                      (isinstance(r['count_hypothesis'], float) and math.isnan(r['count_hypothesis'])))
    add('count_hypothesis_nonzero_rows_integer', n_count_int)
    add('count_hypothesis_nonzero_rows_undetermined', n_count_nan)

    # --- Duplicate plot blocks ---
    dup = qc['duplicates']
    groups = {}
    for key, f in dup.items():
        groups.setdefault(tuple(f['members']), f['status'])
    add('duplicate_plot_groups', len(groups),
        'identical juvenile blocks (all quadrants, all slots); nothing deleted')
    for members, status in sorted(groups.items(), key=lambda kv: kv[0]):
        strong = any(dup[m]['status'] == 'strongly_suspected' for m in members)
        add('duplicate_pair', ' = '.join(f'T{t}.{p}' for t, p in members),
            ('strongly_suspected (adult block cols 5-10 identical too)' if strong else
             'suspected (juvenile block identical, adults differ)'))

    # Transect / zone distribution
    from collections import Counter
    t_counts = Counter(r['transect'] for r in records)
    for t in sorted(VALID_TRANSECTS):
        add(f'transect_{t}_quadrats', t_counts.get(t, 0) // len(MODEL_SPECIES))
    z_counts = Counter(r['model_zone'] for r in records)
    for z in ['A', 'B', 'C']:
        add(f'zone_{z}_quadrats', z_counts.get(z, 0) // len(MODEL_SPECIES))
    return rows


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def write_csv(rows, fields, path):
    """Write list of dicts to CSV (NaN floats as empty cells)."""
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: _fmt(row.get(k)) for k in fields})
    print(f"  Written: {path} ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Extraction audit: tree_type filter analysis
# ---------------------------------------------------------------------------

def audit_tree_type_filter(excel_path=None):
    """Audit rows excluded by the tree_type != 'A' filter: per tree_type value, whether the
    rows carry quadrat IDs or non-zero densities of YEAR. Writes results/extraction_audit_{YEAR}.csv."""
    try:
        import xlrd
    except ImportError:
        print("  WARNING: xlrd not available; skipping extraction audit.")
        return None
    if excel_path is None:
        excel_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), EXCEL_FILE)
    if not os.path.exists(excel_path):
        print(f"  WARNING: Excel file not found ({excel_path}); skipping extraction audit.")
        return None
    wb = xlrd.open_workbook(excel_path)
    if SHEET_NAME not in wb.sheet_names():
        print(f"  WARNING: Sheet '{SHEET_NAME}' not found; skipping extraction audit.")
        return None
    sh = wb.sheet_by_name(SHEET_NAME)
    groups = defaultdict(lambda: {'count': 0, 'quadrat_ids': 0, 'nonzero_densities': 0})
    for r in range(HEADER_ROWS, sh.nrows):
        tree_type = _cell_str(sh.cell_value(r, 0))
        if tree_type == 'A':
            continue
        g = groups[tree_type if tree_type else '(empty)']
        g['count'] += 1
        if _cell_str(sh.cell_value(r, 1)):
            g['quadrat_ids'] += 1
        for _sp_col, dens_col in JUV_PAIRS:
            d = _safe_float(sh.cell_value(r, dens_col))
            if d is not None and d != 0.0:
                g['nonzero_densities'] += 1
                break
    audit_rows = []
    total_count = total_quadrats = total_densities = 0
    for tree_type in sorted(groups.keys()):
        g = groups[tree_type]
        has_quadrats = g['quadrat_ids'] > 0
        has_densities = g['nonzero_densities'] > 0
        qc_flag = ('potential_data' if has_quadrats and has_densities else
                   'non_sample' if not has_quadrats and not has_densities else 'ambiguous')
        audit_rows.append({'tree_type': tree_type, 'count': g['count'], 'has_quadrats': has_quadrats,
                           'has_densities': has_densities, 'qc_flag': qc_flag})
        total_count += g['count']
        total_quadrats += g['quadrat_ids']
        total_densities += g['nonzero_densities']
    audit_rows.append({'tree_type': 'TOTAL', 'count': total_count, 'has_quadrats': total_quadrats > 0,
                       'has_densities': total_densities > 0,
                       'qc_flag': f'{total_quadrats} rows w/quadrats, {total_densities} rows w/densities'})
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, f'extraction_audit_{YEAR}.csv')
    write_csv(audit_rows, ['tree_type', 'count', 'has_quadrats', 'has_densities', 'qc_flag'], out_path)
    print(f"\n  Extraction Audit: excluded tree_type values:")
    print(f"  {'tree_type':<12s} {'count':>6s} {'quadrats?':>10s} {'densities?':>11s} {'qc_flag':<16s}")
    print(f"  {'-'*58}")
    for row in audit_rows:
        print(f"  {row['tree_type']:<12s} {row['count']:>6d} {str(row['has_quadrats']):>10s} "
              f"{str(row['has_densities']):>11s} {row['qc_flag']:<16s}")
    return audit_rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

REC_FIELDS = ['transect', 'plot', 'quadrat_letter', 'quadrat_id', 'veg_class', 'model_zone', 'species',
              'density_per_m2', 'is_missing', 'unit_area_m2_hypothesis', 'count_hypothesis',
              'duplicate_status', 'duplicate_pair', 'duplicate_order']
ZS_FIELDS = ['zone', 'species', 'n_quadrats', 'mean_density', 'sd_density', 'se_density',
             'ci95_lo', 'ci95_hi', 'sum_of_densities_not_a_total']
ZT_FIELDS = ['zone', 'n_quadrats', 'mean_density', 'sd_density', 'se_density',
             'ci95_lo', 'ci95_hi', 'sum_of_densities_not_a_total']


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    excel_path = os.path.join(script_dir, EXCEL_FILE)

    print(f"{'='*70}")
    print(f"  {YEAR} Juvenile Extraction (transects 1-5; recorded densities)")
    print(f"  Source: {EXCEL_FILE}")
    print(f"  Timestamp: {datetime.now().isoformat(timespec='seconds')}")
    print(f"{'='*70}")

    records, qc = extract(excel_path)
    print(f"\n  Extracted {len(records)} records "
          f"({qc['rows_included']} quadrats x {len(MODEL_SPECIES)} species)")
    records.sort(key=lambda r: (r['transect'], r['plot'], r['quadrat_letter'], r['species']))

    zs_summary = zone_species_summary(records)
    zt_summary = zone_totals(records)
    qc_rows = build_qc_rows(records, qc, zs_summary, zt_summary)
    qc_rows.append({'check': 'script_version', 'result': '2.0', 'detail': 'recorded densities; hypothesis columns; duplicate flags'})
    qc_rows.append({'check': 'workbook_path', 'result': excel_path, 'detail': ''})
    qc_rows.append({'check': 'run_timestamp', 'result': datetime.now().isoformat(timespec='seconds'), 'detail': ''})

    print(f"\nWriting outputs:")
    write_csv(records, REC_FIELDS, os.path.join(script_dir, f'juvenile_quadrat_records_{YEAR}_t1_5.csv'))
    write_csv(zs_summary, ZS_FIELDS, os.path.join(script_dir, f'juvenile_zone_species_summary_{YEAR}_t1_5.csv'))
    write_csv(zt_summary, ZT_FIELDS, os.path.join(script_dir, f'juvenile_zone_totals_{YEAR}_t1_5.csv'))
    write_csv(qc_rows, ['check', 'result', 'detail'], os.path.join(script_dir, f'juvenile_qc_report_{YEAR}_t1_5.csv'))

    print(f"\n{'='*70}")
    print(f"  Zone x Species Summary ({YEAR}, Transects 1-5; recorded densities, #/m2)")
    print(f"{'='*70}")
    print(f"  {'Zone':<6s} {'Species':<14s} {'n':>5s} {'mean_d':>10s} {'SD':>10s} {'SE':>8s} {'95%CI':>20s}")
    print(f"  {'-'*76}")
    for s in zs_summary:
        ci_str = f"[{s['ci95_lo']}, {s['ci95_hi']}]" if s['ci95_lo'] != '' else 'n/a'
        print(f"  {s['zone']:<6s} {s['species']:<14s} {s['n_quadrats']:>5d} {s['mean_density']:>10.4f} "
              f"{s['sd_density']:>10.4f} {s['se_density']:>8.4f} {ci_str:>20s}")
    print(f"\n  {'Zone':<6s} {'n':>5s} {'mean_total_d':>14s} {'SD':>10s} {'SE':>8s}")
    print(f"  {'-'*48}")
    for t in zt_summary:
        print(f"  {t['zone']:<6s} {t['n_quadrats']:>5d} {t['mean_density']:>14.4f} {t['sd_density']:>10.4f} {t['se_density']:>8.4f}")

    print(f"\n  Unit-area families (quadrat rows, transects 1-5): "
          + ', '.join(f'{k} {v}' for k, v in sorted(qc['unit_family_counts'].items())))
    print(f"  Duplicate plot blocks:")
    seen = set()
    for key, f in sorted(qc['duplicates'].items(), key=lambda kv: kv[1]['members']):
        m = tuple(f['members'])
        if m in seen:
            continue
        seen.add(m)
        print(f"    {' = '.join(f'T{t}.{p}' for t, p in m)}: "
              f"{'strongly_suspected' if any(qc['duplicates'][x]['status'] == 'strongly_suspected' for x in m) else 'suspected'}")

    print(f"\n  QC Checks:")
    for q in qc_rows:
        if str(q['result']).startswith(('PASS', 'FAIL')):
            print(f"    [{q['result']}] {q['check']}: {q['detail']}")

    try:
        audit_tree_type_filter(excel_path)
    except Exception as e:
        print(f"  WARNING: Extraction audit failed: {e}")
    print(f"\nDone.")


if __name__ == '__main__':
    main()
