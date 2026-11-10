"""Model runs with changed constants (the transport-trait treatments) and their result tables.

run_with_overrides(overrides, n_days, seed) replaces module-level constants of tidal_dispersal_model
for one run and restores them afterwards. TREATMENTS holds the calibrated model ('baseline') and the
equalised-trait treatments: draft, mobility or lateral spread set to the mean of the calibrated species
for every species, one at a time ('*_eq_mean') or all three ('all_transport_eq'). run_treatments runs
treatment x seed in parallel and stores one table per run as knockout_<treatment>_seed<seed>.csv.
run_grid does the same for a list of (tag, overrides) and also writes the summary of each run
(cascade_metrics.summarise_run) next to its table.

The parent-tree inventory is read from the environment variable KNOCKOUT_TREES (set by the run scripts);
the default is the stand shipped with this package.
"""
import copy
import csv
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

import tidal_dispersal_model as tdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TREE_CSV = os.environ.get('KNOCKOUT_TREES') or os.path.join(
    SCRIPT_DIR, 'results', 'reconstruction', 'stand3T_T4T3T1T6_single', 'demo_reconstructed_stand3T.csv')


def run_with_overrides(overrides, n_days, seed, scenario='a',
                       seed_source_mode=None):
    """Temporarily replace tdm globals, run the scenario, restore originals.

    overrides : dict mapping a tdm module-level attribute name (e.g.
                'DRAFT', 'MOBILITY', 'SIGMA_LATERAL') to its full
                replacement value. Dict-valued attributes must be supplied
                as a complete dict with every species key the model expects.
    """
    snap = {attr: copy.deepcopy(getattr(tdm, attr)) for attr in overrides}
    try:
        for attr, value in overrides.items():
            setattr(tdm, attr, value)
        # the source-neutral control supplies its own tree file through SCENARIOS
        tree_csv = None if 'SCENARIOS' in overrides else TREE_CSV
        return tdm.run_scenario(
            scenario, n_days=n_days, seed=seed, csv_dir=SCRIPT_DIR,
            seed_source_mode=seed_source_mode, tree_csv=tree_csv,
        )
    finally:
        for attr, value in snap.items():
            setattr(tdm, attr, value)


# ---------------------------------------------------------------------------
# Style constants (match existing figure conventions)
# ---------------------------------------------------------------------------
# 2026-09-20: 4 species (Sonneratia). The species list comes from the model (single source of
# truth); Sonneratia's colour '#3B6CC0' was validated against the other three, marker 'D'.
SPECIES_COLORS = {
    'Avicennia':  '#E07A5F',
    'Rhizophora': '#3D8B7A',
    'Ceriops':    '#7B2D8E',
    'Sonneratia': '#3B6CC0',
}
SPECIES_ORDER = list(tdm.SPECIES)
CALIBRATED_SPECIES = list(tdm.CALIBRATED_SPECIES)   # the three species calibrated on the 1999 zonation
SPECIES_SHORT = {'Avicennia': 'Avi', 'Ceriops': 'Cer', 'Rhizophora': 'Rhi', 'Sonneratia': 'Son'}
SPECIES_MARKERS = {'Avicennia': 'o', 'Ceriops': 's', 'Rhizophora': '^', 'Sonneratia': 'D'}

ZONE_BOUNDARIES = [1.4, 1.8]
ZONE_LABELS = {1.4: 'C | B', 1.8: 'B | A'}
Z_GRID = np.linspace(1.0, 2.2, 500)

# ---------------------------------------------------------------------------
# Read live constants from model and compute grand means
# ---------------------------------------------------------------------------
# 2026-09-20: 4 species (Sonneratia). Every grand mean below ("cross-species mean", the equalised
# level of the knockout / cascade / equal-CV designs): target = mean of the calibrated species;
# Sonneratia (uncalibrated, set by rule) blends toward the same target but does not define it. The
# target is APPLIED to all four species, so 'all species equal' stays literal.
DRAFT_VALS = dict(tdm.DRAFT)
MOBILITY_VALS = dict(tdm.MOBILITY)
SIGMA_LAT_VALS = dict(tdm.SIGMA_LATERAL)


def calibrated_mean(vals):
    """Equalisation target of a species-keyed trait: the mean over CALIBRATED_SPECIES (2026-09-20)."""
    return float(np.mean([vals[sp] for sp in CALIBRATED_SPECIES]))


DRAFT_MEAN = calibrated_mean(DRAFT_VALS)
MOBILITY_MEAN = calibrated_mean(MOBILITY_VALS)
SIGMA_LAT_MEAN = calibrated_mean(SIGMA_LAT_VALS)


def _full_species_dict(base, values):
    """A complete replacement dict for a species-keyed tdm attribute: a copy of the live dict
    (every key the model expects) with the species in `values` overwritten. Every override dict
    built in this module goes through here, so no species key can be dropped (2026-09-20)."""
    out = dict(base)
    out.update(values)
    return out


# ---------------------------------------------------------------------------
# Treatment definitions
# ---------------------------------------------------------------------------
# Each treatment specifies full replacement dicts for the tdm attribute(s)
# being equalized. Equalization sets every modelled species (2026-09-20: all
# four, Sonneratia included) to the target = mean of the calibrated species.
_DRAFT_EQ    = _full_species_dict(DRAFT_VALS,     {sp: DRAFT_MEAN     for sp in SPECIES_ORDER})
_MOBILITY_EQ = _full_species_dict(MOBILITY_VALS,  {sp: MOBILITY_MEAN  for sp in SPECIES_ORDER})
_SPREAD_EQ   = _full_species_dict(SIGMA_LAT_VALS, {sp: SIGMA_LAT_MEAN for sp in SPECIES_ORDER})

TREATMENTS = {
    'baseline':         {},
    'draft_eq_mean':    {'DRAFT':         _DRAFT_EQ},
    'mobility_eq_mean': {'MOBILITY':      _MOBILITY_EQ},
    'spread_eq_mean':   {'SIGMA_LATERAL': _SPREAD_EQ},
    'all_transport_eq': {'DRAFT':         _DRAFT_EQ,
                         'MOBILITY':      _MOBILITY_EQ,
                         'SIGMA_LATERAL': _SPREAD_EQ},
}


# ---------------------------------------------------------------------------
# Result tables
# ---------------------------------------------------------------------------

def _cache_path(results_dir, treatment, seed):
    return os.path.join(results_dir,
                        f'knockout_{treatment}_seed{seed}.csv')


def _save_results_csv(results, path):
    if not results:
        return
    fieldnames = ['species', 'x', 'y', 'zone', 'source', 'day', 'outcome',
                  'cause', 'survival_days', 'burial_mm', 'settle_method',
                  'first_grounding_t', 'first_grounding_x',
                  'first_grounding_y', 'first_grounding_z', 'terminal_z',
                  'tidal_state_at_settle']
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,
                                extrasaction='ignore')
        writer.writeheader()
        for r in results:
            writer.writerow(r)


_NUMERIC_FIELDS = {'x', 'y', 'day', 'survival_days', 'burial_mm',
                    'first_grounding_t', 'first_grounding_x',
                    'first_grounding_y', 'first_grounding_z', 'terminal_z'}


def _load_results_csv(path):
    results = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for k in _NUMERIC_FIELDS:
                v = row.get(k, '')
                if v != '':
                    try:
                        row[k] = float(v)
                    except (ValueError, TypeError):
                        pass
            results.append(row)
    return results


def _run_single(treatment, overrides, seed, n_days, cache_path,
                seed_source_mode=None):
    """Worker function for parallel execution (runs in subprocess)."""
    results = run_with_overrides(overrides, n_days=n_days, seed=seed,
                                 scenario='a',
                                 seed_source_mode=seed_source_mode)
    _save_results_csv(results, cache_path)
    return treatment, seed, len(results), cache_path


def run_treatments(treatments, seeds, n_days, results_dir, skip_runs=False,
                   seed_source_mode=None, max_workers=1):
    """Run all treatment x seed combinations with caching and parallelism.

    Returns {(treatment, seed): results_list}
    """
    all_results = {}
    pending = []

    # Load cached, collect pending
    for treatment in treatments:
        overrides = TREATMENTS[treatment]
        for seed in seeds:
            cache = _cache_path(results_dir, treatment, seed)
            if os.path.exists(cache):
                print(f"  Loading cached: {treatment} seed={seed}")
                all_results[(treatment, seed)] = _load_results_csv(cache)
            elif not skip_runs:
                pending.append((treatment, overrides, seed, cache))
            else:
                print(f"  SKIP (no cache): {treatment} seed={seed}")

    if not pending:
        return all_results

    print(f"  {len(pending)} runs to execute with {max_workers} workers...")

    if max_workers <= 1:
        for treatment, overrides, seed, cache in pending:
            print(f"  Running: {treatment} seed={seed} ({n_days} days)...")
            results = run_with_overrides(
                overrides, n_days=n_days, seed=seed, scenario='a',
                seed_source_mode=seed_source_mode,
            )
            _save_results_csv(results, cache)
            all_results[(treatment, seed)] = results
            print(f"    -> {len(results)} propagules")
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for treatment, overrides, seed, cache in pending:
                fut = executor.submit(
                    _run_single, treatment, overrides, seed, n_days, cache,
                    seed_source_mode)
                futures[fut] = (treatment, seed)
                print(f"  Submitted: {treatment} seed={seed}")

            for fut in as_completed(futures):
                treatment, seed = futures[fut]
                try:
                    _, _, n_prop, cache = fut.result()
                    print(f"  Done: {treatment} seed={seed} -> {n_prop} propagules")
                    all_results[(treatment, seed)] = _load_results_csv(cache)
                except Exception as e:
                    print(f"  FAILED: {treatment} seed={seed}: {e}")

    return all_results


# ---------------------------------------------------------------------------
# Runs with a summary file per run (calibrated model replicates, drag profiles)
# ---------------------------------------------------------------------------

def cache_path(results_dir, tag, seed, prefix='calib'):
    return os.path.join(results_dir, f'{prefix}_{tag}_seed{seed}.csv')


def _worker(tag, overrides, seed, n_days, cache, stand_csv):
    """Run one variant, cache the results, write its summary (runs in a subprocess)."""
    import cascade_metrics as cm
    t0 = time.time()
    results = run_with_overrides(overrides, n_days=n_days, seed=seed, scenario='a')
    _save_results_csv(results, cache)
    meta = {'tag': tag, 'seed': seed, 'n_days': n_days, 'overrides': overrides,
            'trees': stand_csv, 'minutes': round((time.time() - t0) / 60, 1)}
    summ = cm.summarise_run(results, stand_csv, meta=meta)
    with open(cm.summary_path_for(cache), 'w') as f:
        json.dump(summ, f, indent=1)
    return tag, seed, cm.compact_line(summ)


def _summarise_worker(tag, overrides, seed, n_days, cache, stand_csv):
    import cascade_metrics as cm
    meta = {'tag': tag, 'seed': seed, 'n_days': n_days, 'overrides': overrides, 'trees': stand_csv}
    summ = cm.summarise_cache(cache, stand_csv, force=True, meta=meta)
    return tag, seed, cm.compact_line(summ)


def run_grid(grid, seeds, n_days, results_dir, stand_csv, workers=1, skip_runs=False, prefix='calib',
             stand_for_tag=None):
    """Run or summarise every (tag, overrides) x seed; stand_for_tag(tag) may give a
    variant-specific stand CSV for the summary (default: stand_csv)."""
    import cascade_metrics as cm
    os.makedirs(results_dir, exist_ok=True)
    jobs = []
    for tag, ov in grid:
        for seed in seeds:
            cache = cache_path(results_dir, tag, seed, prefix)
            summ = cm.summary_path_for(cache)
            stand = stand_for_tag(tag) if stand_for_tag else stand_csv
            if os.path.exists(summ):
                with open(summ) as f:
                    print(f'  have    {tag} seed={seed} | {cm.compact_line(json.load(f))}')
            elif os.path.exists(cache):
                jobs.append((_summarise_worker, (tag, ov, seed, n_days, cache, stand)))
            elif not skip_runs:
                jobs.append((_worker, (tag, ov, seed, n_days, cache, stand)))
            else:
                print(f'  SKIP    {tag} seed={seed} (no cache)')
    if not jobs:
        return
    print(f'  {len(jobs)} jobs on {workers} workers')
    log = open(os.path.join(results_dir, f'{prefix}_runs.log'), 'a')
    if workers <= 1:
        for fn, args in jobs:
            tag, seed, line = fn(*args)
            msg = f'  done    {tag} seed={seed} | {line}'
            print(msg, flush=True)
            log.write(msg + '\n')
            log.flush()
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(fn, *args): args[0] for fn, args in jobs}
            for fut in as_completed(futs):
                try:
                    tag, seed, line = fut.result()
                    msg = f'  done    {tag} seed={seed} | {line}'
                except Exception as e:  # keep the pass going
                    msg = f'  FAILED  {futs[fut]}: {e!r}'
                print(msg, flush=True)
                log.write(msg + '\n')
                log.flush()
    log.close()


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------
