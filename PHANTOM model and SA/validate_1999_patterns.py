"""
Pattern-Oriented Validation Against 1999 Field Data
====================================================

Zone-level validation pipeline (6 criteria) with sole project import
from tidal_dispersal_model.

Usage:
    python validate_1999_patterns.py                        # defaults
    python validate_1999_patterns.py --days 365 --seed 42   # full run
    python validate_1999_patterns.py --days 30 --seed 42    # smoke test
    python validate_1999_patterns.py --full                  # + replicates, density ratios
"""

import os
import csv
import json
import math
import time
import shutil
import argparse
from collections import defaultdict, Counter

from tidal_dispersal_model import run_scenario, zone_at, X_MIN, X_MAX, Y_MIN, Y_MAX
from tidal_dispersal_model import SPECIES as _MODEL_SPECIES, CALIBRATED_SPECIES as _MODEL_CALIBRATED
from transect_geometry import assign_to_transect, load_observed_quadrats_per_zone

try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

# ===================================================================
# Constants
# ===================================================================

# 2026-09-20: 4 species (Sonneratia). SPECIES (all modelled species) is used for the density
# tables and the mortality statistics; the acceptance criteria Z1, Z2, Z4, Z5 are evaluated
# on CALIBRATED_SPECIES (the three species of the 1999 calibration target), so the '/3'
# denominators and the thresholds keep their meaning. Sonneratia's constants are set by rule.
SPECIES = list(_MODEL_SPECIES)
CALIBRATED_SPECIES = list(_MODEL_CALIBRATED)
N_CAL = len(CALIBRATED_SPECIES)
ZONES = ['A', 'B', 'C']

# 2026-09-20: the juvenile quadrat area is NOT documented (data audit); the model-side
# densities use the ASSUMED virtual-sampling footprint of transect_geometry (25 m2, the
# paper's 5 x 5 m quadrat). Zone shares are footprint-free; absolute densities are conditional.
from transect_geometry import MODEL_QUADRAT_FOOTPRINT_M2  # noqa: E402
QUADRAT_AREA = MODEL_QUADRAT_FOOTPRINT_M2
# Surveyed quadrats per position zone (1999), read from the re-zoned quadrat records
# (was the vegetation-class count {'A': 83, 'B': 170, 'C': 82} before 2026-09-04)
N_QUADRATS_PER_ZONE = load_observed_quadrats_per_zone()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OBS_CSV_1999 = os.path.join(SCRIPT_DIR, 'juvenile_zone_species_summary_1999_t1_5.csv')

# --- Acceptance thresholds ---
THRESH_SPEARMAN = 0.6
THRESH_PEAK_ZONE_MATCHES = 2
THRESH_MORT_C_LO, THRESH_MORT_C_HI = 40.0, 70.0
THRESH_MORT_A_LO, THRESH_MORT_A_HI = 10.0, 30.0
THRESH_MEAN_JSD = 0.2
THRESH_JSD_ZONES_PASS = 2


# ===================================================================
# Math Helpers
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


def model_zone_sampled_densities(results):
    """Compute zone x species densities using transect-filtered, quadrat-matched area."""
    established = [r for r in results if r['outcome'] == 'established']
    counts = Counter()
    for r in established:
        t = assign_to_transect(r['x'], r['y'])
        if t is None:
            continue  # skip propagules outside transect polygons
        counts[(r['zone'], r['species'])] += 1

    densities = {}
    for z in ZONES:
        sampled_area = N_QUADRATS_PER_ZONE[z] * QUADRAT_AREA
        for sp in SPECIES:
            key = (z, sp)
            densities[key] = counts[key] / sampled_area if sampled_area > 0 else 0.0
    return densities


# ===================================================================
# Observed Data Loaders
# ===================================================================

def load_observed_zone_densities(csv_path=None):
    """Load 1999 observed densities: {(zone, species): mean_density}."""
    if csv_path is None:
        csv_path = OBS_CSV_1999
    obs = {}
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['species'] not in SPECIES:   # Bruguiera (not modelled) is skipped, 2026-09-20
                continue
            key = (row['zone'], row['species'])
            obs[key] = float(row['mean_density'])
    return obs


# ===================================================================
# Mortality Helpers
# ===================================================================

def mortality_stats(results):
    """Seedling mortality % per zone among rooted propagules.

    Propagules that failed to root on unsuitable substrate (cause 'substrate',
    Eq. 10) never became seedlings and are excluded from both numerator and
    denominator, so the statistic stays comparable with the 40-70 % / 10-30 %
    seedling mortality of Balke et al. (2011).
    """
    stats = {}
    for z in ZONES:
        rooted = [r for r in results if r['zone'] == z and r['outcome'] in ('established', 'dead')
                  and r.get('cause') != 'substrate']
        dead = sum(1 for r in rooted if r['outcome'] == 'dead')
        stats[z] = 100.0 * dead / len(rooted) if rooted else 0.0
    return stats


def mortality_cause_analysis(results):
    """Per-zone mortality cause shares."""
    analysis = {}
    for z in ZONES:
        dead = [r for r in results if r['zone'] == z and r['outcome'] == 'dead'
                and r.get('cause') != 'substrate']   # rooting failures are not seedling deaths
        lost = [r for r in results if r['zone'] == z and r['outcome'] == 'lost']
        all_dead = dead + lost
        total = len(all_dead)
        if total == 0:
            analysis[z] = {'siltation': 0.0, 'predation': 0.0, 'washed_out': 0.0,
                           'total_dead': 0}
            continue
        causes = defaultdict(int)
        for r in all_dead:
            causes[r.get('cause', 'unknown')] += 1
        analysis[z] = {
            'siltation': causes.get('siltation', 0) / total,
            'predation': causes.get('predation', 0) / total,
            'washed_out': causes.get('washed_out', 0) / total,
            'total_dead': total,
        }
    return analysis


# ===================================================================
# Acceptance Metrics (zone-level)
# ===================================================================

def zone_spearman_rho(model_dens, obs_dens, species=None):
    """9-cell zone x species Spearman rank correlation (calibrated species by default)."""
    species = CALIBRATED_SPECIES if species is None else species
    cells = [(z, sp) for z in ZONES for sp in species]
    return _spearman([model_dens.get(k, 0.0) for k in cells],
                     [obs_dens.get(k, 0.0) for k in cells])


def species_peak_zone_match(model_dens, obs_dens, species=None):
    """Check each species' peak zone against observed (calibrated species by default)."""
    species = CALIBRATED_SPECIES if species is None else species
    zone_rank = {z: i for i, z in enumerate(ZONES)}
    details = {}
    rank_distances = {}
    matches = 0
    near_misses = 0
    for sp in species:
        m_vals = {z: model_dens.get((z, sp), 0.0) for z in ZONES}
        o_vals = {z: obs_dens.get((z, sp), 0.0) for z in ZONES}
        m_peak = max(ZONES, key=lambda z: m_vals[z])
        o_peak = max(ZONES, key=lambda z: o_vals[z])
        match = m_peak == o_peak
        rd = abs(zone_rank[m_peak] - zone_rank[o_peak])
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


def zone_mortality_gradient(results):
    """Check mortality gradient C > B > A with range targets."""
    mort = mortality_stats(results)
    ordering_correct = mort['C'] > mort['B'] > mort['A']
    a_in_range = THRESH_MORT_A_LO <= mort['A'] <= THRESH_MORT_A_HI
    c_in_range = THRESH_MORT_C_LO <= mort['C'] <= THRESH_MORT_C_HI
    passed = ordering_correct and a_in_range and c_in_range
    return {
        'mortality': mort,
        'ordering_correct': ordering_correct,
        'a_in_range': a_in_range,
        'c_in_range': c_in_range,
        'pass': passed,
    }


def within_zone_species_rank(model_dens, obs_dens, species=None):
    """Per-zone species rank correlation (Spearman on the calibrated species by default)."""
    species = CALIBRATED_SPECIES if species is None else species
    per_zone = {}
    n_agreeing = 0
    for z in ZONES:
        m = [model_dens.get((z, sp), 0.0) for sp in species]
        o = [obs_dens.get((z, sp), 0.0) for sp in species]
        mr = _rank(m)
        orank = _rank(o)
        n = len(m)
        m_mean = sum(mr) / n
        o_mean = sum(orank) / n
        cov = sum((a - m_mean) * (b - o_mean) for a, b in zip(mr, orank))
        var_m = sum((a - m_mean) ** 2 for a in mr)
        var_o = sum((b - o_mean) ** 2 for b in orank)
        denom = math.sqrt(var_m * var_o)
        rho = cov / denom if denom > 0 else 0.0
        agree = rho >= 0.5
        if agree:
            n_agreeing += 1
        per_zone[z] = {'spearman': rho, 'agree': agree}
    return {
        'per_zone': per_zone,
        'n_agreeing_zones': n_agreeing,
        'pass': n_agreeing >= 2,
    }


def zone_composition_jsd(model_dens, obs_dens, species=None):
    """Species-proportion similarity per zone via JSD (calibrated species by default)."""
    species = CALIBRATED_SPECIES if species is None else species
    per_zone = {}
    for z in ZONES:
        m = [model_dens.get((z, sp), 0.0) for sp in species]
        o = [obs_dens.get((z, sp), 0.0) for sp in species]
        jsd = _jsd(m, o)
        per_zone[z] = {'jsd': jsd}
    jsd_values = [per_zone[z]['jsd'] for z in ZONES]
    mean_jsd = sum(jsd_values) / len(jsd_values) if jsd_values else 0.0
    n_passing = sum(1 for v in jsd_values if v < THRESH_MEAN_JSD)
    return {
        'per_zone': per_zone,
        'mean_jsd': mean_jsd,
        'n_zones_pass': n_passing,
        'pass': n_passing >= THRESH_JSD_ZONES_PASS,
    }


# ===================================================================
# Reviewer Diagnostics (--full mode)
# ===================================================================

def replicate_stability(n_reps, n_days, csv_dir, obs_dens_zone):
    """Run model with multiple seeds, return per-seed metric summary.

    Returns dict: {zone_rho: {mean, sd, values}, ...}
    Returns None if n_reps == 0.
    """
    if n_reps == 0:
        return None
    seeds = [42 + i * 7 for i in range(n_reps)]
    zone_rhos = []
    zone_jsds = []

    for si, seed in enumerate(seeds):
        print(f"    Replicate {si + 1}/{n_reps}: seed={seed}", end='', flush=True)
        t0 = time.time()
        results = run_scenario('a', n_days=n_days, seed=seed, csv_dir=csv_dir)
        elapsed = time.time() - t0

        m_zone = model_zone_sampled_densities(results)

        zr = zone_spearman_rho(m_zone, obs_dens_zone)
        zj = zone_composition_jsd(m_zone, obs_dens_zone)['mean_jsd']

        zone_rhos.append(zr)
        zone_jsds.append(zj)
        print(f"  z_rho={zr:.3f} ({elapsed:.0f}s)")

    def _stats(vals):
        n = len(vals)
        mean = sum(vals) / n
        sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / max(1, n - 1)) if n > 1 else 0.0
        return {'mean': round(mean, 4), 'sd': round(sd, 4), 'values': [round(v, 4) for v in vals]}

    return {
        'zone_spearman_rho': _stats(zone_rhos),
        'zone_mean_jsd': _stats(zone_jsds),
        'seeds': seeds,
    }


def density_ratio_table(m_dens_zone, obs_dens_zone):
    """Compute model/observed density ratios per zone.

    Returns dict: {zone: {model_total, obs_total, ratio}}
    """
    table = {}
    for z in ZONES:
        m_total = sum(m_dens_zone.get((z, sp), 0.0) for sp in SPECIES)
        o_total = sum(obs_dens_zone.get((z, sp), 0.0) for sp in SPECIES)
        table[z] = {
            'model_total': round(m_total, 4),
            'obs_total': round(o_total, 4),
            'ratio': round(m_total / o_total, 2) if o_total > 0 else None,
        }
    return table


# ===================================================================
# Acceptance Evaluator
# ===================================================================

def evaluate_acceptance(results, m_dens_zone, obs_dens_zone):
    """Evaluate all 6 zone-level criteria. Returns (criteria_list, summary_dict)."""
    criteria = []

    rho = zone_spearman_rho(m_dens_zone, obs_dens_zone)
    criteria.append({
        'criterion': '9-cell Spearman rho',
        'value': f"{rho:.3f}",
        'threshold': f">= {THRESH_SPEARMAN}",
        'pass': rho >= THRESH_SPEARMAN,
        'type': 'guardrail',
        'note': f'Rank correlation across 3 zones x {N_CAL} calibrated species',
    })

    peak_z = species_peak_zone_match(m_dens_zone, obs_dens_zone)
    criteria.append({
        'criterion': 'Peak zone matches',
        'value': f"{peak_z['matches']}/{N_CAL}",
        'threshold': f">= {THRESH_PEAK_ZONE_MATCHES}/{N_CAL}",
        'pass': peak_z['matches'] >= THRESH_PEAK_ZONE_MATCHES,
        'type': 'guardrail',
        'note': 'Species with correct highest-density zone',
    })

    mort_grad = zone_mortality_gradient(results)
    mort = mort_grad['mortality']
    criteria.append({
        'criterion': 'Mortality gradient C>B>A',
        'value': f"A={mort['A']:.1f}%, C={mort['C']:.1f}%",
        'threshold': (f"A in [{THRESH_MORT_A_LO},{THRESH_MORT_A_HI}], "
                      f"C in [{THRESH_MORT_C_LO},{THRESH_MORT_C_HI}], C>B>A"),
        'pass': mort_grad['pass'],
        'type': 'guardrail',
        'note': 'Balke et al. (2011) targets',
    })

    wz = within_zone_species_rank(m_dens_zone, obs_dens_zone)
    criteria.append({
        'criterion': 'Within-zone species rank',
        'value': f"{wz['n_agreeing_zones']}/3 zones agree",
        'threshold': '>= 2/3 zones',
        'pass': wz['pass'],
        'type': 'guardrail',
        'note': f'Per-zone {N_CAL}-species Spearman >= 0.5 (calibrated species)',
    })

    comp_jsd = zone_composition_jsd(m_dens_zone, obs_dens_zone)
    criteria.append({
        'criterion': 'Composition JSD',
        'value': f"mean={comp_jsd['mean_jsd']:.3f}",
        'threshold': f"< {THRESH_MEAN_JSD} for >= {THRESH_JSD_ZONES_PASS}/3 zones",
        'pass': comp_jsd['pass'],
        'type': 'guardrail',
        'note': 'Jensen-Shannon divergence of species proportions',
    })

    mort_causes = mortality_cause_analysis(results)
    c_causes = mort_causes.get('C', {})
    pred_dominant = c_causes.get('predation', 0.0) > c_causes.get('siltation', 0.0)
    criteria.append({
        'criterion': 'Predation dominant in Zone C',
        'value': (f"pred={c_causes.get('predation', 0):.1%}, "
                  f"silt={c_causes.get('siltation', 0):.1%}"),
        'threshold': 'predation > siltation in Zone C',
        'pass': pred_dominant,
        'type': 'guardrail',
        'note': 'Clarke & Kerrigan (2002)',
    })

    # Summary
    n_total = len(criteria)
    n_pass = sum(1 for c in criteria if c['pass'])

    summary = {
        'passed': n_pass,
        'total': n_total,
        'criteria_species': list(CALIBRATED_SPECIES),   # 2026-09-20: Z1/Z2/Z4/Z5 on the calibrated species
        'reported_species': list(SPECIES),               # density tables and mortality include Sonneratia
    }

    return criteria, summary


# ===================================================================
# Output Layer
# ===================================================================

def _write_aggregated_csv(m_dens_zone, obs_dens_zone, out_path):
    """Write aggregated zone×species CSV with obs, model, bias, ratio columns.

    Format matches validation_aggregated_1999.csv layout.
    """
    rows = []

    # Domain totals
    for sp in SPECIES:
        obs_total = sum(obs_dens_zone.get((z, sp), 0.0) for z in ZONES)
        mod_total = sum(m_dens_zone.get((z, sp), 0.0) for z in ZONES)
        bias = mod_total - obs_total
        ratio = mod_total / obs_total if obs_total > 0 else None
        rows.append(('domain_total', sp, obs_total, mod_total, bias, ratio))

    # Zone × species
    for z in ZONES:
        for sp in SPECIES:
            obs_val = obs_dens_zone.get((z, sp), 0.0)
            mod_val = m_dens_zone.get((z, sp), 0.0)
            bias = mod_val - obs_val
            ratio = mod_val / obs_val if obs_val > 0 else None
            rows.append(('zone_species', f'{z}_{sp}', obs_val, mod_val, bias, ratio))

    # Summary GOF stats
    zs_obs = [obs_dens_zone.get((z, sp), 0.0) for z in ZONES for sp in SPECIES]
    zs_mod = [m_dens_zone.get((z, sp), 0.0) for z in ZONES for sp in SPECIES]
    n = len(zs_obs)
    sq_err = sum((m - o) ** 2 for m, o in zip(zs_mod, zs_obs))
    rmse = math.sqrt(sq_err / n) if n > 0 else 0.0
    mae = sum(abs(m - o) for m, o in zip(zs_mod, zs_obs)) / n if n > 0 else 0.0
    mean_obs = sum(zs_obs) / n if n > 0 else 0.0
    ss_tot = sum((o - mean_obs) ** 2 for o in zs_obs)
    r2 = 1 - sq_err / ss_tot if ss_tot > 0 else 0.0

    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['scale', 'cell', 'obs', 'model', 'bias', 'ratio'])
        for row in rows:
            scale, cell, obs, mod, bias, ratio = row
            ratio_str = f'{ratio:.4f}' if ratio is not None else ''
            writer.writerow([scale, cell, f'{obs:.6f}', f'{mod:.6f}',
                             f'{bias:.6f}', ratio_str])
        writer.writerow(['summary', 'zone_species_gof', '', '',
                         f'RMSE={rmse:.6f} MAE={mae:.6f} R2={r2:.4f}', ''])
    print(f"  Written: {out_path}")


def _update_latest_alias(src, dst, new_days, force=False):
    """Copy src to dst (latest alias), guarded by run length.

    If dst already exists and was produced by a longer run, skip the
    overwrite unless --force is set.  For JSON files the guard reads
    parameters.days directly; for CSV/PNG it looks for the sibling JSON.
    """
    if os.path.exists(dst) and not force:
        existing_days = None
        # Try reading days from the destination itself (JSON) or sibling JSON
        json_to_check = dst if dst.endswith('.json') else None
        if json_to_check is None:
            # Look for sibling JSON in same dir, then parent dir (for figures/)
            for search_dir in [os.path.dirname(dst),
                               os.path.dirname(os.path.dirname(dst))]:
                sibling = os.path.join(search_dir,
                                       'validation_1999_patterns.json')
                if os.path.exists(sibling):
                    json_to_check = sibling
                    break
        if json_to_check and os.path.exists(json_to_check):
            try:
                with open(json_to_check) as f:
                    existing_days = json.load(f).get('parameters', {}).get('days')
            except Exception:
                pass
        if existing_days is not None and existing_days > new_days:
            print(f"  Skipped alias {os.path.basename(dst)}: "
                  f"existing {existing_days}d > new {new_days}d (use --force to override)")
            return
    shutil.copy2(src, dst)
    print(f"  Updated alias: {os.path.basename(dst)}")


def write_results_json(criteria, summary, provenance, out_path,
                       replicates=None, density_ratios=None):
    """Write compact validation JSON."""
    params = provenance.get('parameters', {})
    run_tag = f"{params.get('days', '?')}d_seed{params.get('seed', '?')}"
    payload = {
        'validation': '1999_pattern_agreement',
        'version': '4.0',
        'run_tag': run_tag,
        'parameters': params,
        'summary': summary,
        'acceptance_criteria': criteria,
        'provenance': {
            'n_propagules': provenance.get('n_propagules', 0),
            'n_established': provenance.get('n_established', 0),
            'elapsed_s': provenance.get('elapsed_s', 0),
            'obs_zone_csv': os.path.basename(OBS_CSV_1999),
        },
    }
    if replicates is not None:
        payload['replicate_stability'] = replicates
    if density_ratios is not None:
        payload['density_ratios'] = density_ratios
    with open(out_path, 'w') as f:
        json.dump(payload, f, indent=2, default=str)
    print(f"  Written: {out_path}")


def print_report(criteria, summary, replicates=None, density_ratios=None):
    """Compact console report with optional diagnostics."""
    print("\n" + "=" * 65)
    print("  ACCEPTANCE SUMMARY")
    print("=" * 65)

    print(f"  {'Criterion':<35s} {'Value':<18s} {'Result':>6s}")
    print("  " + "-" * 62)
    for c in criteria:
        status = 'PASS' if c['pass'] else 'FAIL'
        print(f"  {c['criterion']:<35s} {c['value']:<18s} {status:>6s}")
    print("  " + "-" * 62)
    print(f"  Criteria: {summary['passed']}/{summary['total']} passed")

    if replicates:
        zr = replicates['zone_spearman_rho']
        print(f"\n  Replicate stability ({len(replicates['seeds'])} seeds):")
        print(f"    Zone rho: {zr['mean']:.3f} +/- {zr['sd']:.3f}  {zr['values']}")

    if density_ratios:
        print(f"\n  Density ratios (model/observed):")
        for z in ZONES:
            d = density_ratios[z]
            r_str = f"{d['ratio']:.2f}x" if d['ratio'] is not None else "N/A"
            print(f"    Zone {z}: model={d['model_total']:.4f}  obs={d['obs_total']:.4f}  ratio={r_str}")


# ===================================================================
# Validation Figure (2x2 panels)
# ===================================================================

SPECIES_COLORS = {'Avicennia': '#2171B5', 'Rhizophora': '#D95F02', 'Ceriops': '#7570B3',
                  'Sonneratia': '#1B9E77'}   # green, 2026-09-20 (palette is Blues/Oranges/Purples)
ZONE_MARKERS = {'A': 'o', 'B': 's', 'C': '^'}
AMBER_COLOR = '#F39C12'


def _get_criterion_value(data, name):
    """Extract value string for a named criterion."""
    for c in data.get('acceptance_criteria', []):
        if c['criterion'] == name:
            return c['value']
    return '?'


def build_agg_from_dicts(m_dens_zone, obs_dens_zone):
    """Build agg list from in-memory density dicts (no CSV round-trip)."""
    rows = []
    for z in ZONES:
        for sp in SPECIES:
            rows.append({
                'zone': z, 'species': sp,
                'obs': obs_dens_zone.get((z, sp), 0.0),
                'model': m_dens_zone.get((z, sp), 0.0),
            })
    return rows


def _panel_rank_parity(ax, data, agg):
    """Panel (a): Zone x Species rank agreement scatter."""
    obs_vals = [r['obs'] for r in agg]
    mod_vals = [r['model'] for r in agg]
    obs_rank = _rank(obs_vals)
    mod_rank = _rank(mod_vals)
    for i, r in enumerate(agg):
        ax.scatter(obs_rank[i], mod_rank[i], c=SPECIES_COLORS[r['species']],
                   marker=ZONE_MARKERS[r['zone']], s=100, edgecolors='black',
                   linewidths=0.5, zorder=3)
    ax.plot([0, 10], [0, 10], 'k--', alpha=0.3, linewidth=0.8)
    from scipy.stats import spearmanr
    rho_val = _get_criterion_value(data, '9-cell Spearman rho')
    _, p_val = spearmanr(obs_vals, mod_vals)
    rep = data.get('replicate_stability', {}).get('zone_spearman_rho', {})
    ci_text = ''
    if rep:
        ci_text = f"\n95% CI [{rep['mean'] - 1.96*rep['sd']:.3f}, {min(rep['mean'] + 1.96*rep['sd'], 1.0):.3f}]"
    p_text = f', p = {p_val:.4f}' if p_val >= 0.001 else ', p < 0.001'
    ax.set_xlabel('Observed rank', fontsize=10)
    ax.set_ylabel('Model rank', fontsize=10)
    ax.set_title('(a) Zone\u00d7Species Rank Agreement', fontsize=13, fontweight='bold')
    ax.text(0.05, 0.92, f'\u03c1 = {rho_val}{p_text} (n = {len(obs_vals)}){ci_text}',
            transform=ax.transAxes, fontsize=9, verticalalignment='top')
    sp_handles = [plt.Line2D([0], [0], marker='o', color='w',
                             markerfacecolor=SPECIES_COLORS[s], markersize=8,
                             markeredgecolor='black', markeredgewidth=0.5, label=s)
                  for s in SPECIES]
    zn_handles = [plt.Line2D([0], [0], marker=ZONE_MARKERS[z], color='w',
                             markerfacecolor='gray', markersize=8,
                             markeredgecolor='black', markeredgewidth=0.5, label=f'Zone {z}')
                  for z in ZONES]
    ax.legend(handles=sp_handles + zn_handles, fontsize=7, loc='lower right', ncol=2)


def _panel_log_density(ax, data, agg):
    """Panel (b): Log-scale absolute density scatter."""
    for r in agg:
        ax.scatter(r['obs'], r['model'], c=SPECIES_COLORS[r['species']],
                   marker=ZONE_MARKERS[r['zone']], s=100, edgecolors='black',
                   linewidths=0.5, zorder=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    lims = [0.008, 5]
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.plot(lims, lims, 'k--', alpha=0.3, linewidth=0.8, label='1:1')
    ratios = [r['model'] / r['obs'] for r in agg if r['obs'] > 0]
    median_ratio = np.median(ratios) if ratios else 1.0
    ax.plot(lims, [lims[0] * median_ratio, lims[1] * median_ratio], '--',
            color=AMBER_COLOR, alpha=0.6, linewidth=1, label=f'Median ratio = {median_ratio:.2f}\u00d7')
    ax.set_xlabel('Observed density (m\u207b\u00b2)', fontsize=10)
    ax.set_ylabel('Model density (m\u207b\u00b2)', fontsize=10)
    ax.set_title('(b) Absolute Densities (Log Scale)', fontsize=13, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    sp_handles = [plt.Line2D([0], [0], marker='o', color='w',
                             markerfacecolor=SPECIES_COLORS[s], markersize=8,
                             markeredgecolor='black', markeredgewidth=0.5, label=s)
                  for s in SPECIES]
    zn_handles = [plt.Line2D([0], [0], marker=ZONE_MARKERS[z], color='w',
                             markerfacecolor='gray', markersize=8,
                             markeredgecolor='black', markeredgewidth=0.5, label=f'Zone {z}')
                  for z in ZONES]
    ax2_leg = ax.legend(handles=sp_handles + zn_handles, fontsize=7, loc='lower right', ncol=2)
    ax.add_artist(ax2_leg)
    ax.legend(fontsize=7, loc='upper left')


def _panel_normalized_props(ax, data, agg):
    """Panel (c): Normalized proportions scatter."""
    for zone in ZONES:
        zone_rows = [r for r in agg if r['zone'] == zone]
        obs_total = sum(r['obs'] for r in zone_rows)
        mod_total = sum(r['model'] for r in zone_rows)
        if obs_total == 0 or mod_total == 0:
            continue
        for r in zone_rows:
            obs_prop = r['obs'] / obs_total
            mod_prop = r['model'] / mod_total
            ax.scatter(obs_prop, mod_prop, c=SPECIES_COLORS[r['species']],
                       marker=ZONE_MARKERS[r['zone']], s=100, edgecolors='black',
                       linewidths=0.5, zorder=3)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=0.8)
    ax.set_xlabel('Observed proportion', fontsize=10)
    ax.set_ylabel('Model proportion', fontsize=10)
    ax.set_title('(c) Normalized Proportions', fontsize=13, fontweight='bold')
    obs_props_all, mod_props_all = [], []
    for zone in ZONES:
        zone_rows = [r for r in agg if r['zone'] == zone]
        ot = sum(r['obs'] for r in zone_rows)
        mt = sum(r['model'] for r in zone_rows)
        if ot > 0 and mt > 0:
            for r in zone_rows:
                obs_props_all.append(r['obs'] / ot)
                mod_props_all.append(r['model'] / mt)
    if obs_props_all:
        obs_arr = np.array(obs_props_all)
        mod_arr = np.array(mod_props_all)
        ss_res = np.sum((obs_arr - mod_arr) ** 2)
        ss_tot = np.sum((obs_arr - np.mean(obs_arr)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        ax.text(0.05, 0.92, f'R\u00b2 = {r2:.3f}', transform=ax.transAxes, fontsize=9,
                verticalalignment='top')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)


def _panel_zone_density_bars(ax, agg):
    """Panel (d): Zone-level grouped bars, observed vs model (absolute density)."""
    bar_w = 0.35
    tick_positions = []
    tick_labels = []
    group = len(SPECIES) + 1   # one slot per species plus a gap between zones (2026-09-20: 4 species)
    for iz, zone in enumerate(ZONES):
        base = iz * group
        for j, sp in enumerate(SPECIES):
            obs_val = mod_val = 0.0
            for r in agg:
                if r['zone'] == zone and r['species'] == sp:
                    obs_val = r['obs']
                    mod_val = r['model']
            x_obs = base + j - bar_w / 2
            x_mod = base + j + bar_w / 2
            ax.bar(x_obs, obs_val, bar_w, color=SPECIES_COLORS[sp],
                   edgecolor='black', linewidth=0.5, alpha=0.9)
            ax.bar(x_mod, mod_val, bar_w, color=SPECIES_COLORS[sp],
                   edgecolor='black', linewidth=0.5, alpha=0.5, hatch='///')
        tick_positions.append(base + (len(SPECIES) - 1) / 2)
        tick_labels.append(f'Zone {zone}')

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_ylabel('Density (seedlings m\u207b\u00b2)', fontsize=10)
    ax.set_title('(c) Zone-Level Density', fontsize=13, fontweight='bold')

    sp_patches = [Patch(facecolor=SPECIES_COLORS[s], edgecolor='black',
                        linewidth=0.5, label=s) for s in SPECIES]
    obs_patch = Patch(facecolor='white', edgecolor='black', linewidth=0.5,
                      alpha=0.9, label='Observed')
    mod_patch = Patch(facecolor='white', edgecolor='black', linewidth=0.5,
                      alpha=0.5, hatch='///', label='Model')
    ax.legend(handles=sp_patches + [obs_patch, mod_patch], fontsize=7,
              loc='upper right', ncol=2)


def fig_validation(data, agg, out_path):
    """Generate 1x3 validation figure: rank parity, log density, zone bars."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    _panel_rank_parity(axes[0], data, agg)
    _panel_log_density(axes[1], data, agg)
    _panel_zone_density_bars(axes[2], agg)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


# ===================================================================
# Main
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Pattern validation against 1999 field data')
    parser.add_argument('--days', '-d', type=int, default=365)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output-dir', type=str,
                        default=os.path.join(SCRIPT_DIR, 'results'))
    parser.add_argument('--full', action='store_true',
                        help='Include replicate stability and density ratios')
    parser.add_argument('--replicates', type=int, default=5,
                        help='Number of replicate seeds for --full (default: 5)')
    parser.add_argument('--force', action='store_true',
                        help='Force overwrite of latest alias even if existing run used more days')
    parser.add_argument('--trees', default=None,
                        help='tree inventory CSV for scenario (a) (default: the model scenario file, demo.csv)')
    args = parser.parse_args()
    if args.trees:
        import tidal_dispersal_model as _tdm
        _t = os.path.abspath(args.trees)
        for _k in ('a', 'b'):
            _n, _b, _i, _d, _ = _tdm.SCENARIOS[_k]
            _tdm.SCENARIOS[_k] = (_n, _b, _i, _d, _t)
        print(f'  Trees: {_t}')

    os.makedirs(args.output_dir, exist_ok=True)
    csv_dir = SCRIPT_DIR

    n_steps = 6 if args.full else 5

    print("=" * 65)
    print("1999 PATTERN VALIDATION")
    print("=" * 65)
    print(f"  days={args.days}  seed={args.seed}  full={args.full}")

    # [1] Load observed data
    print(f"\n[1/{n_steps}] Loading observed 1999 data...")
    if not os.path.exists(OBS_CSV_1999):
        print(f"  Error: {OBS_CSV_1999} not found.")
        return
    obs_dens_zone = load_observed_zone_densities()
    print(f"  Zone densities: {len(obs_dens_zone)} cells")

    # [2] Run model
    print(f"\n[2/{n_steps}] Running model (scenario 'a', {args.days} days, seed={args.seed})...")
    t0 = time.time()
    results = run_scenario('a', n_days=args.days, seed=args.seed, csv_dir=csv_dir)
    elapsed = time.time() - t0
    n_est = sum(1 for r in results if r['outcome'] == 'established')
    n_dead = sum(1 for r in results if r['outcome'] == 'dead')
    print(f"  {len(results)} propagules, {n_est} established, {n_dead} dead ({elapsed:.0f}s)")

    # [3] Compute model densities
    print(f"\n[3/{n_steps}] Computing model densities...")
    m_dens_zone = model_zone_sampled_densities(results)
    print(f"  Zone: {len(m_dens_zone)} cells")

    # [4] Evaluate 6 criteria
    print(f"\n[4/{n_steps}] Evaluating acceptance criteria...")
    criteria, summary = evaluate_acceptance(results, m_dens_zone, obs_dens_zone)

    # Optional reviewer diagnostics
    replicates = None
    density_ratios = None

    if args.full:
        # [5] Replicate stability
        print(f"\n[5/{n_steps}] Replicate stability ({args.replicates} seeds)...")
        replicates = replicate_stability(
            args.replicates, args.days, csv_dir, obs_dens_zone)

        # Density ratios (computed from primary run, no extra cost)
        density_ratios = density_ratio_table(m_dens_zone, obs_dens_zone)

    # Write outputs
    step_write = n_steps - 1 if not args.full else n_steps - 1
    print(f"\n[{step_write}/{n_steps}] Writing outputs...")
    provenance = {
        'parameters': {
            'days': args.days, 'seed': args.seed, 'scenario': 'a',
            'full': args.full, 'replicates': args.replicates if args.full else 0,
        },
        'n_propagules': len(results),
        'n_established': n_est,
        'elapsed_s': round(elapsed, 1),
    }

    # Build run tag for artifact isolation
    tag = f"{args.days}d_seed{args.seed}"

    # Tagged JSON (always written)
    json_tagged = os.path.join(args.output_dir, f'validation_1999_patterns_{tag}.json')
    write_results_json(criteria, summary, provenance, json_tagged,
                       replicates=replicates, density_ratios=density_ratios)

    # Latest JSON alias (guarded)
    json_latest = os.path.join(args.output_dir, 'validation_1999_patterns.json')
    _update_latest_alias(json_tagged, json_latest, args.days, args.force)

    # Tagged CSV (always written)
    csv_tagged = os.path.join(args.output_dir, f'validation_aggregated_1999_{tag}.csv')
    _write_aggregated_csv(m_dens_zone, obs_dens_zone, csv_tagged)

    # Latest CSV alias (guarded)
    csv_latest = os.path.join(args.output_dir, 'validation_aggregated_1999.csv')
    _update_latest_alias(csv_tagged, csv_latest, args.days, args.force)

    # Generate validation figure
    print(f"\n[{n_steps}/{n_steps}] Generating validation figure...")
    if _HAS_MPL:
        try:
            fig_dir = os.path.join(args.output_dir, 'figures')
            os.makedirs(fig_dir, exist_ok=True)
            agg = build_agg_from_dicts(m_dens_zone, obs_dens_zone)
            data = {
                'acceptance_criteria': criteria,
                'replicate_stability': replicates or {},
            }
            # Tagged figure
            fig_tagged = os.path.join(fig_dir, f'fig_validation_{tag}.png')
            fig_validation(data, agg, fig_tagged)
            # Latest figure alias (guarded)
            fig_latest = os.path.join(fig_dir, 'fig_validation.png')
            _update_latest_alias(fig_tagged, fig_latest, args.days, args.force)
        except Exception as e:
            print(f"  Warning: figure generation failed: {e}")
    else:
        print("  Skipped: matplotlib not available")

    print_report(criteria, summary, replicates=replicates,
                 density_ratios=density_ratios)
    print("\nDone.")


if __name__ == '__main__':
    main()
