"""Virtual quadrat survey of the model: sample the simulated recruits the way the field survey
sampled the forest, so that the simulated zone shares carry the same quadrat-sampling uncertainty
as the observed ones.

Observed side (2026-09-20, data audit docs/juvenile_data_audit.md): the survey records a DENSITY
(#/m2) per quadrat; the quadrat area is not documented and no verified counts exist. The observed
statistics are therefore computed from `density_per_m2` per quadrat, nothing else.

Model side: for every replicate run of the full model, the established recruits are counted inside
a circle of the ASSUMED footprint MODEL_QUADRAT_FOOTPRINT_M2 centred on each quadrat position of the
juvenile records (zone = the quadrat's position zone) and converted to a density,
count / footprint. The footprint (25 m2 = the paper's 5 x 5 m quadrat per PCQM quadrant; sensitivity
12.5 m2, the unit the 1997 densities are consistent with) is an assumption about the sampling
support, never applied to the observations; zone SHARES do not depend on it, the mean densities do.

Statistics (identical operator on both sides):
  * zone shares per species = RELATIVE DENSITY ACROSS ZONES: mean density per quadrat per zone,
    normalised over zones (`shares_from_densities`); uncertainty by resampling quadrats within
    zones, one random replicate run per bootstrap draw (`bootstrap_shares`);
  * `mean_density_all_quadrats`: the mean density over ALL quadrats of the year (observed with its
    quadrat-bootstrap interval; simulated per replicate run), the supply-calibration statistic.
Without this operator the model is read as a census of ~26,000 recruits per run and its replicate
range is 1-2 share points, an uncertainty of a different kind from the survey's (35-228 quadrats
with clustered densities).

Densities are cached in {out_dir}/virtual_survey_{year}.npz (rerun with --refresh after new model
runs; the cache is rebuilt automatically when its species axis, species names, quadrat ids, run
list or FOOTPRINT differ from the current ones).

Usage: python virtual_survey.py [--year 1999] [--refresh] [--footprint 25] [--caches GLOB] [--out-dir DIR]
       (prints the shares and intervals and the mean density over all quadrats)
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import cascade_metrics as cm  # noqa: E402
from tidal_dispersal_model import SPECIES as _MODEL_SPECIES  # noqa: E402
from transect_geometry import MODEL_QUADRAT_FOOTPRINT_M2 as _TG_FOOTPRINT  # noqa: E402

SPECIES = list(_MODEL_SPECIES)   # 2026-09-20: 4 species (Sonneratia); arrays are sized by len(SPECIES)
ZONES = ['C', 'B', 'A']
# ASSUMED model-side sampling footprint per survey quadrat (m2): the paper's 5 x 5 m quadrat per
# PCQM quadrant (Dahdouh-Guebas et al. 2004, p. 81); the juvenile quadrat size is not documented.
# Sensitivity: 12.5. Shared with transect_geometry so every model-side density uses one assumption.
MODEL_QUADRAT_FOOTPRINT_M2 = float(_TG_FOOTPRINT)
CACHE_GLOB = os.path.join(SCRIPT_DIR, 'results', 'knockout_cascade', 'cascade_D1T1S1E1_seed*.csv')
OUT_DIR = os.path.join(SCRIPT_DIR, 'results', 'evaluation')


def quadrat_radius(footprint=MODEL_QUADRAT_FOOTPRINT_M2):
    """Radius (m) of the circle with the footprint's area."""
    return float(np.sqrt(footprint / np.pi))


QUADRAT_RADIUS = quadrat_radius()   # default footprint, for reference


def load_quadrats(year):
    """Distinct quadrats of the survey year: ids, zones, centres (model frame)."""
    seen = {}
    with open(cm.obs_quadrat_path(year), newline='') as f:
        for r in csv.DictReader(f):
            if r['quadrat_id'] not in seen:
                seen[r['quadrat_id']] = (r['model_zone'], float(r['x']), float(r['y']))
    # 2026-09-21 FIX: drop quadrats outside the model box (26 T5 quadrats); the model cannot place a
    # recruit there, so scoring them against it is a one-sided penalty. See cascade_metrics.
    import tidal_dispersal_model as _t
    ids = [q for q in seen
           if _t.X_MIN <= seen[q][1] <= _t.X_MAX and _t.Y_MIN <= seen[q][2] <= _t.Y_MAX]
    zones = np.array([seen[q][0] for q in ids])
    xy = np.array([[seen[q][1], seen[q][2]] for q in ids])
    return ids, zones, xy


def quadrat_transects(ids):
    """Transect number of each quadrat id ('3.11a' -> 3); 2026-09-20, for the transect hold-out view."""
    return np.array([int(q.split('.')[0]) for q in ids])


def observed_densities(year, ids):
    """Recorded density (#/m2) per quadrat x species (n_quadrats x n_species), same order as ids.
    Species not in SPECIES (Bruguiera) are skipped."""
    tab = defaultdict(dict)
    with open(cm.obs_quadrat_path(year), newline='') as f:
        for r in csv.DictReader(f):
            if r['species'] in SPECIES:
                tab[r['quadrat_id']][r['species']] = float(r['density_per_m2'])
    return np.array([[tab[q].get(sp, 0.0) for sp in SPECIES] for q in ids])


def model_densities(year, refresh=False, caches=None, footprint=MODEL_QUADRAT_FOOTPRINT_M2, out_dir=None):
    """Simulated density (#/m2) per quadrat x species for every replicate run:
    array(n_runs, n_quadrats, n_species) = recruits inside the footprint circle / footprint."""
    ids, zones, xy = load_quadrats(year)
    paths = sorted(glob.glob(caches or CACHE_GLOB))
    out_dir = out_dir or OUT_DIR
    out = os.path.join(out_dir, f'virtual_survey_{year}.npz')
    if os.path.exists(out) and not refresh:
        d = np.load(out, allow_pickle=True)
        # 2026-09-20: a cache built for another species list, run list or footprint is rebuilt
        arr = d['densities'] if 'densities' in d.files else None
        same = (arr is not None and arr.ndim == 3 and arr.shape[2] == len(SPECIES)
                and 'species' in d.files and list(d['species']) == SPECIES
                and 'footprint_m2' in d.files and abs(float(d['footprint_m2']) - footprint) < 1e-9
                and list(d['ids']) == ids and list(d['runs']) == [os.path.basename(p) for p in paths])
        if same:
            return arr, ids, zones
        print(f'  {os.path.basename(out)} is stale (species, runs or footprint differ): rebuilding')
    if not paths:
        sys.exit(f'no run caches match {caches or CACHE_GLOB}. Production figures must read the FINAL run '
                 'cache; this stage does not fall back to any archived or integration-test cache: '
                 'generate the runs, or pass an explicit path.')
    radius = quadrat_radius(footprint)
    runs = []
    for p in paths:
        df = pd.read_csv(p, usecols=['species', 'x', 'y', 'outcome'])
        df = df[df.outcome == 'established']
        dens = np.zeros((len(ids), len(SPECIES)))
        for j, sp in enumerate(SPECIES):
            pts = df.loc[df.species == sp, ['x', 'y']].values
            if len(pts) == 0:
                continue
            tree = cKDTree(pts)
            dens[:, j] = np.array([len(v) for v in tree.query_ball_point(xy, radius)], dtype=float) / footprint
        runs.append(dens)
        print(f'  {os.path.basename(p)}: recruits in the {footprint:g} m2 quadrat circles '
              f'{(dens.sum(axis=0) * footprint).astype(int).tolist()}')
    dens = np.array(runs)
    os.makedirs(out_dir, exist_ok=True)
    np.savez(out, densities=dens, ids=np.array(ids), zones=zones, runs=np.array([os.path.basename(p) for p in paths]),
             species=np.array(SPECIES), footprint_m2=np.array(footprint))
    return dens, ids, zones


def shares_from_densities(dens, zones):
    """Zone shares per species from per-quadrat densities (n_quadrats x n_species): mean density
    per zone, normalised over zones (relative density across zones). Returns array(3 zones x
    n_species); a species without any density (Sonneratia in a run or draw) gets zeros."""
    mean = np.array([dens[zones == z].mean(axis=0) if (zones == z).any() else np.zeros(len(SPECIES))
                     for z in ZONES])
    tot = mean.sum(axis=0)
    return np.where(tot > 0, mean / np.where(tot > 0, tot, 1), 0.0)


def bootstrap_shares(dens_runs, zones, n_boot, rng):
    """Joint bootstrap: for each draw a random replicate run and quadrats resampled within zones.
    dens_runs: array(n_runs, n_quadrats, n_species) (a single observed table may be passed as (1, n, n_species)).
    Returns array(n_boot, 3 zones, n_species)."""
    idx = {z: np.where(zones == z)[0] for z in ZONES}
    out = np.empty((n_boot, len(ZONES), len(SPECIES)))
    for b in range(n_boot):
        c = dens_runs[rng.integers(0, len(dens_runs))]
        mean = np.array([c[rng.choice(idx[z], len(idx[z]))].mean(axis=0) if len(idx[z]) else np.zeros(len(SPECIES))
                         for z in ZONES])
        tot = mean.sum(axis=0)
        out[b] = np.where(tot > 0, mean / np.where(tot > 0, tot, 1), 0.0)
    return out


def _bootstrap_means(dens_runs, zones, n_boot, rng):
    """Bootstrap of the mean density over all quadrats (quadrats resampled within zones, one
    random run per draw): array(n_boot, n_species)."""
    idx = [np.where(zones == z)[0] for z in ZONES if (zones == z).any()]
    out = np.empty((n_boot, len(SPECIES)))
    for b in range(n_boot):
        c = dens_runs[rng.integers(0, len(dens_runs))]
        pick = np.concatenate([rng.choice(ix, len(ix)) for ix in idx])
        out[b] = c[pick].mean(axis=0)
    return out


def virtual_survey(year, n_boot=10000, seed=42, refresh=False, transects=None, footprint=MODEL_QUADRAT_FOOTPRINT_M2,
                   caches=None, out_dir=None):
    """Observed and simulated zone shares (relative density across zones) with quadrat-bootstrap
    uncertainty for one survey year.
    Returns dict with 'obs' {sp: (est(3), lo(3), hi(3))}, 'sim' {sp: (median(3), lo(3), hi(3))},
    'obs_boot' {sp: array(n_boot, 3)}, 'sim_boot' {sp: array(n_boot, 3)}, 'n_runs', 'n_quadrats',
    'footprint_m2', and 'mean_density' {sp: {'obs': (mean, lo, hi), 'sim': (median, lo, hi)}} = the
    mean density over all selected quadrats (observed: bootstrap interval; simulated: median and
    range over the replicate runs).
    transects (2026-09-20): restrict to the quadrats of these transect numbers, e.g. [1, 3] for the
    transects withheld from the stand reconstruction (juvenile hold-out view); None = all."""
    rng = np.random.default_rng(seed)
    dens, ids, zones = model_densities(year, refresh=refresh, caches=caches, footprint=footprint, out_dir=out_dir)
    obs = observed_densities(year, ids)
    if transects is not None:
        keep = np.isin(quadrat_transects(ids), list(transects))
        dens, obs, zones = dens[:, keep, :], obs[keep], zones[keep]
        ids = [q for q, k in zip(ids, keep) if k]
    obs_est = shares_from_densities(obs, zones)
    obs_b = bootstrap_shares(obs[None], zones, n_boot, rng)
    sim_est = np.median([shares_from_densities(c, zones) for c in dens], axis=0)
    sim_b = bootstrap_shares(dens, zones, n_boot, rng)
    obs_mean_b = _bootstrap_means(obs[None], zones, n_boot, rng)
    sim_means = dens.mean(axis=1)                       # n_runs x n_species
    res = {'obs': {}, 'sim': {}, 'obs_boot': {}, 'sim_boot': {}, 'n_runs': len(dens),
           'n_quadrats': {z: int((zones == z).sum()) for z in ZONES}, 'transects': transects,
           'footprint_m2': footprint, 'mean_density': {}}
    for j, sp in enumerate(SPECIES):
        res['obs'][sp] = (obs_est[:, j], np.percentile(obs_b[:, :, j], 2.5, axis=0), np.percentile(obs_b[:, :, j], 97.5, axis=0))
        res['sim'][sp] = (sim_est[:, j], np.percentile(sim_b[:, :, j], 2.5, axis=0), np.percentile(sim_b[:, :, j], 97.5, axis=0))
        res['obs_boot'][sp] = obs_b[:, :, j]
        res['sim_boot'][sp] = sim_b[:, :, j]
        res['mean_density'][sp] = {
            'obs': (float(obs[:, j].mean()), float(np.percentile(obs_mean_b[:, j], 2.5)), float(np.percentile(obs_mean_b[:, j], 97.5))),
            'sim': (float(np.median(sim_means[:, j])), float(sim_means[:, j].min()), float(sim_means[:, j].max()))}
    return res


def mean_density_all_quadrats(year, species, refresh=False, footprint=MODEL_QUADRAT_FOOTPRINT_M2, caches=None,
                              out_dir=None, n_boot=4000, seed=42, transects=None):
    """The calibration statistic of one species: mean density (#/m2) over all quadrats of the
    survey year, observed (recorded densities; 95 % bootstrap interval, quadrats resampled within
    zones) and simulated (virtual survey under the assumed footprint; per replicate run, with the
    median and range). Returns {'obs': (mean, lo, hi), 'sim': (median, lo, hi), 'sim_runs': array,
    'n_quadrats', 'n_nonzero_quadrats_obs', 'footprint_m2', 'n_runs'}."""
    rng = np.random.default_rng(seed)
    dens, ids, zones = model_densities(year, refresh=refresh, caches=caches, footprint=footprint, out_dir=out_dir)
    obs = observed_densities(year, ids)
    if transects is not None:
        keep = np.isin(quadrat_transects(ids), list(transects))
        dens, obs, zones = dens[:, keep, :], obs[keep], zones[keep]
    j = SPECIES.index(species)
    ob = _bootstrap_means(obs[None], zones, n_boot, rng)[:, j]
    sim_runs = dens[:, :, j].mean(axis=1)
    return {'obs': (float(obs[:, j].mean()), float(np.percentile(ob, 2.5)), float(np.percentile(ob, 97.5))),
            'sim': (float(np.median(sim_runs)), float(sim_runs.min()), float(sim_runs.max())),
            'sim_runs': sim_runs, 'n_quadrats': int(len(zones)), 'n_nonzero_quadrats_obs': int((obs[:, j] > 0).sum()),
            'footprint_m2': footprint, 'n_runs': int(len(dens))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--year', default='1999')
    ap.add_argument('--n-boot', type=int, default=10000)
    ap.add_argument('--refresh', action='store_true')
    ap.add_argument('--transects', type=int, nargs='*', default=None, help='restrict to these transects (e.g. 1 3)')
    ap.add_argument('--footprint', type=float, default=MODEL_QUADRAT_FOOTPRINT_M2,
                    help='ASSUMED model-side quadrat footprint in m2 (default 25 = 5 x 5 m; sensitivity 12.5)')
    ap.add_argument('--caches', default=None, help=f'glob of the replicate run caches (default {CACHE_GLOB})')
    ap.add_argument('--out-dir', default=None, help=f'directory of the npz cache (default {OUT_DIR})')
    args = ap.parse_args()
    r = virtual_survey(args.year, args.n_boot, refresh=args.refresh, transects=args.transects, footprint=args.footprint,
                       caches=args.caches, out_dir=args.out_dir)
    print(f'\n{args.year}: {r["n_runs"]} replicate runs, assumed footprint {args.footprint:g} m2 '
          f'(circle radius {quadrat_radius(args.footprint):.3f} m); zone shares = relative density across zones')
    print(f'{"species":<11} {"zone":<5} {"obs":>6} {"[95 % CI]":>16}   {"sim":>6} {"[95 % CI]":>16}')
    for sp in SPECIES:
        for i, z in enumerate(ZONES):
            o, ol, oh = r['obs'][sp]
            s, sl, sh = r['sim'][sp]
            print(f'{sp:<11} {z:<5} {100*o[i]:6.1f} [{100*ol[i]:6.1f}, {100*oh[i]:6.1f}]   '
                  f'{100*s[i]:6.1f} [{100*sl[i]:6.1f}, {100*sh[i]:6.1f}]')
    print(f'\nmean density over all {sum(r["n_quadrats"].values())} quadrats (#/m2): observed [95 % CI] | '
          f'simulated median [range over runs], conditional on the footprint')
    for sp in SPECIES:
        o, s = r['mean_density'][sp]['obs'], r['mean_density'][sp]['sim']
        print(f'{sp:<11} {o[0]:.4f} [{o[1]:.4f}, {o[2]:.4f}] | {s[0]:.4f} [{s[1]:.4f}, {s[2]:.4f}]')


if __name__ == '__main__':
    main()
