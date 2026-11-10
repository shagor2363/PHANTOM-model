"""Vegetation-drag sensitivity on the final calibrated model (2026-09-22).

Puts the zone-specific movement attenuation factor on the same footing as the transport traits
(draft, mobility, lateral spread) that the knockout sweep varies: one factor changed at a time,
everything else at its calibrated value, on the four-transect stand with T2 withheld.

Profiles (C, B, A). The factor is a flow-RETENTION factor, so C < B < A means drag is highest in
zone C and lowest in zone A, the ordering the user requires:

    calibrated    0.4 / 0.7 / 1.0   the fitted profile
    shallow_1     0.6 / 0.8 / 1.0   weaker gradient
    shallow_2     0.8 / 0.9 / 1.0   nearly flat, ordering still intact
    steep         0.2 / 0.6 / 1.0   stronger gradient (the profile used 2026-09-20 to 09-21)
    flat          1.0 / 1.0 / 1.0   NULL: no gradient at all, isolates how much the GRADIENT matters
                                    rather than how much attenuation matters
    pre_reversal  1.0 / 0.6 / 0.2   the pre-2026-09-20 orientation, drag highest LANDWARD

Reported per profile: 1999 zone shares of every species, the 1999 zonation error, the share of
boundary propagules that cross zone C into zone B (the mechanism the profile acts through), the
boundary share of recruits per species, the zone C recruit share and the landward clamp.

The calibrated profile is not re-run: the pipeline's full-model baseline (cascade_D1T1S1E1) has
identical constants and is reused.

Usage:
    python drag_sensitivity.py run   [--seeds 42,49,56] [--workers 12]
    python drag_sensitivity.py report
"""
import argparse
import glob
import json
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TREES = os.path.join(SCRIPT_DIR, 'results', 'reconstruction', 'stand3T_T4T3T1T6_single',
                     'demo_reconstructed_stand3T.csv')
os.environ['KNOCKOUT_TREES'] = TREES

import transport_treatments as cg  # noqa: E402
import tidal_dispersal_model as tdm  # noqa: E402
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results', 'drag_sensitivity')
BASELINE_CACHES = os.path.join(SCRIPT_DIR, 'results', 'knockout_cascade', 'cascade_D1T1S1E1_seed*.csv')

PROFILES = [('steep', (0.2, 0.6, 1.0)), ('shallow_1', (0.6, 0.8, 1.0)), ('shallow_2', (0.8, 0.9, 1.0)),
            ('flat', (1.0, 1.0, 1.0))]


def grid():
    return [(f'drag_{name}_C{c}_B{b}_A{a}', {'ZONE_DRAG': {'C': c, 'B': b, 'A': a}})
            for name, (c, b, a) in PROFILES]


def mechanism(cache):
    """Boundary C->B crossing, boundary share of recruits per species, zone C recruit share."""
    import pandas as pd
    d = pd.read_csv(cache, usecols=['species', 'source', 'outcome', 'zone',
                                    'first_grounding_x', 'first_grounding_y'])
    b = d[d.source == 'boundary']
    z0 = pd.Series([tdm.zone_at(x, y) for x, y in zip(b.first_grounding_x, b.first_grounding_y)], index=b.index)
    c2b = float(((z0 == 'C') & (b.zone == 'B')).sum() / max(int((z0 == 'C').sum()), 1))
    e = d[d.outcome == 'established']
    bf = {sp: float((e[e.species == sp].source == 'boundary').mean()) for sp in tdm.SPECIES
          if (e.species == sp).any()}
    return c2b, float((e.zone == 'C').mean()), bf


def report(seeds):
    import cascade_metrics as cm
    import numpy as np
    import pandas as pd
    obs = {'Avicennia': 0.54, 'Ceriops': 0.80, 'Rhizophora': 0.37, 'Sonneratia': 0.0}
    entries = [('calibrated (0.4/0.7/1.0)', sorted(glob.glob(BASELINE_CACHES))[:len(seeds)])]
    for tag, _ in grid():
        name = tag.split('_', 1)[1]
        entries.append((name, [cg.cache_path(RESULTS_DIR, tag, s, prefix='calib') for s in seeds]))
    rows = []
    for label, caches in entries:
        caches = [c for c in caches if os.path.exists(cm.summary_path_for(c))]
        if not caches:
            print('  missing', label)
            continue
        errs, clamps, shares, c2bs, zcs = [], [], [], [], []
        for c in caches:
            s = json.load(open(cm.summary_path_for(c)))
            sc = cm.score_summary(s, w1997=0.0)
            errs.append(sc['zonation_error_1999']); clamps.append(sc['clamp_frac'])
            shares.append(s['strata']['all']['sampled_shares'])
            a, b, _ = mechanism(c)
            c2bs.append(a); zcs.append(b)
        r = {'profile': label, 'n': len(caches), 'err99': round(float(np.mean(errs)), 4),
             'clamp': round(float(np.mean(clamps)), 3), 'C_to_B': round(float(np.mean(c2bs)), 3),
             'zoneC_recruits': round(float(np.mean(zcs)), 3)}
        for sp in tdm.SPECIES:
            v = [sh[sp]['A'] for sh in shares if sp in sh]
            if v:
                r[f'{sp[:3]}_A'] = round(float(np.mean(v)), 2)
        rows.append(r)
    df = pd.DataFrame(rows)
    pd.set_option('display.width', 250)
    print(df.to_string(index=False))
    print('\nobserved 1999 zone A shares (in-domain quadrats): ' +
          ', '.join(f'{k[:3]} {v}' for k, v in obs.items()))
    print('C_to_B = share of boundary propagules first grounded in zone C that end in zone B;')
    print('this is the mechanism the drag profile acts through, so it should track the zonation response.')
    os.makedirs(RESULTS_DIR, exist_ok=True)
    df.to_csv(os.path.join(RESULTS_DIR, 'drag_sensitivity.csv'), index=False)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['run', 'report'])
    ap.add_argument('--seeds', default='42,49,56')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--days', type=int, default=365)
    a = ap.parse_args()
    seeds = [int(s) for s in a.seeds.split(',')]
    os.makedirs(RESULTS_DIR, exist_ok=True)
    if a.command == 'run':
        g = grid()
        print(f'{len(g)} drag profiles x {len(seeds)} seeds (the calibrated profile reuses the pipeline baseline)')
        cg.run_grid(g, seeds, a.days, RESULTS_DIR, TREES, workers=a.workers, prefix='calib')
    report(seeds)


if __name__ == '__main__':
    main()
