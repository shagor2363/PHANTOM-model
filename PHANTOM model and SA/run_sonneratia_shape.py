"""Sonneratia 'shape runs': the calibrated model with Sonneratia's supply amplified, for drawing its
terminal-elevation distributions at a usable sample size (2026-09-23).

Why. Sonneratia's supply is FITTED to the mean recorded 1999 juvenile density (0.02 % of the juvenile
density), which puts its boundary influx at 0.13 propagules per spring day: ~30 boundary and ~70
parent-tree propagules per year, 29 boundary strandings over three seeds. No distribution can be drawn
from that, and the user ruled a Sonneratia-free boundary panel unacceptable to reviewers.

Why this is legitimate. Propagules are passive and do not interact (ODD 2.3): drag is a zone property,
establishment is evaluated per propagule with no density dependence, the substrate is per cell. The
terminal-elevation distribution of a source stratum, and its recruit share, are therefore independent
of the supply rate; the rate only sets the sample size. Amplifying Sonneratia's supply by SON_FACTOR
changes nothing for the other three species and nothing in Sonneratia's shape, only its counts. These
caches must NEVER be used for densities, recruit fractions or the observed-vs-simulated comparison;
they serve the distribution figures only (endpoints grid, stranded-to-recruits), and captions must say
so.

Runs: treatments 'baseline' (calibrated transport) and 'all_transport_eq' (the two columns of the
endpoints grid), seeds 42/49/56, 365 d, on the pipeline's stand, cached in results/knockout_son_shape/
with the same file names as results/knockout_final so that the figure scripts can swap the Sonneratia
rows in (`--son-shape-dir`).

Usage: python run_sonneratia_shape.py [--seeds 42,49,56] [--workers 6] [--factor 1000]
"""
import argparse
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
STAND = os.path.join(SCRIPT_DIR, 'results', 'reconstruction', 'stand3T_T4T3T1T6_single',
                     'demo_reconstructed_stand3T.csv')
OUT_DIR = os.path.join(SCRIPT_DIR, 'results', 'knockout_son_shape')
SON_FACTOR = 1000.0          # 0.1344 -> 134.4 per spring day (Rhizophora: 160); F_BASE 0.00384 -> 3.84
TREATMENTS_USED = ['baseline', 'all_transport_eq']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--seeds', default='42,49,56')
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--factor', type=float, default=SON_FACTOR)
    ap.add_argument('--days', type=int, default=365)
    ap.add_argument('--trees', default=STAND)
    args = ap.parse_args()

    # the stand must be fixed BEFORE the knockout module is imported (it reads KNOCKOUT_TREES at import;
    # the worker processes inherit the environment)
    os.environ['KNOCKOUT_TREES'] = os.path.abspath(args.trees)
    sys.path.insert(0, SCRIPT_DIR)
    import transport_treatments as tt
    import tidal_dispersal_model as tdm

    son = 'Sonneratia'
    supply = {
        'BOUNDARY_SPRING': {**tdm.BOUNDARY_SPRING, son: tdm.BOUNDARY_SPRING[son] * args.factor},
        'BOUNDARY_NEAP':   {**tdm.BOUNDARY_NEAP,   son: tdm.BOUNDARY_NEAP[son] * args.factor},
        'F_BASE':          {**tdm.F_BASE,          son: tdm.F_BASE[son] * args.factor},
    }
    print(f"Sonneratia supply x{args.factor:g}: boundary {supply['BOUNDARY_SPRING'][son]:.1f} / "
          f"{supply['BOUNDARY_NEAP'][son]:.1f} per day, F_BASE {supply['F_BASE'][son]:.3f}")
    for t in TREATMENTS_USED:
        tt.TREATMENTS[t] = {**tt.TREATMENTS[t], **supply}   # overrides travel to the workers as arguments
    os.makedirs(OUT_DIR, exist_ok=True)
    seeds = [int(s) for s in args.seeds.split(',')]
    tt.run_treatments(TREATMENTS_USED, seeds, args.days, OUT_DIR, max_workers=args.workers)
    with open(os.path.join(OUT_DIR, 'README.txt'), 'w') as f:
        f.write(f"Sonneratia shape runs (2026-09-23): supply x{args.factor:g} for Sonneratia only; stand {args.trees};\n"
                f"treatments {TREATMENTS_USED}; seeds {seeds}; {args.days} d. Distribution figures only; never densities.\n")
    print('done')


if __name__ == '__main__':
    main()
