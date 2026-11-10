"""Run everything: model runs of the experiments and the figures.

    python run_all.py                      all runs and figures (128 one-year runs, about 6 GB in results/)
    python run_all.py --rebuild-inputs     first rebuild the shipped inputs from the raw data (needs R, about 15 min)
    python run_all.py --figures-only       redraw the figures from existing runs
    python run_all.py --days 60 --workers 4    quick check that the chain executes (figures not meaningful)

Stages (sections of the methods in brackets):
  A  calibrated model, 20 seeds                                        [evaluation, establishment]
  B  transport traits equalised: all three, and each one alone         [2.4.2 deconstruction]
  C  one-at-a-time sweep of draft, mobility and lateral spread         [2.4.2 OAT]
  D  Sonneratia supply x 1000, used for distribution curves only       [2.4.2]
  E  vegetation-drag profiles                                          [2.4.1]
  F  figures -> results/figures
"""
import argparse
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
STAND = os.path.join(ROOT, 'results', 'reconstruction', 'stand3T_T4T3T1T6_single', 'demo_reconstructed_stand3T.csv')
SEEDS = [42, 49, 56]
SEEDS20 = list(range(42, 176, 7))


def py(*args):
    """Run one script of the package and stop at the first failure."""
    print('\n>>> python ' + ' '.join(args), flush=True)
    subprocess.run([sys.executable, *args], cwd=ROOT, check=True)


def rebuild_inputs(workers):
    py('extract_juvenile.py', '1997')                      # workbook -> juvenile quadrat tables
    py('extract_juvenile.py', '1999')
    py('zone_quadrats_by_position.py')                     # quadrats zoned by position in the model frame
    py('reconstruction/make_substrate_map.py')             # forest floor / treeless band
    py('reconstruction/run_stand3T.py', '--retain', '4', '3', '1', '6', '--band-breaks', '0', '60', '120',
       '--members', '5', '--method', 'single', '--density-mode', 'recorded', '--tag', '_single',
       '--skip-validate', '--workers', str(workers))       # adult stand, T2 withheld


def runs(days, workers):
    import transport_treatments as tt
    seeds, w, d = ','.join(map(str, SEEDS)), str(workers), str(days)
    cascade = os.path.join(ROOT, 'results', 'knockout_cascade')
    final = os.path.join(ROOT, 'results', 'knockout_final')
    os.makedirs(final, exist_ok=True)
    print('\n=== A. calibrated model, 20 seeds', flush=True)
    tt.run_grid([('D1T1S1E1', {})], SEEDS20, days, cascade, STAND, workers=workers, prefix='cascade')
    for s in SEEDS:
        shutil.copyfile(os.path.join(cascade, f'cascade_D1T1S1E1_seed{s}.csv'),
                        os.path.join(final, f'knockout_baseline_seed{s}.csv'))
    print('\n=== B. transport traits equalised', flush=True)
    tt.run_treatments(['all_transport_eq', 'draft_eq_mean', 'mobility_eq_mean', 'spread_eq_mean'],
                      SEEDS, days, final, max_workers=workers)
    print('\n=== C. one-at-a-time sweep', flush=True)
    py('equal_cv_experiment.py', 'run', '--seeds', seeds, '--workers', w, '--days', d)
    print('\n=== D. Sonneratia supply x 1000', flush=True)
    py('run_sonneratia_shape.py', '--seeds', seeds, '--workers', w, '--days', d, '--trees', STAND)
    print('\n=== E. vegetation-drag profiles', flush=True)
    py('drag_sensitivity.py', 'run', '--seeds', seeds, '--workers', w, '--days', d)


def figures():
    print('\n=== F. figures', flush=True)
    os.makedirs(os.path.join(ROOT, 'results', 'figures'), exist_ok=True)
    os.makedirs(os.path.join(ROOT, 'results', 'evaluation'), exist_ok=True)
    py('fig_layout_stand.py', '--out', os.path.join('results', 'figures'))
    py('figure_endpoints_main.py')
    py('figure_establishment_filter.py')
    py('equal_cv_experiment.py', 'figure', '--seeds', ','.join(map(str, SEEDS)))
    py('figure_drag_zonation.py')
    py('virtual_survey.py', '--year', '1999', '--refresh')
    py('virtual_survey.py', '--year', '1997', '--refresh')
    py('figure_evaluation_stacked.py')
    py('evaluate_recruits_vs_adults.py', '--withheld', '2')
    py('figure_evaluation_1999_change_adults.py')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workers', type=int, default=12, help='parallel model runs')
    ap.add_argument('--days', type=int, default=365, help='length of every model run')
    ap.add_argument('--rebuild-inputs', action='store_true', help='rebuild juvenile tables, substrate map and stand from the raw data')
    ap.add_argument('--figures-only', action='store_true', help='skip the model runs')
    args = ap.parse_args()
    os.environ['KNOCKOUT_TREES'] = STAND     # the parent-tree inventory of every run
    os.environ['MPLBACKEND'] = 'Agg'
    if args.rebuild_inputs:
        rebuild_inputs(args.workers)
    if not args.figures_only:
        runs(args.days, args.workers)
    figures()
    print('\n=== done: figures in results/figures')


if __name__ == '__main__':
    main()
