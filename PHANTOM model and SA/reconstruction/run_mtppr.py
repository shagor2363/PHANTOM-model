"""Stage 3 orchestrator: run MTPPR (tile x seed) jobs as parallel Rscript processes.

Usage:
    python reconstruction/run_mtppr.py --geom results/reconstruction/geom \
        --out results/reconstruction/mtppr --seeds 1 2 3 4 5 --workers 12
    python reconstruction/run_mtppr.py ... --tiles 202 303 --max-steps 20000   (subset / pass-through)

Unknown `--key value` options are forwarded to mtppr_tile.R as `--key=value`.
Finished jobs (summary JSON without error) are skipped unless --no-resume.
"""
import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recon_config as C  # noqa: E402

TILE_R = os.path.join(C.RECON_DIR, 'mtppr_tile.R')


def engine_sha256(path=C.ENGINE_R):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


def check_engine(path=C.ENGINE_R):
    """2026-09-20 B0: the vendored engine must be byte-identical to upstream commit 6982c39."""
    got = engine_sha256(path)
    if got != C.ENGINE_SHA256:
        raise RuntimeError(f'MTPPR engine {path} has sha256 {got}, expected {C.ENGINE_SHA256}: the engine must not be modified')
    return got


def build_jobs(geom_dir, out_dir, seeds, tiles=None, resume=True):
    tl = pd.read_csv(os.path.join(geom_dir, 'tiles.csv'))
    if tiles:
        tl = tl[tl['tile_id'].isin(tiles)]
    tl = tl.sort_values('expected_n', ascending=False)
    jobs = []
    for r in tl.to_dict(orient='records'):   # dict access: 'lambda' is a keyword and namedtuples rename it
        for s in seeds:
            stem = os.path.join(out_dir, f'tile{int(r["tile_id"])}_seed{int(s)}')
            if resume and os.path.exists(stem + '_summary.json'):
                with open(stem + '_summary.json') as f:
                    j = json.load(f)
                if 'error' not in j and os.path.exists(stem + '_points.csv'):
                    continue
            # 2026-09-20 B0: the tile's band density (lambda_used) is the count-control target
            jobs.append(dict(tile_id=int(r['tile_id']), seed=int(s), expected_n=float(r['expected_n']),
                             lambda_target=float(r.get('lambda', 0.0)),
                             expected_n_eligible=float(r.get('expected_n_eligible', r['expected_n']))))
    return jobs


CONFIG_ARGS = {
    'max-steps': C.MAX_STEPS, 'min-steps': C.MIN_STEPS, 'steps-per-point': C.STEPS_PER_POINT,
    'conv-start': C.CONV_START, 'conv-window': C.CONV_WINDOW, 'conv-tol': C.CONV_TOL,
    'trace-every': C.TRACE_EVERY, 'rmax': C.R_MAX, 'dr': C.DR, 'bw': C.BW, 'margin': C.TILE_MARGIN,
}


def engine_args(extra):
    """Config defaults for mtppr_tile.R, overridden by user-supplied --key=value extras."""
    merged = {k: str(v) for k, v in CONFIG_ARGS.items()}
    for a in extra:
        k, _, v = a.lstrip('-').partition('=')
        merged[k] = v if v != '' else '1'
    return [f'--{k}={v}' for k, v in merged.items()]


def count_control_args(job, count_control=None):
    """2026-09-20 B0 ADAPTATION: pinned count at the tile's band density (C.COUNT_CONTROL)."""
    cc = C.COUNT_CONTROL if count_control is None else count_control
    if not cc:
        return []
    lam = job.get('lambda_target', 0.0)
    assert lam > 0, f'tile {job["tile_id"]}: count control needs a positive lambda_target'
    args = ['--pin-n=1', f'--lambda-target={lam:.10g}']
    if C.REF_MOMENT_RESCALING:
        args.append('--ref-density-scale=1')   # reference moment rescaling (count-control adaptation), see recon_config
    return args


def run_job(job, geom_dir, out_dir, log_dir, extra, count_control=None):
    check_engine()
    cmd = [C.RSCRIPT, '--vanilla', TILE_R, f'--tile={job["tile_id"]}', f'--seed={job["seed"]}',
           f'--geom={geom_dir}', f'--out={out_dir}', '--quiet=1'] + engine_args(extra) + count_control_args(job, count_control)
    log = os.path.join(log_dir, f'tile{job["tile_id"]}_seed{job["seed"]}.log')
    t0 = time.time()
    rc = -1
    for attempt in (1, 2):
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            rc = p.returncode
            with open(log, 'a') as f:
                f.write(f'--- attempt {attempt} rc={rc} ---\n{p.stdout}\n{p.stderr}\n')
        except subprocess.TimeoutExpired:
            rc = -9
            with open(log, 'a') as f:
                f.write(f'--- attempt {attempt} timeout ---\n')
        if rc == 0:
            break
    return dict(job, returncode=rc, seconds=time.time() - t0, command=' '.join(cmd))


def run_jobs(jobs, geom_dir, out_dir, log_dir, n_workers=C.N_WORKERS, extra=(), count_control=None):
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    results = []
    if not jobs:
        return results
    sha = check_engine()
    with open(os.path.join(out_dir, 'engine_sha256.txt'), 'w') as f:
        f.write(f'{sha}  {os.path.relpath(C.ENGINE_R, C.ROOT)}\n')
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(run_job, j, geom_dir, out_dir, log_dir, extra, count_control): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            results.append(r)
            flag = 'ok' if r['returncode'] == 0 else f'FAILED rc={r["returncode"]}'
            print(f'[{i}/{len(jobs)}] tile {r["tile_id"]} seed {r["seed"]} {flag} ({r["seconds"]:.0f} s, '
                  f'{time.time() - t0:.0f} s elapsed)', flush=True)
    return results


def collect_summaries(out_dir):
    rows, copies = [], []
    for path in sorted(glob.glob(os.path.join(out_dir, 'tile*_seed*_summary.json'))):
        with open(path) as f:
            j = json.load(f)
        if 'error' in j:
            continue
        # 2026-09-20 superposition: per-copy summaries (key 'copy' >= 1) are listed separately; the
        # combined tile summaries (superpose.superpose_tile) carry copy 1's energies
        if int(j.get('copy', 0) or 0) >= 1:
            copies.append(j)
            continue
        rows.append({k: v for k, v in j.items() if k != 'copies'})
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(['tile_id', 'seed']).reset_index(drop=True)
        df.to_csv(os.path.join(out_dir, 'mtppr_summaries.csv'), index=False)
    if copies:
        pd.DataFrame(copies).sort_values(['tile_id', 'seed', 'copy']).to_csv(os.path.join(out_dir, 'mtppr_copy_summaries.csv'), index=False)
    return df


def select_best(summaries):
    """Lowest final energy per tile; ties broken by count closest to the expected count."""
    best = {}
    for tid, g in summaries.groupby('tile_id'):
        g = g.assign(dn=(g['n_final'] - g['expected_n']).abs())
        g = g.sort_values(['e_final', 'dn'])
        best[int(tid)] = int(g['seed'].iloc[0])
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--geom', default=C.GEOM_DIR)
    ap.add_argument('--out', default=C.MTPPR_DIR)
    ap.add_argument('--logs', default=None)
    ap.add_argument('--seeds', type=int, nargs='+', default=list(range(1, C.N_SEEDS + 1)))
    ap.add_argument('--tiles', type=int, nargs='*', default=None)
    ap.add_argument('--workers', type=int, default=C.N_WORKERS)
    ap.add_argument('--no-resume', action='store_true')
    ap.add_argument('--no-count-control', action='store_true', help='emergent count (engine default) instead of C.COUNT_CONTROL')
    args, unknown = ap.parse_known_args()
    extra = []
    i = 0
    while i < len(unknown):
        a = unknown[i]
        if a.startswith('--') and '=' not in a and i + 1 < len(unknown) and not unknown[i + 1].startswith('--'):
            extra.append(f'{a}={unknown[i + 1]}')
            i += 2
        else:
            extra.append(a)
            i += 1
    log_dir = args.logs or os.path.join(os.path.dirname(os.path.abspath(args.out)), 'logs')
    jobs = build_jobs(args.geom, args.out, args.seeds, args.tiles, resume=not args.no_resume)
    print(f'{len(jobs)} jobs to run with {args.workers} workers; extra args: {extra}')
    results = run_jobs(jobs, args.geom, args.out, log_dir, args.workers, extra,
                       count_control=(False if args.no_count_control else None))
    failed = [r for r in results if r['returncode'] != 0]
    summ = collect_summaries(args.out)
    if len(summ):
        best = select_best(summ)
        with open(os.path.join(args.out, 'selection.json'), 'w') as f:
            json.dump(best, f, indent=2)
        cols = ['tile_id', 'seed', 'n_fixed', 'n_final', 'expected_n', 'e_initial', 'e_final', 'steps', 'stop_reason', 'seconds']
        print(summ[cols].round(4).to_string(index=False))
        print('best seeds:', best)
    if failed:
        print(f'{len(failed)} jobs failed: {[(r["tile_id"], r["seed"]) for r in failed]}')
        sys.exit(1)


if __name__ == '__main__':
    main()
