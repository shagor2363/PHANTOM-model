"""Three-transect, band-wise MTPPR stand ensemble (2026-09-20 B0/B1; replaces run_holdout.py for
this design, run_holdout.py is kept).

Design (plan 'and-with-the-evaluation-synthetic-octopus', Phase B): the adult stand of the whole
domain is reconstructed from the RETAINED parallel transects only (default T4, T2, T6; T3 and T1
withheld, T5 never used).  Corridors are the midlines between the retained transects, bands are
fringe-distance bands with explicit breaks (default [0, 60, 120] m, last band open), the band
density is the censored PCQM maximum-likelihood estimate (prep_pcqm.band_density, DENSITY_MODE
'pcqm_mle') and every tile is reconstructed with Chris Wudel's unmodified MTPPR v2.0.0 engine
under the COUNT-CONTROL adaptation (mtppr_tile.R: pinned count = lambda_band x tile area,
add/delete proposals off; verified by b0_check.py before any ensemble run).  Each ensemble member
k has its own anchor realisation (prep seed 42 + k - 1) paired with its own MTPPR seed k, so the
ensemble carries positional and optimiser uncertainty; the representative stand is the member with
the lowest total final energy (predefined criterion).

Usage:
    python reconstruction/run_stand3T.py --retain 4 2 6 --band-breaks 0 60 120 --members 5 --workers 12
    python reconstruction/run_stand3T.py --members 1 --prep-only              # geometry of member 1 only
    python reconstruction/run_stand3T.py --members 1 --only-tile 201          # one tile (checks)
    python reconstruction/run_stand3T.py --band-breaks 0 80 --tag _2band      # two-band sensitivity
    python reconstruction/run_stand3T.py --drop-duplicates first --tag _dupfirst

Outputs under results/reconstruction/stand3T_T4T2T6{tag}/:
    geom_m{k}/                        prep of member k (records, anchors, band_density, tiles, ...)
    mtppr_m{k}/                       MTPPR tiles of member k (seed k)
    reconstructed_stand_m{k}.csv      post-processed stand of member k (--clip-treeless)
    demo_reconstructed_m{k}.csv       demo.csv schema
    reconstructed_stand_stand3T.csv, demo_reconstructed_stand3T.csv   copies of the representative member
    stand_meta.json                   settings, per-member energies, realised counts per tile / band vs
                                      targets, representative_member
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recon_config as C  # noqa: E402
import prep_pcqm  # noqa: E402
import run_mtppr  # noqa: E402
import postprocess as P  # noqa: E402


def stand_tag(retain):
    order = [t for t in C.CORRIDOR_ORDER if t in retain]
    return ''.join(f'T{t}' for t in order)


def realised_counts(df, tiles, dens=None, summaries=None):
    """Realised (post-processed, clipped) count per tile and band against the eligible-area targets,
    the target's plot-bootstrap 95 % interval (lambda_MLE CI x modelled share; `dens` = band_density.csv)
    and, for the superposition, the cross-copy close-pair diagnostic per tile (`summaries`)."""
    n_tile = df.groupby('tile_id').size()
    ci = {}
    if dens is not None:
        for d in dens.to_dict(orient='records'):
            ci[int(d['band'])] = (float(d['lambda_mle_lo'] * d['share_modelled']), float(d['lambda_mle_hi'] * d['share_modelled']))
    cross = {}
    if summaries is not None and 'n_cross_copy_pairs_lt_hardcore' in summaries.columns:
        for s in summaries.to_dict(orient='records'):
            cross[int(s['tile_id'])] = (s.get('n_cross_copy_pairs_lt_hardcore'), s.get('share_trees_in_cross_copy_pairs'), s.get('k'), s.get('n_copies'), s.get('p_keep_last'))
    rows = []
    for t in tiles.to_dict(orient='records'):   # dict access: 'lambda' is a keyword, namedtuples rename it
        n = int(n_tile.get(int(t['tile_id']), 0))
        lo, hi = ci.get(int(t['band']), (np.nan, np.nan))
        dens_e = n / t['eligible_area'] if t['eligible_area'] > 0 else np.nan
        cr = cross.get(int(t['tile_id']), (np.nan,) * 5)
        rows.append(dict(tile_id=int(t['tile_id']), transect=int(t['transect']), band=int(t['band']), area=float(t['area']),
                         eligible_area=float(t['eligible_area']), lambda_target=float(t['lambda']), lambda_lo=lo, lambda_hi=hi,
                         expected_n_eligible=float(t['expected_n_eligible']), n_realised=n,
                         ratio=(n / t['expected_n_eligible'] if t['expected_n_eligible'] > 0 else np.nan),
                         deviation_pct=(100 * (n / t['expected_n_eligible'] - 1) if t['expected_n_eligible'] > 0 else np.nan),
                         density_eligible=dens_e, within_target_ci=bool(lo <= dens_e <= hi) if np.isfinite(lo) else None,
                         n_cross_copy_pairs_lt_hardcore=cr[0], share_trees_in_cross_copy_pairs=cr[1], k=cr[2], n_copies=cr[3], p_keep_last=cr[4]))
    tl = pd.DataFrame(rows)
    bands = tl.groupby('band').agg(area=('area', 'sum'), eligible_area=('eligible_area', 'sum'),
                                   expected_n_eligible=('expected_n_eligible', 'sum'), n_realised=('n_realised', 'sum'),
                                   lambda_target=('lambda_target', 'first'), lambda_lo=('lambda_lo', 'first'), lambda_hi=('lambda_hi', 'first'),
                                   k=('k', 'first'), n_copies=('n_copies', 'first')).reset_index()
    bands['ratio'] = bands['n_realised'] / bands['expected_n_eligible']
    bands['deviation_pct'] = 100 * (bands['ratio'] - 1)
    bands['density_eligible'] = bands['n_realised'] / bands['eligible_area']
    bands['within_target_ci'] = (bands['density_eligible'] >= bands['lambda_lo']) & (bands['density_eligible'] <= bands['lambda_hi'])
    return tl, bands


def run_member(k, out, retain, band_breaks, prep_seed, mtppr_seed, args):
    geom_k = os.path.join(out, f'geom_m{k}')
    mtppr_k = os.path.join(out, f'mtppr_m{k}')
    excl = [t for t in C.STAND_TRANSECTS if t not in retain]
    if not os.path.exists(os.path.join(geom_k, 'tiles.csv')) or args.redo_prep:
        prep_pcqm.run(geom_k, seed=prep_seed, exclude_transect=excl, check=False, verbose=False,
                      density_mode=args.density_mode, band_breaks=band_breaks, drop_duplicates=args.drop_duplicates)
    with open(os.path.join(geom_k, 'prep_meta.json')) as f:
        pm = json.load(f)
    assert sorted(pm['transects']) == sorted(retain), pm['transects']
    rec = pd.read_csv(os.path.join(geom_k, 'records.csv'))
    assert set(rec['transect']) == set(retain), f'records of other transects in member {k}'
    print(f'member {k}: prep seed {prep_seed}, MTPPR seed {mtppr_seed}, anchors {pm["n_anchors"]} '
          f'(by transect {pm["n_anchors_by_transect"]}), adults {pm["n_adults"]}, plots dropped {pm["duplicates"]["plots_dropped"]}',
          flush=True)
    member = dict(member=k, prep_seed=prep_seed, mtppr_seed=mtppr_seed, geom=os.path.relpath(geom_k, C.ROOT),
                  mtppr=os.path.relpath(mtppr_k, C.ROOT), n_anchors=pm['n_anchors'], n_adults=pm['n_adults'],
                  n_anchors_by_transect=pm['n_anchors_by_transect'], plots_dropped=pm['duplicates']['plots_dropped'])
    if args.prep_only:
        return member
    tiles = pd.read_csv(os.path.join(geom_k, 'tiles.csv'))
    run_ids = args.only_tile or tiles[tiles['n_ref'] > 0]['tile_id'].tolist()
    if args.method == 'superposition':
        import superpose as SP   # rejected superposition scenario; superpose.py is not part of this package
        # 2026-09-20 candidate: emergent-count copies superposed (superpose.py); k per band recorded
        if not args.skip_mtppr:
            _, kt = SP.superpose_geom(geom_k, mtppr_k, mtppr_seed, os.path.join(out, 'logs'), args.workers, tiles=run_ids)
        else:
            kt = SP.k_table(geom_k)
        member['k_per_band'] = kt[['band', 'lambda_target', 'lambda_belt', 'k_lo', 'k_hat', 'k_hi', 'n_copies']].round(4).to_dict(orient='records')
    elif not args.skip_mtppr:
        jobs = run_mtppr.build_jobs(geom_k, mtppr_k, [mtppr_seed], tiles=run_ids, resume=True)
        print(f'member {k}: {len(jobs)} MTPPR jobs (tiles {run_ids})', flush=True)
        res = run_mtppr.run_jobs(jobs, geom_k, mtppr_k, os.path.join(out, 'logs'), args.workers,
                                 count_control=(False if args.no_count_control else None))
        failed = [r for r in res if r['returncode'] != 0]
        if failed:
            raise RuntimeError(f'member {k}: {len(failed)} failed jobs {[(r["tile_id"], r["seed"]) for r in failed]}')
    summ = run_mtppr.collect_summaries(mtppr_k)
    summ = summ[summ['seed'] == mtppr_seed]
    if 'copy' in summ.columns:
        # keep the combined tile summaries only. Superposition writes them with 'copy' missing (the
        # per-copy rows carry copy >= 1); the single-realisation path writes copy = 0. BUG FIXED
        # 2026-09-21: the old test kept only NaN, so every single-realisation summary was dropped,
        # leaving selection.json empty, the energies 0.0 and the representative member defaulted to
        # member 1 instead of the lowest-energy one (affected stand3T_T4T2T6_single as well).
        keep = summ['copy'].isna() | (summ['copy'].fillna(0).astype(int) == 0)
        if 'method' in summ.columns:
            keep = keep | (summ['method'] == 'superposition')
        summ = summ[keep]
    with open(os.path.join(mtppr_k, 'selection.json'), 'w') as f:   # one seed per member (validate.py reads it)
        json.dump({str(int(t)): int(mtppr_seed) for t in summ['tile_id']}, f, indent=2)
    member['total_energy_final'] = float(summ['e_final'].sum())
    member['total_energy_initial'] = float(summ['e_initial'].sum())
    member['tiles_run'] = sorted(int(t) for t in summ['tile_id'])
    member['count_control'] = bool(summ['count_control'].all()) if 'count_control' in summ.columns else False
    member['ref_moment_rescaling'] = bool((summ['ref_density_scale'] > 0).all()) if 'ref_density_scale' in summ.columns else False
    if args.only_tile:
        return member
    geom = P.load_geom(geom_k)
    df, meta = P.assemble(geom, mtppr_k, seed=mtppr_seed, clip_treeless=True)
    df.to_csv(os.path.join(out, f'reconstructed_stand_m{k}.csv'), index=False)
    df[P.DEMO_COLS].to_csv(os.path.join(out, f'demo_reconstructed_m{k}.csv'), index=False)
    dens = pd.read_csv(os.path.join(geom_k, 'band_density.csv'))
    tl, bands = realised_counts(df, tiles, dens, summ)
    tl['within_10pct'] = (tl['ratio'] - 1).abs() <= 0.10
    bands['within_10pct'] = (bands['ratio'] - 1).abs() <= 0.10
    tl.to_csv(os.path.join(out, f'realised_counts_m{k}.csv'), index=False)
    member['postprocess'] = meta
    member['realised_tiles'] = tl.to_dict(orient='records')
    member['realised_bands'] = bands.to_dict(orient='records')
    member['n_trees'] = int(len(df))
    member['by_species'] = df['species'].value_counts().to_dict()
    # per-member validation (band tables from meta, marks KS, summary functions, seams, zones)
    if not args.skip_validate:
        cmd = [sys.executable, os.path.join(C.RECON_DIR, 'validate.py'), '--recon-dir', out, '--geom', geom_k, '--mtppr', mtppr_k,
               '--stand', os.path.join(out, f'reconstructed_stand_m{k}.csv'), '--tag', f'_m{k}']
        p = subprocess.run(cmd, capture_output=True, text=True)
        with open(os.path.join(out, 'logs', f'validate_m{k}.log'), 'w') as f:
            f.write(p.stdout + '\n' + p.stderr)
        member['validate_rc'] = p.returncode
        fit = os.path.join(out, f'tables_m{k}', 'summary_function_fit.csv')
        if os.path.exists(fit):
            member['summary_function_fit'] = pd.read_csv(fit).groupby('band')[['g_mad', 'kmm_r_mad', 'kmm_h_mad', 'mingling_mad', 'G_mad']].mean().round(4).reset_index().to_dict(orient='records')
        ks = os.path.join(out, f'tables_m{k}', 'species_marks_ks.csv')
        if os.path.exists(ks):
            member['marks_ks'] = pd.read_csv(ks)[['species', 'mark', 'ks_stat', 'ks_p']].round(4).to_dict(orient='records')
    return member


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--recon-dir', default=C.OUT_ROOT)
    ap.add_argument('--retain', type=int, nargs='+', default=C.RETAINED_TRANSECTS_3T)
    ap.add_argument('--band-breaks', type=float, nargs='+', default=C.BAND_BREAKS_3T)
    ap.add_argument('--members', type=int, default=C.N_SEEDS)
    ap.add_argument('--prep-seeds', type=int, nargs='+', default=None, help='anchor seeds per member (default 42, 43, ...)')
    ap.add_argument('--mtppr-seeds', type=int, nargs='+', default=None, help='MTPPR seeds per member (default 1, 2, ...)')
    ap.add_argument('--drop-duplicates', choices=list(prep_pcqm.DROP_DUPLICATE_MODES), default='none')
    ap.add_argument('--density-mode', choices=list(C.DENSITY_MODES), default=None)
    ap.add_argument('--workers', type=int, default=C.N_WORKERS)
    ap.add_argument('--only-tile', type=int, nargs='*', default=None, help='reconstruct only these tiles (checks)')
    ap.add_argument('--skip-mtppr', action='store_true', help='reuse the tiles, only post-process')
    ap.add_argument('--prep-only', action='store_true', help='prep geometry only (no MTPPR, no assembly)')
    ap.add_argument('--redo-prep', action='store_true')
    ap.add_argument('--no-count-control', action='store_true')
    ap.add_argument('--method', choices=['superposition', 'single'], default='superposition' if C.SUPERPOSITION else 'single',
                    help='superposition candidate (superpose.py) or one MTPPR run per tile (count control per recon_config)')
    ap.add_argument('--skip-validate', action='store_true')
    ap.add_argument('--tag', default='', help="suffix of the output dir, e.g. '_2band', '_dupfirst'")
    args = ap.parse_args()
    retain = sorted(set(args.retain))
    for t in retain:
        assert t in C.STAND_TRANSECTS, f'T{t} is not a stand transect'
    out = os.path.join(args.recon_dir, f'stand3T_{stand_tag(retain)}{args.tag}')
    os.makedirs(out, exist_ok=True)
    prep_seeds = args.prep_seeds or [C.ANCHOR_SEED + i for i in range(args.members)]
    mtppr_seeds = args.mtppr_seeds or [i + 1 for i in range(args.members)]
    assert len(prep_seeds) >= args.members and len(mtppr_seeds) >= args.members
    members = []
    for k in range(1, args.members + 1):
        members.append(run_member(k, out, retain, args.band_breaks, prep_seeds[k - 1], mtppr_seeds[k - 1], args))
    os.makedirs(os.path.join(out, 'logs'), exist_ok=True)
    meta = dict(retained_transects=retain, withheld_transects=[t for t in C.STAND_TRANSECTS if t not in retain],
                band_breaks=args.band_breaks, density_mode=args.density_mode or C.DENSITY_MODE, method=args.method,
                count_control=bool(args.method == 'single' and C.COUNT_CONTROL and not args.no_count_control),
                ref_moment_rescaling=bool(args.method == 'single' and C.REF_MOMENT_RESCALING and not args.no_count_control),
                representative_criterion='lowest total final energy of copy 1 (superposition) / of the tile runs (single) over the member\'s tiles',
                drop_duplicates=args.drop_duplicates, engine_sha256=run_mtppr.engine_sha256(), members=members)
    done = [m for m in members if 'total_energy_final' in m and 'n_trees' in m]
    if done:
        rep = min(done, key=lambda m: m['total_energy_final'])
        meta['representative_member'] = rep['member']
        meta['representative_criterion'] = 'lowest total final energy over the member\'s tiles (predefined)'
        for src, dst in ((f'reconstructed_stand_m{rep["member"]}.csv', 'reconstructed_stand_stand3T.csv'),
                         (f'demo_reconstructed_m{rep["member"]}.csv', 'demo_reconstructed_stand3T.csv')):
            shutil.copyfile(os.path.join(out, src), os.path.join(out, dst))
    with open(os.path.join(out, 'stand_meta.json'), 'w') as f:
        json.dump(meta, f, indent=2, default=lambda o: int(o) if isinstance(o, (np.integer,)) else float(o))
    print(json.dumps({k: v for k, v in meta.items() if k != 'members'}, indent=2))
    for m in members:
        print({k: v for k, v in m.items() if k in ('member', 'prep_seed', 'mtppr_seed', 'n_anchors', 'total_energy_final', 'n_trees', 'by_species')})


if __name__ == '__main__':
    main()
