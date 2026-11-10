"""Stage 4: assemble MTPPR tile outputs into the reconstructed stand (model frame).

Usage:
    python reconstruction/postprocess.py --geom results/reconstruction/geom \
        --mtppr results/reconstruction/mtppr --out results/reconstruction [--all-seeds] [--extend-landward]
        [--clip-fringe] [--fringe-surveyed-start] [--clip-treeless] [--clip-substrate]

Outputs in --out:
    reconstructed_stand.csv          best seed per tile; demo.csv columns + provenance columns
    demo_reconstructed.csv           exact demo.csv schema (plant,x,y,r_stem,h_stem,species,r_crown,r_root)
    reconstructed_stand_seed{s}.csv  one replicate stand per seed (with --all-seeds)
    stand_meta.json                  counts and clipping / seam statistics
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recon_config as C  # noqa: E402
import recon_geometry as G  # noqa: E402

DEMO_COLS = ['plant', 'x', 'y', 'r_stem', 'h_stem', 'species', 'r_crown', 'r_root']


def load_geom(geom_dir):
    g = dict(
        tiles=pd.read_csv(os.path.join(geom_dir, 'tiles.csv')),
        core_rects=pd.read_csv(os.path.join(geom_dir, 'core_rects.csv')),
        ref=pd.read_csv(os.path.join(geom_dir, 'reference.csv')),
        anchors=pd.read_csv(os.path.join(geom_dir, 'anchors.csv')),
    )
    with open(os.path.join(geom_dir, 'prep_meta.json')) as f:
        meta = json.load(f)
    with open(os.path.join(geom_dir, 'allometry.json')) as f:
        g['allometry'] = json.load(f)
    g['segments'] = {int(t): [tuple(s) for s in v] for t, v in meta['tile_segments'].items()}
    g['belt_segments'] = {int(t): [tuple(s) for s in v] for t, v in meta['belt_segments'].items()}
    g['frames'] = {int(t): v for t, v in meta['frames'].items()}
    g['corridors'] = {int(t): tuple(v) for t, v in meta['corridors'].items()}
    g['meta'] = meta
    return g


def load_tile_outputs(mtppr_dir, tiles, selection=None, seed=None):
    parts = []
    for r in tiles.itertuples(index=False):
        s = seed if seed is not None else selection[str(int(r.tile_id))]
        path = os.path.join(mtppr_dir, f'tile{int(r.tile_id)}_seed{int(s)}_points.csv')
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        p = pd.read_csv(path)
        if 'in_tile' in p.columns:
            p = p[p['in_tile'] == 1].drop(columns=['in_tile'])
        p['tile_id'] = int(r.tile_id)
        p['transect'] = int(r.transect)
        p['band'] = int(r.band)
        p['seed'] = int(s)
        parts.append(p)
    return pd.concat(parts, ignore_index=True)


def tag_sources(df, ref, tiles):
    """Attach reference provenance to fixed rows; everything else is MTPPR-generated.
    Reference points outside every tile (e.g. T5's tails) are appended as fixed rows."""
    key = pd.MultiIndex.from_arrays([ref['xp'].round(6), ref['yp'].round(6)])
    lookup = pd.DataFrame(dict(source=ref['source'].values, anchor_id=ref['anchor_id'].values,
                               ref_id=ref['ref_id'].values, r_crown_ref=ref['r_crown'].values,
                               ref_tile=ref['tile_id'].values, ref_transect=ref['transect'].values), index=key)
    dkey = pd.MultiIndex.from_arrays([df['xp'].round(6), df['yp'].round(6)])
    hit = lookup.reindex(dkey)
    df = df.copy()
    fixed = df['is_fixed'].values == 1
    df['source'] = np.where(fixed, hit['source'].values, 'mtppr')
    df['anchor_id'] = np.where(fixed, hit['anchor_id'].values, -1)
    df['ref_id'] = np.where(fixed, hit['ref_id'].values, -1)
    df['r_crown_ref'] = np.where(fixed, hit['r_crown_ref'].values, np.nan)
    df['ref_transect'] = np.where(fixed, hit['ref_transect'].values, -1)
    unmatched = int(pd.isna(hit['source'].values[fixed]).sum())
    assert unmatched == 0, f'{unmatched} fixed points not found in the reference'
    for tid, gg in df[fixed].groupby('tile_id'):
        n_ref = int((ref['tile_id'] == tid).sum())
        assert len(gg) == n_ref, f'tile {tid}: {len(gg)} fixed rows vs {n_ref} reference rows in the box'
    orphans = ref[ref['tile_id'] < 0]
    if len(orphans):
        o = pd.DataFrame(dict(xp=orphans['xp'].values, yp=orphans['yp'].values,
                              r_stem_s=orphans['r_stem_s'].values, h_stem_s=orphans['h_stem_s'].values,
                              species=orphans['species'].values, is_fixed=1, tile_id=-1,
                              transect=orphans['transect'].values, band=orphans['band'].values,
                              seed=df['seed'].iloc[0] if len(df) else 0, source=orphans['source'].values,
                              anchor_id=orphans['anchor_id'].values, ref_id=orphans['ref_id'].values,
                              r_crown_ref=orphans['r_crown'].values, ref_transect=orphans['transect'].values))
        df = pd.concat([df, o], ignore_index=True)
    return df


def clip_to_envelope(df, geom, extend_landward=False, clip_envelope=False):
    frames, corridors, segments = geom['frames'], geom['corridors'], geom['segments']
    xp, yp = df['xp'].values, df['yp'].values
    x, y = G.to_model(xp, yp)
    par = {t: f for t, f in frames.items() if f['parallel']}
    fringe, land = G.envelope(par, xp, segments)
    if extend_landward or not clip_envelope:
        land = np.full_like(land, np.inf)
    if not clip_envelope:
        fringe = np.full_like(fringe, -np.inf)
    u, corr = G.fringe_distance(xp, yp, par, corridors)
    df = df.copy()
    df['x'], df['y'], df['u'], df['corridor'] = x, y, u, corr
    inside_dom = G.in_domain(x, y)
    inside_env = (yp >= fringe) & (yp < land)
    keep = inside_dom & (inside_env | (df['is_fixed'].values == 1))
    stats = dict(n_before=int(len(df)), n_outside_domain=int((~inside_dom).sum()),
                 n_outside_envelope=int((inside_dom & ~inside_env & (df['is_fixed'].values == 0)).sum()),
                 n_fixed_outside_domain=int((~inside_dom & (df['is_fixed'].values == 1)).sum()))
    return df[keep].reset_index(drop=True), stats


def fringe_line(frames, xp, surveyed_start=False):
    """Seaward forest edge yp at cross positions xp: piecewise-linear across x through
    the seaward ends (x0, y0) of the parallel transects, constant beyond the outer transects.
    With surveyed_start=False T4's end is the extrapolated tape start (its unsurveyed 0-120 m
    stays inside the stand); with surveyed_start=True it is the start of T4's surveyed part
    (y0 + tape_offset, the only evidence for T4), the line used since 2026-09-06."""
    order = G.parallel_ids(frames)
    xs = np.array([frames[t]['x0'] for t in order], dtype=float)
    ys = np.array([frames[t]['y0'] + (frames[t]['tape_offset'] if surveyed_start else 0.0) for t in order],
                  dtype=float)
    return np.interp(np.asarray(xp, dtype=float), xs, ys)


def clip_to_fringe(df, geom, surveyed_start=False):
    """Drop synthetic trees seaward of the mapped forest fringe (sand / water on the
    Dahdouh-Guebas et al. 2004 map); recorded (fixed) trees are always kept."""
    frames = geom['frames']
    fr = fringe_line(frames, df['xp'].values, surveyed_start)
    seaward = (df['yp'].values < fr) & (df['is_fixed'].values == 0)
    stats = dict(n_seaward_of_fringe=int(seaward.sum()),
                 seaward_by_species=df.loc[seaward, 'species'].value_counts().to_dict())
    return df[~seaward].reset_index(drop=True), stats


def clip_to_treeless(df):
    """Drop synthetic trees inside the treeless band of the substrate map (seaward of the forest
    edge through the transect seaward ends, following the mapped sand-bank boundary west of T3;
    'treeless' grid of results/reconstruction/substrate_map.npz); recorded (fixed) trees kept."""
    from tidal_dispersal_model import SUBSTRATE_MAP_FILE
    d = np.load(SUBSTRATE_MAP_FILE)
    mask, cell = d['treeless'], float(d['cell'])
    mx, my = G.to_model(df['xp'].values, df['yp'].values)
    iy = np.clip(((my - float(d['y_min'])) / cell).astype(int), 0, mask.shape[0] - 1)
    ix = np.clip(((mx - float(d['x_min'])) / cell).astype(int), 0, mask.shape[1] - 1)
    inside = mask[iy, ix]
    drop = inside & (df['is_fixed'].values == 0)
    stats = dict(n_in_treeless_band=int(drop.sum()),
                 treeless_by_species=df.loc[drop, 'species'].value_counts().to_dict(),
                 fixed_in_treeless_band=int((inside & (df['is_fixed'].values == 1)).sum()))
    return df[~drop].reset_index(drop=True), stats


def clip_to_substrate(df):
    """Drop synthetic trees standing on cells that the substrate map of the model frame
    (results/reconstruction/substrate_map.npz, Eq. 10) classes as sand bank, water / seagrass
    or terrestrial, so that the parent stand and the rooting rule share one physiography;
    recorded (fixed) trees are always kept."""
    from collections import Counter
    from tidal_dispersal_model import substrate_at, SUBSTRATE_SUITABILITY
    mx, my = G.to_model(df['xp'].values, df['yp'].values)
    sub = np.array([substrate_at(x, y) for x, y in zip(mx, my)])
    unsuitable = np.array([SUBSTRATE_SUITABILITY.get(c, 1.0) < 1.0 for c in sub], dtype=bool)
    drop = unsuitable & (df['is_fixed'].values == 0)
    stats = dict(n_on_unsuitable_substrate=int(drop.sum()),
                 unsuitable_by_class=dict(Counter(sub[drop])),
                 unsuitable_by_species=df.loc[drop, 'species'].value_counts().to_dict(),
                 fixed_on_unsuitable_substrate=int((unsuitable & (df['is_fixed'].values == 1)).sum()))
    return df[~drop].reset_index(drop=True), stats


def seam_cleanup(df, hardcore=C.HARDCORE):
    """Remove non-fixed points closer than the hard core to a point of another tile."""
    xy = df[['xp', 'yp']].values
    tree = cKDTree(xy)
    pairs = np.array(list(tree.query_pairs(hardcore)), dtype=int)
    drop = set()
    n_cross = 0
    if len(pairs):
        tid = df['tile_id'].values
        fx = df['is_fixed'].values == 1
        for i, j in pairs:
            if tid[i] == tid[j]:
                continue
            n_cross += 1
            if i in drop or j in drop:
                continue
            if fx[i] and fx[j]:
                continue
            if fx[i]:
                drop.add(j)
            elif fx[j]:
                drop.add(i)
            else:
                drop.add(max(i, j))
    keep = np.ones(len(df), dtype=bool)
    keep[list(drop)] = False
    return df[keep].reset_index(drop=True), dict(n_cross_tile_pairs=n_cross, n_seam_removed=len(drop))


def add_marks(df, allometry):
    df = df.copy()
    df['r_stem'] = df['r_stem_s'] * C.R_STEM_SCALE
    df['h_stem'] = df['h_stem_s'] * C.H_STEM_SCALE
    dbh = df['r_stem'] * 200.0
    crown = np.empty(len(df))
    for i, (sp, d) in enumerate(zip(df['species'], dbh)):
        cf = allometry['crown'].get(sp) or allometry['crown']['pooled']
        crown[i] = max(0.3, cf['c0'] + cf['c1'] * d)
    df['r_crown'] = np.where(df['r_crown_ref'].notna(), df['r_crown_ref'], crown)
    root = np.array([C.ROOT_ALLOM[sp][0] * d ** C.ROOT_ALLOM[sp][1] for sp, d in zip(df['species'], dbh)])
    df['r_root'] = np.round(root, 3)
    return df


def assign_ids(df):
    df = df.sort_values(['species', 'y', 'x']).reset_index(drop=True)
    ids = []
    counters = {}
    for sp in df['species']:
        counters[sp] = counters.get(sp, 0) + 1
        ids.append(f'{C.PLANT_PREFIX[sp]}_{counters[sp]:09d}')
    df.insert(0, 'plant', ids)
    return df


def assemble(geom, mtppr_dir, selection=None, seed=None, extend_landward=False, clip_envelope=False,
             clip_fringe=False, clip_substrate=False, fringe_surveyed_start=False, clip_treeless=False):
    df = load_tile_outputs(mtppr_dir, geom['tiles'], selection, seed)
    df = tag_sources(df, geom['ref'], geom['tiles'])
    df, clip_stats = clip_to_envelope(df, geom, extend_landward, clip_envelope)
    if clip_fringe:
        df, fringe_stats = clip_to_fringe(df, geom, fringe_surveyed_start)
        clip_stats.update(fringe_stats)
    if clip_treeless:
        df, treeless_stats = clip_to_treeless(df)
        clip_stats.update(treeless_stats)
    if clip_substrate:
        df, sub_stats = clip_to_substrate(df)
        clip_stats.update(sub_stats)
    df, seam_stats = seam_cleanup(df)
    df = add_marks(df, geom['allometry'])
    df = assign_ids(df)
    extra = ['source', 'transect', 'ref_transect', 'band', 'corridor', 'tile_id', 'seed', 'anchor_id', 'xp', 'yp', 'u', 'is_fixed']
    df = df[DEMO_COLS + extra]
    for c in ('x', 'y', 'xp', 'yp', 'u'):
        df[c] = df[c].round(4)
    for c in ('r_stem', 'h_stem', 'r_crown'):
        df[c] = df[c].round(4)
    meta = dict(n_total=int(len(df)), by_source=df['source'].value_counts().to_dict(),
                by_species=df['species'].value_counts().to_dict(),
                by_band=df.groupby('band').size().to_dict(), **clip_stats, **seam_stats)
    return df, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--geom', default=C.GEOM_DIR)
    ap.add_argument('--mtppr', default=C.MTPPR_DIR)
    ap.add_argument('--out', default=C.OUT_ROOT)
    ap.add_argument('--seeds', type=int, nargs='*', default=None, help='replicate stands to write (default: all in selection)')
    ap.add_argument('--all-seeds', action='store_true')
    ap.add_argument('--extend-landward', action='store_true')
    ap.add_argument('--clip-envelope', action='store_true', help='drop synthetic trees outside the surveyed fringe/landward envelope (default: fill the whole domain)')
    ap.add_argument('--clip-fringe', action='store_true',
                    help='drop synthetic trees seaward of the mapped forest fringe (line through the transect seaward ends); landward fill kept')
    ap.add_argument('--fringe-surveyed-start', action='store_true',
                    help='with --clip-fringe: take T4 at the start of its surveyed part (tape 120) instead of the extrapolated tape start')
    ap.add_argument('--clip-treeless', action='store_true',
                    help="drop synthetic trees in the 'treeless' band of the substrate map (seaward of the forest edge)")
    ap.add_argument('--clip-substrate', action='store_true',
                    help='drop synthetic trees on sand bank / water / terrestrial cells of the substrate map (Eq. 10)')
    ap.add_argument('--tag', default='')
    args = ap.parse_args()
    geom = load_geom(args.geom)
    with open(os.path.join(args.mtppr, 'selection.json')) as f:
        selection = json.load(f)
    os.makedirs(args.out, exist_ok=True)
    tag = args.tag
    df, meta = assemble(geom, args.mtppr, selection=selection, extend_landward=args.extend_landward, clip_envelope=args.clip_envelope,
                        clip_fringe=args.clip_fringe, clip_substrate=args.clip_substrate, fringe_surveyed_start=args.fringe_surveyed_start, clip_treeless=args.clip_treeless)
    df.to_csv(os.path.join(args.out, f'reconstructed_stand{tag}.csv'), index=False)
    df[DEMO_COLS].to_csv(os.path.join(args.out, f'demo_reconstructed{tag}.csv'), index=False)
    from tidal_dispersal_model import F_BASE
    df[df['species'].isin(F_BASE)][DEMO_COLS].to_csv(os.path.join(args.out, f'demo_reconstructed{tag}_model_species.csv'), index=False)
    meta_all = dict(best=meta, selection=selection, extend_landward=args.extend_landward, clip_envelope=args.clip_envelope,
                    clip_fringe=args.clip_fringe, clip_substrate=args.clip_substrate, fringe_surveyed_start=args.fringe_surveyed_start, clip_treeless=args.clip_treeless)
    seeds = args.seeds if args.seeds else (list(range(1, C.N_SEEDS + 1)) if args.all_seeds else [])
    for s in seeds:
        try:
            ds, ms = assemble(geom, args.mtppr, seed=s, extend_landward=args.extend_landward, clip_envelope=args.clip_envelope,
                              clip_fringe=args.clip_fringe, clip_substrate=args.clip_substrate, fringe_surveyed_start=args.fringe_surveyed_start, clip_treeless=args.clip_treeless)
        except FileNotFoundError as e:
            print('skip seed', s, e)
            continue
        ds.to_csv(os.path.join(args.out, f'reconstructed_stand{tag}_seed{s}.csv'), index=False)
        ds[DEMO_COLS].to_csv(os.path.join(args.out, f'demo_reconstructed{tag}_seed{s}.csv'), index=False)
        meta_all[f'seed{s}'] = ms
    with open(os.path.join(args.out, f'stand_meta{tag}.json'), 'w') as f:
        json.dump(meta_all, f, indent=2, default=int)
    print(json.dumps(meta, indent=2, default=int))


if __name__ == '__main__':
    main()
