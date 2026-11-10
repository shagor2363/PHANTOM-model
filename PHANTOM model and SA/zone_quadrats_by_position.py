"""Assign the observed juvenile quadrats to model zones by position (the model's own
elevation rule) instead of by vegetation class, and regenerate the zone summaries.

For each quadrat record: position from the digitised transect geometry (PCQM point + quarter
offset), elevation and zone from `tidal_dispersal_model.elevation` / `zone_at`.  The
vegetation-class zone is kept in a `veg_zone` column; `model_zone` becomes the position zone.

Usage: python zone_quadrats_by_position.py [--years 1999 1997] [--dry-run]

Rewrites (with a *_vegzone.csv backup of the originals):
    juvenile_quadrat_records_{year}_t1_5.csv
    juvenile_zone_species_summary_{year}_t1_5.csv
    juvenile_zone_totals_{year}_t1_5.csv
"""
import argparse
import csv
import importlib.util
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reconstruction'))
import recon_config as C  # noqa: E402
import transect_geometry as tg  # noqa: E402
from tidal_dispersal_model import elevation, zone_at, along_gradient  # noqa: E402

ROOT = C.ROOT


def _load_extractor(year):
    path = os.path.join(ROOT, 'extract_juvenile.py')
    spec = importlib.util.spec_from_file_location(f'extract_{year}', path)
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, [path, str(year)]     # extract_juvenile.py reads the survey year from its argument
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.argv = argv
    return mod


def quadrat_positions(df):
    xs, ys = [], []
    for r in df.itertuples(index=False):
        t, k, q = int(r.transect), int(r.plot), r.quadrat_letter
        fr = tg.transect_frame(t)
        px, py = tg.transect_position(t, float(tg.plot_u(t, k)))
        sx, sy = C.QUARTER_SIGNS[q]
        # quadrat centre: 2.5 m into the quarter along and across the transect
        xs.append(float(px) + sy * 2.5 * fr['along'][0] + sx * 2.5 * fr['perp'][0])
        ys.append(float(py) + sy * 2.5 * fr['along'][1] + sx * 2.5 * fr['perp'][1])
    return xs, ys


def rezone(year, dry_run=False):
    rec_path = os.path.join(ROOT, f'juvenile_quadrat_records_{year}_t1_5.csv')
    df = pd.read_csv(rec_path)
    if 'veg_zone' not in df.columns:
        df['veg_zone'] = df['model_zone']
    xs, ys = quadrat_positions(df)
    df['x'], df['y'] = xs, ys
    df['d_gradient'] = [along_gradient(x, y) for x, y in zip(xs, ys)]
    df['elevation'] = [elevation(x, y) for x, y in zip(xs, ys)]
    df['model_zone'] = [zone_at(x, y) for x, y in zip(xs, ys)]
    q = df.drop_duplicates('quadrat_id')
    print(f'{year}: {len(q)} quadrats; position zones {q["model_zone"].value_counts().sort_index().to_dict()}; '
          f'agreement with vegetation-class zones {(q["model_zone"] == q["veg_zone"]).mean():.2f}')
    print(pd.crosstab(q['veg_class'], q['model_zone']).to_string())
    if dry_run:
        return df
    mod = _load_extractor(year)
    records = df.to_dict('records')
    zs = mod.zone_species_summary(records)
    zt = mod.zone_totals(records)
    for path in (rec_path,
                 os.path.join(ROOT, f'juvenile_zone_species_summary_{year}_t1_5.csv'),
                 os.path.join(ROOT, f'juvenile_zone_totals_{year}_t1_5.csv')):
        bak = path.replace('.csv', '_vegzone.csv')
        if os.path.exists(path) and not os.path.exists(bak):
            shutil.copy(path, bak)
    df.to_csv(rec_path, index=False)
    pd.DataFrame(zs).to_csv(os.path.join(ROOT, f'juvenile_zone_species_summary_{year}_t1_5.csv'), index=False)
    pd.DataFrame(zt).to_csv(os.path.join(ROOT, f'juvenile_zone_totals_{year}_t1_5.csv'), index=False)
    print(f'  rewrote records + zone summaries for {year} (originals kept as *_vegzone.csv)')
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--years', nargs='+', default=['1999', '1997'])
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    for y in args.years:
        if os.path.exists(os.path.join(ROOT, f'juvenile_quadrat_records_{y}_t1_5.csv')):
            rezone(y, args.dry_run)
        else:
            print(f'{y}: no quadrat records file, skipped')


if __name__ == '__main__':
    main()
