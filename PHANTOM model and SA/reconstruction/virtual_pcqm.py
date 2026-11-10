"""Virtual PCQM survey of a stand and comparison with the real survey records.

The survey re-applies the field protocol to a reconstructed (or reference) stand:
at every PCQM point of a transect (plot centres along its axis), the nearest tree in
each 90-degree quarter of the transect's local frame is recorded with its distance and
marks; quarters with no tree within DMAX are empty.
"""
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import ks_2samp

import recon_config as C
import recon_geometry as G


def survey_points(stand, points, dmax=C.DMAX):
    """2026-09-20 B0: PCQM survey at explicit sample points (frames-free; used for the synthetic
    stands and for plot subsets).  points: DataFrame with transect, plot, u, px, py, ux, uy, nx, ny
    (transect-local axes, as recon_geometry.pcqm_points gives them).  stand: DataFrame with xp, yp,
    species, r_stem, h_stem.  One row per quarter: the nearest tree in each 90-degree quarter,
    empty when no tree lies within dmax (dmax = np.inf for the uncapped variant)."""
    xp = stand['xp'].values.astype(float)
    yp = stand['yp'].values.astype(float)
    rows = []
    for p in points.itertuples(index=False):
        dx = xp - float(p.px)
        dy = yp - float(p.py)
        along = dx * p.ux + dy * p.uy
        perp = dx * p.nx + dy * p.ny
        dist = np.hypot(dx, dy)
        for q, (sx, sy) in C.QUARTER_SIGNS.items():
            m = (sx * perp > 0) & (sy * along > 0)
            rec = dict(transect=int(p.transect), plot=int(p.plot), quarter=q, s=float(p.u), u=float(p.u), present=False, d=np.nan,
                       d_cap=dmax, species=None, dbh_cm=np.nan, height_m=np.nan, idx=-1)
            if m.any():
                idx = np.flatnonzero(m)
                j = idx[np.argmin(dist[idx])]
                d = float(dist[j])
                if d <= dmax:
                    rec.update(present=True, d=d, d_cap=d, species=stand['species'].values[j],
                               dbh_cm=float(stand['r_stem'].values[j]) * 200.0,
                               height_m=float(stand['h_stem'].values[j]), idx=int(j))
            rows.append(rec)
    return pd.DataFrame(rows)


def survey_transect(stand, t, frames, dmax=C.DMAX):
    """stand: DataFrame with xp, yp, species, r_stem, h_stem (rotated frame). One row per quarter."""
    return survey_points(stand, G.pcqm_points({t: frames[t]}), dmax)


def cottam_curtis(d_cap):
    d_cap = np.asarray(d_cap, dtype=float)
    d_cap = d_cap[np.isfinite(d_cap)]
    if len(d_cap) == 0:
        return np.nan
    return 1.0 / np.mean(d_cap) ** 2


def jsd(p, q):
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.sum() == 0 or q.sum() == 0:
        return np.nan
    return float(jensenshannon(p / p.sum(), q / q.sum(), base=2) ** 2)


def composition(records, levels=C.SPECIES_LEVELS):
    pres = records[records['present'] == True]
    return np.array([(pres['species'] == sp).sum() for sp in levels], dtype=float)


def compare_records(real, virtual, band_by_plot):
    """Compare a real transect's PCQM records with a virtual survey of the same transect.
    band_by_plot: {plot: band} (band of the tile box holding the PCQM point)."""
    out = {}
    out['n_present_real'] = int(real['present'].sum())
    out['n_present_virtual'] = int(virtual['present'].sum())
    out['lambda_real'] = cottam_curtis(real['d_cap'])
    out['lambda_virtual'] = cottam_curtis(virtual['d_cap'])
    out['density_ratio'] = out['lambda_virtual'] / out['lambda_real'] if out['lambda_real'] else np.nan
    out['jsd_all'] = jsd(composition(real), composition(virtual))
    bands_real = np.array([band_by_plot.get(int(p), -1) for p in real['plot']])
    bands_virt = np.array([band_by_plot.get(int(p), -1) for p in virtual['plot']])
    for k in sorted(set(bands_real)):
        if k < 0:
            continue
        rr = real[bands_real == k]
        vv = virtual[bands_virt == k]
        out[f'jsd_band{k}'] = jsd(composition(rr), composition(vv))
        out[f'lambda_real_band{k}'] = cottam_curtis(rr['d_cap'])
        out[f'lambda_virtual_band{k}'] = cottam_curtis(vv['d_cap'])
    rp = real[real['present'] == True]
    vp = virtual[virtual['present'] == True]
    for col, name in (('d', 'd'), ('dbh_cm', 'dbh'), ('height_m', 'height')):
        a = rp[col].dropna().values
        b = vp[col].dropna().values
        if len(a) > 1 and len(b) > 1:
            ks = ks_2samp(a, b)
            out[f'ks_{name}_stat'] = float(ks.statistic)
            out[f'ks_{name}_p'] = float(ks.pvalue)
            out[f'mean_{name}_real'] = float(np.mean(a))
            out[f'mean_{name}_virtual'] = float(np.mean(b))
        else:
            out[f'ks_{name}_stat'] = np.nan
            out[f'ks_{name}_p'] = np.nan
    return out


def profile_by_plot(records, levels=C.SPECIES_LEVELS):
    pres = records[records['present'] == True]
    tab = pd.crosstab(pres['plot'], pres['species']).reindex(columns=levels, fill_value=0)
    return tab.reindex(range(1, int(records['plot'].max()) + 1), fill_value=0)
