# Parameter provenance and calibration history of `tidal_dispersal_model.py`

Written 2026-10-11 when the model script was refactored. Every dated note that stood next to
a constant in the script before that day is kept here verbatim in substance; the script now
carries only value, unit, a one-line meaning and a status tag ([lit] literature, [cal]
calibrated on the 1999 juvenile zonation, [asm] assumed). The frozen configuration of the
final v6 runs is in `frozen_configuration_2026-09-23.md`; the Sonneratia parameter audit is in
`species_parameters_sonneratia.md`; citation problems are in `citation_corrections_2026-09-22.md`.

## Domain and tidal forcing

- Model frame: the Gazi transect map of Dahdouh-Guebas et al. (2004, Fig. 6) rotated 17 deg so the
  five parallel PCQM transects run along +y (x alongshore, west to east: T4 T3 T2 T1 T6; y seaward
  to landward). The 325 x 225 m box is the tight envelope of the transects (`transect_geometry.py`).
- Elevation (Eq. 1): planar, 1.0 m at the seaward corner to 2.2 m at the landward corner along
  `GRAD_DIR`. `GRAD_ANGLE_DEG` = 15 deg is the orientation of the surveyed mangrove belt relative to
  the transect axis (seaward fringe line of T1, T2, T3, T6 tilted 18.5 deg, landward-end line 10.5 deg,
  belt centre line 15 deg, all toward +x). Zone boundaries 1.4 and 1.8 m are the thirds of the
  gradient extent D = 301 m.
- Tide (Eq. 4): H0 1.6 m, spring amplitude 1.25 m, neap 0.5 m, M2 period 12.42 h, 15-day
  spring-neap cycle (days 1-7 spring, 8-15 neap), cosine amplitude envelope peaking on day 4,
  phase offset -4.50 rad (high water at noon on day 0).

## Transport

| constant | value | history |
|---|---|---|
| `MOBILITY` | Avi 55 / Cer 45 / Rhi 40 / Son 55 m per cycle | Net landward displacement per tidal cycle, after Di Nitto et al. (2013/2014; reference to verify). Recalibrated 2026-09-06 on the map-frame domain (`mobility_experiment.py`, set A55_C45_R40, user's choice between equivalent sets); was 40/40/40 (calibrate_geom2, 2026-09-05) and 70 (Avicennia) before that. Sonneratia = Avicennia, assumed (2026-09-20); sensitivity {40, 55, 70}. |
| `SIGMA_LATERAL` | 12 / 8 / 5 / 12 m per cycle | Cross-shore spread after Van der Stocken et al. (2015; reference to verify). Sonneratia = Avicennia, assumed; sensitivity {8, 12}. |
| `DRAFT` | 0.03 / 0.12 / 0.175 / 0.02 m | Mid-points of the Table 3 literature ranges (Chen et al. 2024). Sonneratia assumed from fruit morphology (depressed-globose berry 2-5 cm, specific gravity ~0.94-1.07, Liu et al. 2025), not a measured draft; sensitivity {0.01, 0.02, 0.03}. |
| `TIDE_SCALE` | spring 1.0 / neap 0.4 | Neap/spring amplitude ratio applied to the displacement. |
| `N_CYCLES` | spring 8 / neap 4 | Transport window. Recalibrated 2026-09-05 (calibrate_geom2 `p2_dA0.20_dB0.6_bs1`); was 16 / 8. |
| `MIN_FLOAT_STEPS` | 1 / 1 | Cycles afloat before the settlement rule applies (freshly detached propagules need time to enter the flow). Recalibrated 2026-09-05; was 8 / 4. |
| `SETTLE_STEEPNESS` | 2.0 per m | Logistic steepness (Balke et al. 2013); was 10, reduced for the 1.2 m elevation range. |
| `SETTLE_TIDAL_COEFF` | 0.6 | Spring amplitude raises z_crit, delaying settlement. |
| `SETTLE_DRAFT_COEFF` | 6.0 | Recalibrated 2026-09-05; was 2.0. |
| `GROUNDING_MARGIN` | 0.02 m | ODD Eq. 7, delta_g. |
| `SUBSTEP_DT` | 2 h | Six substeps per cycle (12 h, not 12.42 h: the transport cycle is the rounded M2 period). |
| `ZONE_DRAG` | C 0.4 / B 0.7 / A 1.0 | Movement retention factor (1 = no attenuation). Mazda et al. (1997) for the concept. REVERSED 2026-09-20 at the user's request (drag strongest in the seaward zone; before: C 1.0 / B 0.6 / A 0.2), then recalibrated 2026-09-22 on the four-transect stand with T2 withheld (`recalibrate_4T.py`; was 0.2 / 0.6 / 1.0 after the reversal). The T2-blind recalibration of 2026-09-24 could not discriminate the profile and kept it. |

## Establishment

| constant | value | history |
|---|---|---|
| `R_BASE` | C 0.65 / B 0.3 / A 0.1 mm per day | Siltation at full inundation; attributed to Okello et al. (2014), a reference that does not exist as cited (see citation corrections). Zone C calibrated 2026-03-16 (`probe_zone_c_correction.py`). |
| `T_FACTOR` | spring 1.3 / neap 0.7 | ODD Eq. 8. |
| `EROSION` | 0.005 / 0.003 / 0.002 mm per day | Ebb erosion subtracted from the daily deposition. |
| `BURIAL_TOL` | Avi 30 / Cer 50 / Rhi 80 / Son 20 mm | Lethal cumulative burial. The three calibrated species' values are attributed to papers using different species or 10-18 year old trees (citation corrections). Sonneratia 20 mm: literature-based inference (reclassified 2026-09-22, was 'assumed'). Balke et al. (2013, J. Appl. Ecol. 50: 740-747) is the only sediment-dynamics experiment on S. alba; it did not measure a critical accretion threshold but infers lethal burial ~ complete burial of the shoot, their seedlings averaging 3.4 cm at day 53, and 3 x 1 cm burial over 3 weeks was largely survivable; 20 mm sits between the 2 cm start and that 3.4 cm height. Congeners are not transferable (Thampanya et al. 2002, Aquat. Bot. 74: 315-325, used the larger S. caseolaris; not the Siam Soc. paper cited in the ODD). |
| `LAMBDA_PRED` | 0.003 | Base encounter rate (Clarke & Kerrigan 2002); calibrated 2026-02-28. |
| `PRED_D` | C 24.0 / B 7.97 / A 4.24 | Effective predator density; zone C calibrated 2026-03-16. |
| `VULN` | 0.85 / 0.55 / 0.25 / 0.85 | Species vulnerability; the ordering is well supported, the numbers are not measured. Sonneratia = Avicennia, assumed; sensitivity {0.55, 0.85}. |
| `TIDE_E` | spring 0.8 / neap 1.2 | Tidal modifier of predation. |
| `AGE_A` | days 1-14: 1.5, 15-30: 1.0, 31-90: 0.7 | Age modifier of predation. |
| `H_REF` | 10 / 28 / 30 / 8 cm | Height at which the height protection term reaches its floor; Sonneratia assumed. |
| `MAX_DAILY_P` | 0.15 | Cap on the daily predation probability. |
| `INIT_HEIGHT` | 1.4 / 8.25 / 18.0 / 2.0 cm | Sonneratia: measured nursery height 30 days after sowing (Trangia 2022), used as the starting height of a newly rooted seedling. |
| `DAILY_GROWTH` | 0.3 / 0.22 / 0.25 / 0.05 cm per day | Sonneratia 0.05: two converging studies. Trangia (2022, IJERD 13: 184-189) measured 1.87-5.41 cm gain over 90 d from a 2 cm transplant = 0.021-0.060 cm per day (soil-pooled 0.045); Balke et al. (2013) reached 3.4 cm at day 53 from germination = 0.05-0.06 cm per day. A change to 0.0264 on 2026-09-22 wrongly treated Balke's 3.4 cm as growth from the 2 cm start and was reverted the same night. Nursery and mesocosm only; no field growth rate exists. |
| `SUBSTRATE_SUITABILITY` | forest 1.0 / sand 0.0 | Eq. 10, simplified substrate rule of 2026-09-06 (replaces the five-class map of 2026-09-05): the model frame has two substrates, forest floor and the treeless sand band seaward of the forest edge (`results/reconstruction/substrate_map.npz`, made by `reconstruction/make_substrate_map.py`; the band is everything seaward of the line through the transect seaward ends, drawn by hand west of T1 so that it narrows towards T4 as the sand bank does on the map). On the sand band a stranded propagule cannot root and dies at day 0 with cause 'substrate', which is excluded from the seedling-mortality statistics. The same band is the only mask on the reconstructed stand. |

## Propagule supply

| constant | value | history |
|---|---|---|
| `F_BASE` | Avi 12.8 / Cer 48.0 / Rhi 8.0 / Son 0.00384 seeds per year at r_min | Annual seeds of a minimum-size tree (ODD Eq. 2; Clarke 1992 for the concept). Rebalanced 2026-03-16: internal halved, boundary doubled, total seed budget preserved. Sonneratia: 12.8 x s with s = 3.0e-04, FITTED 2026-09-22 on the mean recorded 1999 juvenile density over the 309 in-domain quadrats (`son_refit.py`, 5 seeds x 4 levels, realised curve under the stochastic-rounding scheduler). At this level the model gives 6.21e-04 per m2 against the target 6.15e-04. It is the lowest tested level, so a bound rather than a bracketed optimum: the observational bootstrap (all 309 quadrats, zeros retained) spans 0 to 1.49e-03 with 5.1 % zero-detection draws, giving a multiplier interval of ~3.0e-04 to 7.8e-04 with no meaningful lower bound. An effective supply scale conditional on the fixed survival rules, not a fecundity estimate. |
| `GAMMA_FECUNDITY` | 2.0 / 2.5 / 2.5 / 2.0 | Allometric exponent (ODD Eq. 2); Sonneratia = Avicennia, assumed. |
| `MATURITY` | 0.027 / 0.027 / 0.024 / 0.025 m r_stem | Minimum reproductive stem radius. Sonneratia 0.025 m (5 cm DBH): literature upper bound, site-specific. Wang'ondu PhD thesis (Univ. Nairobi, Gazi Bay), Table 3.1: the 11-year reforested stand, mean DBH 5.09 +- 0.82 cm, produced flowers, fruits and propagule fall throughout 2005-2006. No published minimum reproductive size exists, so this is an upper bound, not a threshold. Was 0.035 m (7 cm DBH, assumed), which exceeded a stand demonstrably in fruit. Changed 2026-09-22; note F = F_BASE (r / r_min)^gamma, so r_min moves both eligibility and per-tree output. |
| `BOUNDARY_SPRING`, `BOUNDARY_NEAP` | 448 / 384 / 160 / 0.1344 and 224 / 192 / 64 / 0.0672 propagules per day | Peak-season influx at the seaward edge, modulated by `seasonal_scale()`. Rebalanced 2026-03-16 (doubled) and multiplied by 4 on 2026-09-05 (calibrate_geom2 `p4_bs4`: boundary propagules ~31 % of the cohort; the boundary scale is what separates the combined scenario from internal-only). The 2026-09-22 recalibration grid (`recalibrate_4T.py`) tested a x2 boundary scale; the frozen configuration of 2026-09-23 carries the single supply 448 / 384 / 160 listed here. Sonneratia: 448 x s and 224 x s with the fitted s = 3.0e-04 (same fit as F_BASE). |
| `FRUITING_WINDOW` | see script | Gazi Bay phenology: Wang'ondu et al. (2013), Ochieng & Erftemeijer (2002), Slim et al. (1996). Sonneratia (213, 244, 319, 349): literature, site-specific (placeholder flag removed 2026-09-22). Wang'ondu et al. (2013, Flora 208: 522-531) and the underlying PhD thesis, Gazi Bay, three stands, fortnightly 2005-2006: flower bud primordia Jul-Sep, fruit formation Aug-Nov, peak propagule fall September-November, cycle 4-5 months. The Oct 2005 / Sep 2006 peaks often quoted are fruit on the tree, not fall; a narrowing to an Oct 31 peak end on 2026-09-22 was made on that misreading and reverted the same night. The small Jan-Feb mature-fruit tail the thesis reports is not represented (stated sensitivity). |

## Stochastic release counting (2026-09-22)

Boundary releases per species and day are `stochastic_round(rate x seasonal weight)`, i.e.
floor(lambda) + Bernoulli(lambda - floor(lambda)), which is unbiased (E[N] = lambda). Before
2026-09-22 the count was `int(rate x season)`, which truncates: any expected release below one
propagule became zero, discarding the whole budget. For Sonneratia at its fitted supply that
silently removed ~30 propagules per year (28 % of its internal source) while leaving the
internal source scaled exactly, so the supply multiplier did not scale both sources
proportionally. The rounding is applied once, after every multiplier, to the boundary
releases of all species; the internal Poisson fecundity process is unchanged. It alters the
release-counting rule only; the supply, transport and establishment equations are untouched.

## Release-count generators and reproducibility (fixed 2026-10-11)

The rounding draws are taken from one generator per species, separate from the main
generator, so that adding or changing them cannot shift the transport and mortality draws.
Until 2026-10-11 these generators were seeded with
`abs(hash((seed, 'release', species))) % 2**32`. Python salts the hash of strings per
interpreter process (PEP 456), so the seed differed between runs: the Bernoulli part of every
fractional release count (all shoulder-season days of every species, and every Sonneratia
day) was not reproducible from the run seed. Verified 2026-10-11: two 75-day runs of the
pre-refactor script with seed 42 in separate processes gave 72,632 and 72,636 result rows,
and `PYTHONHASHSEED` is not set on the machine that produced the v6 caches. Consequence: the v6 caches are one valid
realisation of the model but cannot be regenerated bit for bit from their seeds; all
aggregate results (20 seeds, replicate ranges) are unaffected in kind. The generators are
now spawned from `numpy.random.SeedSequence(seed)`, which is deterministic everywhere; the
main generator `default_rng(seed)` is unchanged, so runs without fractional releases are
identical to the old code. Verified 2026-10-11: two 75-day runs of the refactored script in
separate processes are byte-identical (72,629 rows), and with the old hash seeding re-attached
under a pinned `PYTHONHASHSEED` the refactored script reproduces the old script byte for byte
(72,631 rows), so the seeding is the only behavioural change of the refactor.
