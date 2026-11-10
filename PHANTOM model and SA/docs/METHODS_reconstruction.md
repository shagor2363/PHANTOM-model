# Methods note: whole-domain reconstruction of the adult stand at Gazi Bay

Generated 2026-09-04 from `GAZ1_6 (1).XLS` (sheet `GAZ1 env.summary`) and the transect map of
Dahdouh-Guebas et al. (2004, Fig. 6a) with the scripts in `reconstruction/`. The stand used by
the model is in `results/reconstruction/` (three modelled species, recorded-tree density). Two
earlier variants are archived for sensitivity: `results/reconstruction_5sp/` (all five recorded
species, recorded-tree density) and `results/reconstruction_pcqm/` (all five species, belts
filled to the survey's distance-based density). All numbers are read from those folders
(`geom/`, `tables/`, `mtppr/mtppr_summaries.csv`, `stand_meta.json`).

## 1. Field data and what they can support

The adult inventory is a point-centred quarter method (PCQM) survey (Cottam & Curtis 1956;
Cintron & Schaeffer Novelli 1984) along five parallel land-sea transects (T4, T3, T2, T1, T6 from
west to east) and one orthogonal transect (T5) that the authors "covered in order to corroborate
the zonation pattern" and do not show in their Fig. 6 (Dahdouh-Guebas et al. 2004, p. 80). Only
the five parallel transects are used. At every 10 m of tape the nearest adult tree in each of
four quarters was recorded with its distance `d` to the sampling point, girth/DBH (d130), height
and crown diameter: 68 sampling points, 272 quarters and 266 adults (Rhizophora mucronata 142,
Sonneratia alba 56, Ceriops tagal 34, Avicennia marina 30, Bruguiera gymnorhiza 4).

**Species.** The stand is reconstructed for the four species the model simulates: R. mucronata,
A. marina, C. tagal and, since 2026-09-20, S. alba. S. alba dominates the seaward fringe (56
recorded adults on the five parallel transects, 25 more on T5) and recruits only marginally inside
the forest (11 juvenile quadrat records in 1997, 3 in 1999, i.e. 0.6 % and 0.02 % of the juvenile
density; the earlier statement that it has no recruits was wrong). B. gymnorhiza (4 adults, 5
juvenile records per year) is not placed. Their recorded distances are still used as exclusion sectors: a
recorded Sonneratia at 3 m means that no tree of any species, and therefore no reconstructed
tree, stands closer than 3 m in that quarter.

Two further properties of the records determine the design:

1. **Orientation.** Along every transect the vegetation class recorded at the sampling points
   runs `A. mar` -> `Mixed` -> `R. muc` -> `S. alb` with increasing tape distance, the last plots
   carry the remarks "water level at 105 cm depth", "much sand" and "seagrasses on mangrove
   floor", and the GPS fixes give bearings of 108-128 deg, i.e. towards the bay. Tape 0 is the
   landward edge and the last plot the seaward Sonneratia fringe. T4's six plots occupy tape
   120-180 m of a 180 m line (its tape column reads 180 -> 120 and plot 1 is marked
   "INTERSECTION 5.1"), so T4 was sampled only in its landward part.
2. **Density.** PCQM records are a distance sample, not a census. Counting the recorded trees
   per sampled strip gives about 400 adults/ha (all species), whereas the survey's own
   distance-based estimator (Cottam-Curtis, lambda = 1/mean(d)^2, distances capped at 20 m)
   gives 340-2,450 adults/ha along the gradient (mean nearest-tree distance 5.5 m at the fringe,
   2.0 m in the landward band). The stand uses the recorded-tree reading, the same convention as
   the Pambala reconstruction; the distance-based reading is the archived sensitivity variant
   and its consequence is quantified in section 10.

## 2. Model frame and transect layout

The figure was rendered at 300 dpi and digitised at 4.0 px per metre (one PCQM block of 10 m =
40 px; the 300 m close-up = 1,216 px). The seaward end of each drawn strip was read off the
figure; the strips lean 17 deg on the printed map and lie 70-78 m apart. The model frame is the
map mirrored (sea at the bottom) and rotated by 17 deg so that the transects run exactly along
+y, the land-sea gradient they were laid on: x is alongshore, y is seaward -> landward. The
transects sit at x = 14 (T4), 92 (T3), 166 (T2), 235 (T1) and 308 m (T6); their seaward
fringes at y = 42 (T4, extrapolated from its landward part), 82, 38, 5 and 13 m; their landward
ends at y = 223, 212, 199, 175 and 173 m. Plot k of a transect sits at fringe distance
`u = tape_length - 10 (k - 1) - 5` from the seaward end of its tape.

The domain is the tight envelope of the five transects: 325 x 225 m (x 0..325, y 0..225). Bed
elevation is planar (Eq. 1), rising from 1.0 m at the seaward corner (0, 0) to 2.2 m at the
landward corner (325, 225) along a gradient tilted 15 deg from the transect axis toward +x:
the survey's own belt is oblique to the transects (seaward fringe line of T1, T2, T3 and T6
tilted 18.5 deg, landward-end line 10.5 deg, belt centre 15 deg), so the iso-elevation lines
and the three model zones are diagonal bands 100 m wide along the gradient (zone boundaries
from (0, 104) to (325, 17) m and from (0, 208) to (325, 121) m). Zone C is the seaward band
(mudflat and Sonneratia fringe), B the Rhizophora belt, A the landward Ceriops/Avicennia belt. `transect_geometry.py` holds the digitised geometry;
`figures/fig_transect_layout.png` shows the paper's map next to the model frame with the PCQM
quadrats coloured by the recorded species.

The stand is clipped at the mapped forest fringe (added 2026-09-05): synthetic trees seaward of
the piecewise-linear line through the transects' seaward ends (T4 14/42 m, extrapolated; T3
92/82; T2 166/38; T1 235/5; T6 308/13; constant beyond T4 and T6) are removed, because the map
shows sand, sand bank and water with seagrass there (`postprocess.py --clip-fringe`). This removed
167 of 2,276 synthetic trees (125 Rhizophora, 42 Avicennia), all in zone C; the landward fill and
the recorded (anchor) trees are unchanged, leaving 2,109 trees. The unclipped stand is archived as
`demo_reconstructed_filled.csv` (with `reconstructed_stand_filled.csv` and the replicate seeds
`*_filled_seed*.csv`).

The same map provides the substrate raster used by the model's establishment stage (added 2026-09-05,
`make_substrate_map.py`): every pixel of the 300-dpi page render inside the map frame is classified by
the nearest legend colour, transformed into the model frame and binned into 1 m cells by majority vote
(white pixels, i.e. quadrat symbols and the seaward margin, are filled from the nearest classified cell);
classes are grouped into mangrove, sandy stretch, sand bank, water + seagrass and terrestrial. Inside the
frame: zone C 47 % mangrove / 49 % sand bank / 3 % water, zone B 84 % mangrove / 11 % sand bank, zone A
89 % mangrove / 10 % sandy stretch (`figures/fig_substrate_map.png`).

## 3. Stratification and tiles

- **Corridors** are bounded by the midlines between neighbouring transects (x = 53, 129, 200,
  272 m); the outer corridors extend to the domain edges.
- **Bands** are classes of distance from each corridor transect's seaward fringe: 0-40, 40-80,
  80-120 and >= 120 m (a remnant shorter than 20 m is merged into the previous band: T3's belt
  ends at 130 m in band 80-130). Bands are pooled across transects. The pooled adult composition
  supports the breaks: the fringe band is dominated by Sonneratia and holds few modelled trees,
  Rhizophora dominates 40-120 m, Ceriops and Avicennia from 120 m.
- **Tiles** (19) are the corridor x band rectangles. The first band of every corridor extends
  to the seaward domain edge and the last band to the landward edge, so the tiles cover the
  whole 325 x 225 m box; each is reconstructed separately from the pooled reference of its band,
  and the band of any point is the band of the tile that contains it. Each tile is
  reconstructed in a window buffered by 8 m on all sides and clipped to the tile afterwards, so
  that edge effects of the optimiser fall into the discarded margin.
- **Reference belts** are strips along each transect over its sampled length (16 belt
  segments), each as wide as the recorded trees it holds (section 5).
- **Extrapolated area.** The stand fills the box, so the parts of the domain seaward of a
  transect's fringe, landward of its end, and the whole of T4's corridor seaward of its surveyed
  part (8,600 m2) are extrapolations of the neighbouring band's reference; they carry the
  `source = mtppr` label and the surveyed extent is drawn on the maps.

## 4. Stage 1: positioning the recorded trees (anchors)

Each recorded adult of the three species is placed at exactly its recorded distance `d` from
its sampling point, at a uniformly random angle inside its 90 deg quarter (a = seaward-left,
b = seaward-right, c = landward-left, d = landward-right; maximum-entropy choice of the unknown
angle, as in the Pambala reconstruction), constrained to stay within 10 m of the transect line.
An angle is re-drawn (up to 200 times) if the tree would fall inside the exclusion sector of
another sampling point; no unresolvable conflicts occurred. Excluded from placement: four
records without a distance (T3 plot 13) and one record with `d` > 20 m. Result: 203 anchors
(Rhizophora 139, Ceriops 34, Avicennia 30).

**Exclusion sectors.** A recorded nearest tree at distance `d` implies that no other adult
stands closer than `d` within that quarter; the sectors of all 260 recorded trees with a usable
distance (including the 57 Sonneratia and Bruguiera) are enforced on all synthetic trees.
Quarters without a recorded tree within 20 m get no sector, because the recorded trees of
neighbouring points fall inside such a zone.

## 5. Stage 2: reference density

The reference of a band is the set of recorded modelled-species trees in the band's belt
segments, each segment as wide as the trees it holds (half width = largest perpendicular offset
of its anchors + 0.25 m, floored at 2.5 m and capped at 10 m; 2.5-6.2 m here). No trees are
added, so the reference density is the count of recorded trees per belt area:

| band (m from fringe) | recorded trees | belt area (m2) | trees/ha |
|---|---|---|---|
| 0-40 | 17 | 1,468 | 116 |
| 40-80 | 48 | 1,386 | 346 |
| 80-120 | 65 | 1,658 | 392 |
| >= 120 | 73 | 1,316 | 555 |

This is the recipe of the Pambala reconstruction (~450 trees/ha there). For comparison, the
Cottam-Curtis species density of the same records (all-species distance density x modelled
share) is 108, 584, 1,835 and 2,356 trees/ha; a virtual PCQM survey of the belts reproduces all
203 recorded distances exactly.

## 6. Stage 3: MTPPR extrapolation to the unsampled area

Each tile is reconstructed with Multi-trait Point Pattern Reconstruction v2.0.0 (Wudel et al.
2023; `MTPPR-v.2.0.0.R`, used unmodified). Settings:

- Reference: all belt trees of the tile's band (union of the band's belt segments, supplied
  through the engine's `in_window` hook; Ohser-Stoyan translation edge correction evaluated for
  the union of belts and for the buffered tile rectangle).
- Marks: stem radius (r_stem / 0.4), stem height (h_stem / 17) and species (factor); the
  engine's full set of pair-correlation, mark-correlation and mark-product statistics with unit
  weights, Epanechnikov kernel (bandwidth 0.5 m), r = 0.25-8 m in 0.25 m steps, L1 energy.
- Conditioning: the belt trees inside the buffered tile are fixed (copied verbatim; the engine
  never proposes points there); all exclusion sectors act as forbidden zones.
- Optimisation: greedy acceptance with move, swap, mark-exchange, mark-draw, add and delete
  proposals (engine defaults, so the tree count adapts to the reference intensity); stop at
  max(40,000, 200 n) steps, at 400,000 steps, or when the energy improves by less than 2e-6 of its
  value per step over 5,000 steps. Five seeds per tile; the seed with the lowest final energy is
  used, all five are kept as replicate stands.
- 95 jobs (19 tiles x 5 seeds); all converged or reached their budget; final energies
  0.001-0.024 (initial 0.1-0.5).

## 7. Stage 4: assembly

Tile outputs are clipped to their tiles and rotated back to the model frame (the stand fills the
box; nothing else is clipped), cross-tile pairs closer than 0.5 m are resolved by dropping the
synthetic member (none occurred), crown radius comes from the measured crown diameter for
anchors and from a species-specific linear crown-DBH fit otherwise, and root radius from the
power laws used in the Pambala reconstruction (Ceriops with Rhizophora coefficients). Output:
`demo_reconstructed.csv` (model schema; identical to `demo_reconstructed_model_species.csv`) and
`reconstructed_stand.csv` (with provenance: anchor / mtppr, transect, band, tile, seed); one
replicate stand per seed alongside.

## 8. Result

| | trees |
|---|---|
| total adults in 7.3 ha (clipped at the fringe; 2,276 before the clip) | 2,109 |
| measured (anchors) | 203 |
| MTPPR-generated | 1,906 |
| Rhizophora / Ceriops / Avicennia | 1,362 / 421 / 326 |
| density by model zone C / B / A (trees/ha) | 55 / 322 / 464 |

Zone composition follows the survey: zone A (landward) Rhizophora 50 %, Ceriops 34 %, Avicennia
16 %; zone B Rhizophora 75 %; zone C (the Sonneratia fringe, which holds few modelled trees)
Rhizophora 77 %, Avicennia 23 %, no Ceriops. Archived variants: the five-species recorded-tree
stand has 2,660 trees, the five-species distance-density stand 9,514.

## 9. Validation

- **Density and composition.** Tile densities per band (0.013, 0.033, 0.038, 0.048 trees/m2)
  match the band references (0.012, 0.035, 0.039, 0.056) within +9/-13 %; species shares per
  band are preserved from the anchors into the MTPPR-filled area (`fig_recon_profile.png`).
- **Traits.** Two-sample KS tests of r_stem and h_stem between MTPPR-generated and measured
  trees, per species: all p >= 0.14 (`species_marks_ks.csv`, `fig_recon_marks.png`).
- **Spatial structure.** Pair-correlation, mark-correlation, mingling and nearest-neighbour
  functions of the tiles against the reference belts (`fig_recon_summary_band.png`,
  `summary_function_fit.csv`): mean absolute deviations 0.09-0.17 for g(r) and <= 0.11 for
  k_mm(r); the fringe band reference holds only 17 trees, so its functions are noisy.
- **Seams.** Tree density within 2.5 m of tile boundaries is 0.92 (corridor seams) and 0.80
  (band seams) times the interior density (`seam_diagnostics.csv`).
- **Leave-one-transect-out cross-validation.** See section 10.

## 10. Leave-one-transect-out cross-validation

For each transect in turn, the belt references were rebuilt without it, every tile of its
corridor was reconstructed blind from the pooled references of the other transects (5 seeds),
and a virtual PCQM survey along its axis (same plots, quarters and 20 m search radius) was
compared with its real records. Real quarters whose nearest tree was a Sonneratia or Bruguiera
count as "no modelled tree within 20 m" (`tables/loto_results.csv`, `figures/fig_loto_*.png`).

| held out | plots | density: real (95 % CI) vs blind, trees/ha | seeds in CI | nearest modelled tree, real vs blind (m) | composition JSD | KS p (d / DBH / height) |
|---|---|---|---|---|---|---|
| T1 | 17 | 200 (90-710) vs 300 | 5/5 | 2.0 vs 5.8 | 0.02 | 0.00 / 0.47 / 0.11 |
| T2 | 16 | 190 (90-640) vs 270 | 5/5 | 2.2 vs 6.0 | 0.02 | 0.00 / 0.12 / 0.01 |
| T3 | 13 | 240 (100-420) vs 220 | 5/5 | 3.5 vs 6.6 | 0.05 | 0.00 / 0.43 / 0.10 |
| T4 | 6 | 720 (300-2,740) vs 410 | 5/5 | 2.3 vs 5.0 | 0.17 | 0.00 / 0.06 / 0.27 |
| T6 | 16 | 180 (80-530) vs 280 | 5/5 | 2.4 vs 6.0 | 0.04 | 0.00 / 0.14 / 0.11 |

Blind reconstructions reproduce the zonation of every held-out transect
(`fig_loto_profile.png`), the diameter distributions, and the censored PCQM density within its
wide confidence interval. The distance test is the sharper one: where the field survey found a
modelled-species tree as the nearest tree, it stood 2.0-3.5 m from the sampling point, whereas
the blind stand puts the nearest tree at 5.0-6.6 m. This is the expected signature of the
recorded-tree convention (section 5): the stand is 2-4 times sparser than the recorded
nearest-tree distances imply. The archived distance-density variant reproduces those distances
(density within the sampling uncertainty for four of five transects) at 9,514 trees.

## 11. Assumptions and limitations

1. Within a band, the stand is treated as spatially homogeneous alongshore (pooling of the five
   transects). Alongshore differences at the same fringe distance (e.g. the Avicennia patch
   crossed by T2 at 45-115 m versus Rhizophora at T1) are retained only inside the belts; the
   cross-validation quantifies the cost of this assumption.
2. The stand density follows the recorded trees per sampled strip (Pambala convention), which
   is 2-4 times below the density implied by the recorded nearest-tree distances; the model's
   fecundity constants were calibrated on `demo.csv` (556 trees) and need re-checking against
   this stand (2,276 trees).
3. Second-order structure beyond the nearest-neighbour scale cannot be estimated from PCQM
   data; the stand carries the short-range inhibition of a pattern that contains only recorded
   nearest trees.
4. Since 2026-09-20 Sonneratia is reconstructed (387 trees in the clipped stand, 71 % of zone C);
   only Bruguiera (4 recorded adults) is represented through its exclusion sectors alone. The
   S. alba root allometry is the Avicennia proxy (both pneumatophore species). The stem-mark
   scale was raised from 0.4 to 0.45 m so that the largest S. alba stem stays inside [0, 1].
   Adding the Sonneratia records widens the recorded-mode belts of the fringe bands, so the
   three original species hold ~300 fewer trees than in the 3-species stand (2,074 -> 1,771).
5. The stand fills the whole box by design. Seaward of each transect's fringe (up to 82 m in T3's
   corridor) and landward of its end the trees are extrapolated from the fringe and landward
   bands; T4's corridor seaward of its surveyed 60 m is extrapolated from the pooled bands.
6. The transect positions are digitised from a hand-traced overlay on a 1992 aerial photograph;
   the map spacing (70-78 m) is about 45 % wider than the tape record of the orthogonal transect
   T5 suggests, and the GPS fixes suggest a 15-20 % tighter spacing. The stand may therefore be
   stretched alongshore by up to that factor; along the gradient the tape lengths are exact.
7. Young trees (1997/1999 "Y" records) are not part of the stand: they are the model's validation
   target. The orthogonal transect T5 is not used for the stand; its juvenile quadrats remain
   usable for validation and it stays defined in `transect_geometry.py` for that purpose.
8. The model domain, elevation and transect geometry changed with this work, and the juvenile
   quadrats are now assigned to zones by position with the model's own rule
   (`zone_quadrats_by_position.py`: 1999 quadrats A 72, B 228, C 35; the vegetation-class
   assignment is kept for reference). With 100 m zone bands the middle zone holds most of the
   surveyed quadrats; zone C is the Sonneratia fringe (34 of its 35 quadrats are S. alba class)
   and zone A the Avicennia/mixed belt (71 of 72). All scenario, calibration and validation
   runs must be re-done before comparison with earlier results.

## References

- Cintron, G. & Schaeffer Novelli, Y. (1984) Methods for studying mangrove structure. In: Snedaker & Snedaker (eds) The mangrove ecosystem: research methods. UNESCO, Paris, 91-113.
- Cottam, G. & Curtis, J.T. (1956) The use of distance measures in phytosociological sampling. Ecology 37, 451-460.
- Dahdouh-Guebas, F. et al. (2004) Human-impacted mangroves in Gazi (Kenya): predicting future vegetation based on retrospective remote sensing, social surveys, and tree distribution. Marine Ecology Progress Series 272, 77-92.
- Mitchell, K. (2007) Quantitative analysis by the point-centered quarter method. arXiv:1010.3303.
- Ohser, J. & Stoyan, D. (1981) On the second-order and orientation analysis of planar stationary point processes. Biometrical Journal 23, 523-533.
- Pollard, J.H. (1971) On distance estimators of density in randomly distributed forests. Biometrics 27, 991-1002.
- Wudel, C., Schlicht, R. & Berger, U. (2023) Multi-trait point pattern reconstruction of plant ecosystems. Methods in Ecology and Evolution 14, 2668-2679.

### Forest edge and substrate (2026-09-06)

The stand used by the model is the whole-domain reconstruction clipped seaward of a forest edge. East of T1 the forest edge is the line through the seaward ends of the transects; west of T1 it is a line drawn by hand on the stand map (`forest_edge_user.csv`) that leaves the start of T4's surveyed part (tape 120 m, the only evidence for T4) hugging the left edge of the box and bends seaward to run about 10 m above the seaward boundary, so that the treeless band narrows towards T4 as the sand bank does on the printed map (band 1.09 ha, 82 % of it in zone C; 202 synthetic trees removed; 2,074 trees and 996 reproductive parents remain). The same band is the model's substrate map (`make_substrate_map.py` -> `substrate_map.npz`, two classes: forest floor where propagules root and parents stand, sand band where stranded propagules cannot root and no parent stands; Eq. 10 of the model), so parents and recruits obey one schematic physiography. The paper's vegetation classes (Dahdouh-Guebas et al. 2004, Fig. 6a, classified pixel by pixel by the nearest legend colour, transformed with the transect digitisation, binned into 1 m cells by majority vote, patches under 25 m2 merged into their neighbours; 14.7 % of the box lies outside the printed frame) are digitised only to show, in `fig_layout_substrate.py`, where the band sits relative to the mapped sand bank; they are not read by the model. All model-frame panels, the digitisation figure included, are drawn with the sea at the bottom (the model's convention); the printed map, sea at the top, therefore appears as the model panels turned upside down (the model frame is the printed map rotated by 17 degrees and seen from the other side), a display convention that does not affect the model.


## 12. Multi-transect hold-out and the adult-zonation evaluation (2026-09-20)

`run_holdout.py --holdout 3 1` rebuilds the stand from T4, T2 and T6 only: the prep step excludes
T3 and T1 entirely (no anchors, sectors or band-density records; 142 anchors remain), the corridor
and tile geometry is unchanged, and the seven tiles of the T3 and T1 corridors are reconstructed
blind (5 seeds). The blind stand is those tiles, post-processed like the published stand, merged
with the published stand's trees in the other corridors (`results/reconstruction/holdout_T1T3/`).
`adult_holdout_evaluation.py` surveys the blind stands with the field protocol along T3, T1 and
the orthogonal T5 (never used by any reconstruction; 26 of its 32 plots lie inside the box, 21 in
zone B and 5 in zone A, its Sonneratia belt at plots 29-32 lies east of the box), zones the
quarters by position, and compares the observed share of each species' adults per zone with the
simulated share (bootstrap over plots within zones and over the replicate stands, TVD with a 95 %
upper bound against the a priori 15-point margin, per-zone JSD, plot profiles). The distance test
fails by construction in recorded mode (section 10); DBH and height distributions are compared
with KS tests. Outputs: `results/evaluation/adult_holdout_*.csv`, `results/figures/fig_adult_zonation_holdout.png`,
`fig_holdout_map.png`.
