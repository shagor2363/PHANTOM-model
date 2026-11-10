# Methods text for the simulation experiments (2026-11-06)

Figure roles (user, 2026-11-06):
- MAIN: `results/figures/fig_endpoints_main.png` (transport-trait contrasts and stranding elevation) and
  `results/figures/fig_establishment_filter.png` (from stranded propagules to recruits).
- APPENDIX: `results/figures_final/fig_drag_zonation.png` (drag profile, recruit shares, four species) and
  `results/figures/fig_equal_cv.png` (one-at-a-time equal-contrast sweep; identical to `figures_final/fig_mechanisms.png`).

All numbers below were checked against the scripts and caches on 2026-11-06 (`figure_endpoints_main.py`,
`figure_establishment_filter.py`, `figure_drag_zonation.py`, `equal_cv_experiment.py`; caches `results/knockout_final`,
`results/knockout_cascade`, `results/knockout_son_shape`, `results/drag_sensitivity`, `results/equal_cv`).

---

## MAIN TEXT

### 2.X Simulation experiments

All experiments used the calibrated model (Table 2) on the reconstructed stand (2,119 trees). Each run covered 365 days
with both propagule sources active and recorded, for every propagule, its species, its origin (seaward boundary or
parent tree), its terminal position and its fate. Replicate runs differ only in the random seed.

### 2.X.1 Transport-trait contrasts and stranding elevation (Fig. A = fig_endpoints_main)

To test whether species differences in transport traits generate the separation of stranding elevations, we compared
two versions of the model. In the reference version each species kept its own draft, mobility and lateral spread
(Table 2). In the equalised version these three traits were set, simultaneously and for all four species, to the mean of
the three calibrated species (draft 0.108 m, mobility 46.7 m per tidal cycle, lateral spread 8.3 m per tidal cycle).
This removes the between-species contrast while leaving the mean transport capacity unchanged. Nothing else differed
between the versions: the parent stand, propagule supply and phenology, zone-specific vegetation drag, the transport
window, the settlement rule and the establishment submodel were identical. Each version was run with the same three
random seeds. Because every propagule carries its origin, the source configurations compared in Fig. A (boundary
origin, parent-tree origin, both pooled) are subsets of the same runs and not separate source-removal runs. The
response variable was the stranding elevation, i.e. the bed elevation (Eq. 1) at the terminal position of every
propagule that grounded inside the domain, irrespective of whether it later established or died. Propagules washed out
across the seaward boundary (< 0.1 %) were excluded. Elevations were pooled over the three seeds and summarised per
species, origin and model version as a Gaussian kernel density estimate (Scott's rule bandwidth) over 1.0-2.2 m, each
curve integrating to one. The pooled samples of the reference version held 39,600-89,901 boundary-origin and
145,771-286,505 parent-tree-origin propagules per calibrated species. Because the supply of Sonneratia alba was fitted
to its very low recorded juvenile density, the standard runs hold too few of its propagules to estimate a distribution.
Its curves were therefore taken from otherwise identical runs in which its supply from both sources was multiplied by
1000 (31,087 boundary-origin and 225,035 parent-tree-origin propagules). Propagules do not interact in the model, so
this changes the sample size of the species and not its distribution, and leaves the other species unaffected.

### 2.X.2 From stranded propagules to recruits (Fig. B = fig_establishment_filter)

To show how establishment reshapes the stranded distribution into the recruit distribution, we compared, within the
reference version, the elevations of all stranded propagules with those of the subset that established. A stranded
propagule is one that grounded inside the domain. It was counted as established if it survived the 90-day establishment
window (Section X.X): rooting fails on the treeless band seaward of the forest edge, after which the seedling is
exposed every day to burial by net sediment deposition, lethal once the cumulative burial reaches the tolerance of the
species, and to predation. The window is counted from the day of release, since transport lasts at most four days. All
stranded propagules that did not establish had died within the window. For each species, the stranding elevations of
the whole domain were counted in 4 cm classes between 1.0 and 2.2 m, and the class counts were averaged over 20
replicate runs. Both histograms were divided by the total number of stranded propagules of the species. The histogram
of stranded propagules therefore sums to 100 %, and the histogram of recruits sums to the overall establishment
percentage of the species. Zone-specific establishment was calculated as the established percentage of the propagules
stranded within that zone (zone C 1.0-1.4 m, zone B 1.4-1.8 m, zone A 1.8-2.2 m). For S. alba the histograms were
taken from the three supply-amplified runs described above. For each species we also report the share of its stranded
propagules that were of boundary origin and the share of its boundary-origin propagules that stranded on the treeless
band; these come from the standard runs.

---

## APPENDIX

### Appendix A. Sensitivity of recruit zonation to the vegetation drag profile (Fig. S-A = fig_drag_zonation)

Vegetation drag enters transport as a zone-specific retention factor F (1 = no reduction of movement; smaller values =
stronger drag) that multiplies the landward displacement, its random component and the lateral spread of every tidal
cycle. We ran the model with five profiles of F for zones C / B / A, ordered from strong drag to no drag:
0.2 / 0.6 / 1.0, 0.4 / 0.7 / 1.0 (the calibrated profile), 0.6 / 0.8 / 1.0, 0.8 / 0.9 / 1.0 and 1.0 / 1.0 / 1.0. All
other parameters kept their calibrated values, and each profile was run with three seeds. For each species and
propagule origin we calculated the percentage of the established recruits that were located in each zone, over the
whole domain, and averaged it over the seeds. The factors of zones C and B change together across the profiles, so the
experiment tests the strength of the seaward-to-landward drag gradient and not the factor of a single zone. S. alba was
run at its fitted supply here, which yields 0-4 boundary-origin and 24-38 parent-tree-origin recruits per run and none
of boundary origin under the strongest drag; its panels are indicative only.

| Profile (left to right in Fig. S-A) | F zone C | F zone B | F zone A |
|---|---|---|---|
| 1 (strong drag) | 0.2 | 0.6 | 1.0 |
| 2 (calibrated) | 0.4 | 0.7 | 1.0 |
| 3 | 0.6 | 0.8 | 1.0 |
| 4 | 0.8 | 0.9 | 1.0 |
| 5 (no drag) | 1.0 | 1.0 | 1.0 |

### Appendix B. Transport traits compared at equal relative contrast (Fig. S-B = fig_equal_cv, Table S-B)

The calibrated species differ far more in draft than in mobility (coefficient of variation across the three
calibrated species, CV: 0.55 for draft, 0.13 for mobility, 0.34 for lateral spread). Removing the calibrated contrasts
therefore compares contrasts of different size. To compare the traits at equal relative contrast we varied the
between-species contrast of one trait at a time, keeping the other two traits at their calibrated values. For a target
contrast c, the trait value of species i was set to

    T_i(c) = mean + (c / CV_cal) (T_i,cal - mean)

where the mean and CV_cal are taken over the three calibrated species. The cross-species mean is the same at every
level; at c = 0 all species share the mean, at c = CV_cal the calibrated values are recovered, and at c > CV_cal the
calibrated differences are stretched. Ten levels were run for each trait (c = 0, 0.045, 0.089, 0.134, 0.25, 0.344,
0.45, 0.552, 0.65 and 0.75), each with three seeds. The upper limit is the last level below c = 0.76, at which the
draft of Avicennia marina reaches zero. S. alba followed the same rule, with its draft set to zero above c = 0.68. The
response was the separation of stranding elevations, calculated per run as the mean of the pairwise Wasserstein-1
distances between the stranding-elevation distributions of the three calibrated species, separately for
boundary-origin propagules, parent-tree-origin propagules and both pooled. We report the change in separation relative
to the zero-contrast run of the same trait and seed, averaged over the three seeds. Levels above the calibrated
contrast of a trait are extrapolations made for the comparison and are not estimates for the site.

Table S-B. Transport-trait values in the equal-contrast experiment (Avicennia / Ceriops / Rhizophora / Sonneratia).

| Trait (unit) | All species equal (c = 0) | Calibrated values | Calibrated contrast CV_cal | Values at the upper limit (c = 0.75) |
|---|---|---|---|---|
| Draft (m) | 0.108 | 0.030 / 0.120 / 0.175 / 0.020 | 0.55 | 0.002 / 0.124 / 0.199 / 0.000 |
| Mobility (m per tidal cycle) | 46.7 | 55 / 45 / 40 / 55 | 0.13 | 93.4 / 37.3 / 9.3 / 93.4 |
| Lateral spread (m per tidal cycle) | 8.3 | 12 / 8 / 5 / 12 | 0.34 | 16.3 / 7.6 / 1.1 / 16.3 |

Each species moves monotonically from the common mean (c = 0) through its calibrated value (c = CV_cal) to the value at
the upper limit. The column 'All species equal' also gives the values of the equalised model version of Section 2.X.1,
where the three traits are equalised together.

---

## Table 2. Species-specific parameters (current calibrated model, checked against `tidal_dispersal_model.py` on 2026-11-09)

Status: lit = literature value; cal = calibrated on the 1999 juvenile zonation of the three calibrated species;
asm = assumed (Sonneratia alba, never calibrated); fit = Sonneratia supply multiplier s = 3.0e-4 fitted to its
mean recorded 1999 juvenile density. Sources for the literature values: `docs/parameter_provenance.md`,
`docs/species_parameters_sonneratia.md`.

| Parameter | Avicennia | Rhizophora | Ceriops | Sonneratia | Status |
|---|---|---|---|---|---|
| Mobility M (m per tidal cycle) | 55 | 40 | 45 | 55 | cal; Son asm (= Avicennia) |
| Lateral spread sigma (m per tidal cycle) | 12 | 5 | 8 | 12 | lit; Son asm (= Avicennia) |
| Draft (m) | 0.03 | 0.175 | 0.12 | 0.02 | lit mid-range; Son asm (fruit size) |
| Initial height (cm) | 1.4 | 18.0 | 8.25 | 2.0 | lit |
| Daily growth (cm per day) | 0.30 | 0.25 | 0.22 | 0.05 | lit |
| Burial tolerance (mm) | 30 | 80 | 50 | 20 | lit-inferred; Son asm |
| Vulnerability V | 0.85 | 0.25 | 0.55 | 0.85 | lit ordering, values set; Son asm (= Avicennia) |
| H_ref (cm) | 10 | 30 | 28 | 8 | asm |
| F_base (propagules per year, tree at r_min) | 12.8 | 8 | 48 | 12.8 | cal; Son set equal to Avicennia (asm) |
| gamma (allometric exponent) | 2.0 | 2.5 | 2.5 | 2.0 | asm; Son = Avicennia |
| Maturity r_stem (m) | 0.027 | 0.024 | 0.027 | 0.025 | asm; Son lit (Gazi plantation DBH, upper bound) |
| Boundary influx, spring / neap (propagules per day) | 448 / 224 | 160 / 64 | 384 / 192 | 448 / 224 | cal; Son set equal to Avicennia (asm) |
| Survival from dispersal unit to a seedling entering the establishment window | 1 | 1 | 1 | 3.0 x 10^-4 | viviparous species 1 by definition; Son fitted to its mean 1999 juvenile density (3 occupied quadrats; resolved to about a factor of 2) |
| Fruiting window (day of year: shoulder start, peak start, peak end, shoulder end) | 60, 91, 152, 182 | 152, 244, 304, 335 | 335, 1, 59, 90 | 213, 244, 319, 349 | lit (Gazi phenology) |

Common to all species: zone retention factor F = 0.4 / 0.7 / 1.0 (zones C / B / A, cal); transport window 8 spring /
4 neap tidal cycles (cal); settlement draft coefficient 6.0 (cal); siltation R_base 0.65 / 0.3 / 0.1 mm per day
(zones C / B / A); predator density 24.0 / 7.97 / 4.24 (zones C / B / A); base encounter rate 0.003; daily predation
capped at 0.15; 90-day establishment window.

Corrections to the draft of Table 2 supplied on 2026-11-09: mobility Avicennia 70 -> 55 and Ceriops 40 -> 45
(mobility experiment 2026-09-06); boundary influx 112 / 40 / 96 and 56 / 16 / 48 -> 448 / 160 / 384 and
224 / 64 / 192 (x4 recalibration 2026-09-05; the x2 variant of 2026-09-22 was NOT adopted); the title
'All values calibrated for Gazi Bay' is not defensible, use the status column; Sonneratia column added.

Footnote to Table 2: for S. alba the model applies the survival factor to the supply constants at release, so the
effective values held in the code are F_base 0.00384 per year and a boundary influx of 0.1344 / 0.0672 propagules per
day (about 10 per year, drawn by stochastic rounding). Propagules do not interact, so a factor applied at release is
equivalent to a day-0 survival draw at stranding. Calibration procedure: docs/parameter_provenance.md (F_BASE row)
and the account of 2026-11-09 (son_bracket.py, son_refit.py; target 6.15e-4 per m2 over 309 in-domain quadrats;
realised 6.2e-4 at s = 3.0e-4, seed range 3.9-9.1e-4; lowest level tested; 1997 would give s near 1.2e-3).

Sonneratia boundary influx: the daily rates 0.1344 / 0.0672 are the calibrated Avicennia rates x s. Releases are
drawn by stochastic rounding, so a day yields 0 or 1 propagule; the expected annual total is 10.4 and the
three baseline runs realised 11 / 4 / 14 boundary-origin and 74 / 74 / 71 parent-tree Sonneratia propagules
(seeds 42 / 49 / 56, 365 d). Give the annual row in the table and the per-day values in the footnote.
This is the 'Table X' referred to in Sections 2.X and 2.X.1.

---

## 2.X Environmental parameters (text and Table 3, 2026-11-09)

The spatial frame, the elevation field and the tidal forcing were fixed before calibration and were not adjusted to
the recruit data. The 325 x 225 m domain is the envelope of the five parallel survey transects on the published site
map (Dahdouh-Guebas et al., 2004), rotated so that the transects run along the land-sea axis. No elevation model
exists for the site, so bed elevation was represented as a planar gradient from 1.0 m at the seaward corner to 2.2 m
at the landward corner, oriented along the surveyed mangrove belt, and the three zones were defined as equal thirds of
this gradient (boundaries at 1.4 and 1.8 m). The tide was represented by a semi-diurnal harmonic (M2 period 12.42 h)
about a mean level of 1.6 m with amplitudes of 1.25 m at spring and 0.5 m at neap tide over a 15-day spring-neap
cycle. Under this forcing the seaward edge is inundated for 16 h per spring day and continuously at neap, the landward
edge for 9 h per spring day and not at all at neap, which gives the seaward-to-landward hydroperiod gradient that
drives stranding, siltation and predation.

The zone-specific constants of the establishment submodel, the siltation rate at full inundation (0.65 / 0.3 / 0.1
mm per day in zones C / B / A) and the effective predator density (24.0 / 7.97 / 4.24), were fixed during model
development so that burial and predation pressure decrease landward and seedling mortality in the seaward zone lies in
the range reported for exposed mangrove fringes (40-70 %). They were then held constant through all subsequent
calibration steps. Because the siltation rate acts on cumulative burial over the 90-day window, the zone C value is
the constant that decides which species the seaward zone excludes: at this rate the window accumulates about 37 mm,
above the burial tolerance of A. marina and S. alba and below that of C. tagal and R. mucronata. The rate should be
read as the short-term bed-level change experienced by a seedling, not as the net surface accretion of the forest.

The only environmental constant calibrated to the recruit data was the vegetation retention factor F. The seaward and
middle zones were calibrated jointly, with the landward zone fixed at 1 (open flow), by minimising the Jensen-Shannon
divergence between observed and simulated 1999 zone shares of the three calibrated species over four of the five
parallel transects, with the fifth (T2) withheld. The retained profile (0.4 / 0.7 / 1.0) could not be distinguished
from its neighbours within the replicate noise of the objective; the sensitivity of recruit zonation to the full range
of profiles is shown in Appendix A.

Sensitivity sentence, version without new runs: A Morris elementary-effects screening of 36 parameters (10
trajectories, 1,110 runs of 365 days; `results/sensitivity`, March 2026, development version: schematic domain, three
species, old drag orientation) ranked the tidal amplitudes and the zone boundaries among the least influential
parameters for the zone composition error and the zone C siltation rate among the most influential; because the
screening predates the final geometry and constants, we report it as a structural result only.
Version after an OAT on the frozen configuration (proposed, 72 runs, not made): no perturbation within the stated
ranges changed the zone ordering of any species or the exclusion of A. marina and S. alba from zone C.

Table 3. Environmental parameters.

| Parameter | Value | Basis |
|---|---|---|
| Domain extent | 325 x 225 m (0-325 alongshore, 0-225 seaward to landward) | published site map |
| Bed elevation, seaward to landward corner | 1.0 to 2.2 m, planar | assumed, no DEM |
| Zone thresholds C/B and B/A | 1.4 / 1.8 m | thirds of the gradient |
| Mean water level H0 | 1.6 m | fixed a priori |
| Tidal amplitude, spring / neap | 1.25 / 0.5 m | fixed a priori |
| M2 period | 12.42 h | tidal constant |
| Spring-neap cycle | 15 days, days 1-7 spring | fixed a priori |
| Retention factor F (C / B / A) | 0.4 / 0.7 / 1.0 | C and B calibrated, A fixed at 1 |
| Siltation R_base (C / B / A) | 0.65 / 0.3 / 0.1 mm per day | fixed in development to the zone C mortality range |
| Predator density D (C / B / A) | 24.0 / 7.97 / 4.24 | fixed in development to the zone C mortality range |

Corrections to the user's draft of 2026-11-09: domain extent was (-30, -50) to (220, 260) (the pre-2026-09-04
schematic domain); drag was 1.0 / 0.7 / 0.3 (reversed orientation and stale values).
Open before submission: (1) the 40-70 % zone C mortality range has no verified source (not in Balke et al. 2011,
citation audit 2026-09-22); (2) the tide amplitudes have not been compared with Gazi tide data (Kitheka); the model's
spring range is 2.5 m, neap 1.0 m; (3) the implied zone C burial (~150 mm per year) must be presented as seedling-scale
bed-level change, with Balke et al. (2013, Biogeosciences, Gazi cross-shore disturbance gradients) to check as support;
(4) the seaward-high predator density contradicts the Gazi predation studies (mid-intertidal maximum or no pattern) and
must be called a modelling choice.

---

## 2.X Species parameters: how they were set (text, 2026-11-09)

The transport traits of the three calibrated species were taken from the literature where measurements exist. Draft
was set to the mid-point of the published range for each species' propagule, and the lateral spread per tidal cycle to
published cross-shore dispersal estimates; neither was adjusted to the recruit data. Net landward mobility per tidal
cycle is not measured for the species, so it was calibrated: starting from literature-based values, we ran a grid of
mobilities (A. marina 30-70, C. tagal 30-40, R. mucronata 20-40 m per cycle, two seeds each, 365 days;
`mobility_experiment.py`, 2026-09-06) and scored each set on the 1999 zone shares of the three species together with
the three field signatures used as guardrails. The data constrained R. mucronata to 40 m per cycle (lower values
overfilled the seaward zone) and left A. marina between 30 and 50 and C. tagal between 30 and 40 indistinguishable;
among the equivalent sets (confirmation runs 50/45/40, 55/45/40, 55/50/40) we retained 55 / 45 / 40, which preserves the
size-based ordering A. marina > C. tagal > R. mucronata at no cost to the fit; the zonation data cannot detect mobility
differences of this size. The common transport constants (transport window 8 / 4 cycles, cycles afloat before
settlement 1, settlement draft coefficient 6) were calibrated in the same way before the mobility grid
(`calibrate_geom2.py` pass 1, 2026-09-05).

None of the species-specific establishment constants was calibrated. Initial seedling height and daily growth are
literature values; burial tolerances follow the published ordering of sediment tolerance with values inferred from
burial experiments (attribution weak, citation audit 2026-09-22); predation vulnerabilities follow the ordering reported
for the species, with the values chosen to span it; the reference heights at which predation protection saturates are
assumed. The supply constants define the propagule budget of the stand and the seaward boundary. The per-tree output at
the minimum reproductive size (F_base) and the allometric exponent were fixed during development as effective values
(no derivation; ordering Cer > Avi > Rhi does not follow the ODD literature fecundities); the species proportions of the
boundary influx (448 : 384 : 160) were fixed during development as an assumed regional source composition, and the
overall scale of the influx was calibrated jointly with the retention factor (x4, calibrate_geom2 pass 2, 2026-09-05;
x2 tested 2026-09-22 and not adopted), because it controls the share of boundary-origin recruits that separates the
combined-source model from the parent-tree-only model. For S. alba every constant was set by rule (transport,
vulnerability and per-tree output equal to A. marina; draft from fruit size; growth, initial height, burial tolerance
and reproductive size from the literature) and one factor s was fitted to its mean recorded 1999 juvenile density.

Status column for Table 2: mobility calibrated (3 species; Son = Avi asm); lateral spread, draft lit (Son draft asm);
initial height, daily growth lit; burial tolerance lit-inferred; vulnerability lit ordering, values set; H_ref asm;
F_base, gamma fixed in development (effective; Son = Avi); maturity fixed in development (Son lit); boundary influx
proportions fixed in development, overall scale calibrated; s fitted (Son only).
The Morris screen of March 2026 varied mobility, lateral spread, F_base and the boundary rates but NOT the burial
tolerances, vulnerabilities, growth or heights; no sensitivity record exists for those.
User's draft of 2026-11-09 still has Ceriops mobility 40; the model has 45.

---

## Results / appendix text for fig_drag_zonation_all_SUPPL.png and fig_equal_cv.png (2026-11-09)

### Fig. S-A (fig_drag_zonation_all_SUPPL.png): recruit zonation of all recruits under five drag profiles

Text. Fig. S-A shows how the zone composition of all established recruits of each species responds to the vegetation
drag profile, in the whole domain (top row) and inside the survey belts seen by the 1999 survey (bottom row). The five
profiles (Table S-A1) weaken the seaward-to-landward drag gradient from left to right; the calibrated profile is the
second point. Zone C recruits of A. marina and S. alba are absent under every profile because burial in the seaward
zone exceeds their tolerance within the 90-day window; their lines lie on zero. In the whole domain the composition of
every species changes by at most 10 share points across the full range of profiles (Table S-A2): weakening drag moves
a few per cent of C. tagal and R. mucronata recruits from zone C to zone A, and A. marina changes only between the
strongest and the calibrated profile. The reason is that 62-88 % of each species' recruits originate from parent trees,
whose propagules strand within a few metres of the parent and are therefore insensitive to drag; the drag response of
the boundary-origin recruits, which is large (main-text figure), is diluted when the two origins are pooled. Inside
the survey belts the response is larger, up to 16 points for C. tagal (zone A share 84 -> 68 %) and 10 points for
R. mucronata (zone C 2 -> 12 %), because the belts contain a higher share of boundary-origin recruits near the transect
fringes; the belt composition is the quantity the calibration targets, which is why the drag profile could be
calibrated at all, and why neighbouring profiles were hard to separate. S. alba is drawn from 30-37 recruits per run
at its fitted supply and its lines are indicative only. Seed-to-seed ranges are below 1 point in the whole domain and
below 3 points in the belts. The figure therefore shows that the drag profile is a weak lever on the zonation of all
recruits and a strong lever on boundary-origin recruits only; the zonation conclusions of the paper do not depend on
the profile chosen.

Table S-A1. Drag profiles (movement retention factor F; 1 = no reduction of movement, smaller = stronger drag).

| Profile (left to right) | F zone C | F zone B | F zone A |
|---|---|---|---|
| 1, strong drag | 0.2 | 0.6 | 1.0 |
| 2, calibrated | 0.4 | 0.7 | 1.0 |
| 3 | 0.6 | 0.8 | 1.0 |
| 4 | 0.8 | 0.9 | 1.0 |
| 5, no drag | 1.0 | 1.0 | 1.0 |

Table S-A2. Zone shares (%) of all established recruits, C / B / A, mean of three seeds; whole domain | survey belts.
Mean recruits per run in the whole domain in the last column (profile 1 -> 5).

| Species | Profile 1 | 2 (calibrated) | 3 | 4 | 5 | n per run (domain) |
|---|---|---|---|---|---|---|
| A. marina | 0/28/72 \| 0/24/76 | 0/38/62 \| 0/33/67 | 0/37/63 \| 0/30/70 | 0/35/65 \| 0/27/73 | 0/35/65 \| 0/28/72 | 33,500 -> 48,500 |
| C. tagal | 6/43/51 \| 3/13/84 | 6/40/53 \| 7/15/78 | 6/39/55 \| 9/18/73 | 5/38/57 \| 9/22/69 | 4/38/58 \| 7/25/68 | 47,600 -> 55,200 |
| R. mucronata | 12/57/31 \| 2/56/42 | 11/57/32 \| 8/51/41 | 10/58/32 \| 12/48/39 | 9/58/33 \| 12/48/40 | 8/58/34 \| 11/49/40 | 90,500 -> 95,300 |
| S. alba (indicative) | 0/67/33 \| 0/79/21 | 0/67/33 \| 0/79/21 | 0/66/34 \| 0/77/23 | 0/55/45 \| 0/84/16 | 0/60/40 \| 0/84/16 | 30 -> 37 |

Source: results/drag_sensitivity/drag_zonation_shares.csv (stage recruits, stratum all, operators domain / belt).
Note: the n column of the belt operator in the CSV is the total over species, not per species.

### Fig. S-B (fig_equal_cv.png): transport traits compared at equal relative contrast

Text. Fig. S-B shows the separation of stranding elevations between the three calibrated species as a function of the
between-species contrast of one transport trait, the other two traits held at their calibrated values. The x axis is
the target coefficient of variation c of the trait across the three species (Table S-B1 gives the species values at
each level); the y axis is the change in the mean pairwise Wasserstein-1 distance between the species' stranding
elevation distributions relative to the zero-contrast run of the same trait and seed, so every curve starts at zero.
Three results. First, lateral spread has no effect in any source configuration: the curve is zero to the millimetre at
every level, because lateral displacement is alongshore and does not change elevation. Second, in the boundary-origin
stratum (A) both draft and mobility separate the species linearly with contrast, mobility about twice as steeply per
unit contrast (0.21 against 0.12 m per unit CV; Table S-B2). The species, however, differ four times more in draft
(CV 0.55) than in mobility (CV 0.13), so at the contrasts they actually have, draft contributes about twice the
separation of mobility (67 against 31 mm). The two readings are consistent: mobility is the more effective trait per
unit contrast, draft the more important one at the site. Third, in the parent-tree stratum (B) the same ordering holds
at about 40 % of the boundary slopes (26 and 10 mm at the calibrated contrasts), because parent-tree propagules
strand close to where they are released and transport traits can only shift them slightly. In the full configuration
(C) the pooled separation is not the sum of A and B: increasing the draft contrast lowers it by up to 10 mm, and the
mobility contrast lowers it up to c = 0.35 and raises it beyond. Each species' pooled distribution is a mixture of a
seaward boundary-origin mode and a landward parent-tree mode in species-specific proportions (boundary share 0.38 /
0.29 / 0.12 for A. marina / C. tagal / R. mucronata); the trait contrasts move the boundary mode of A. marina landward
into the parent-tree modes of the other species, so the within-stratum separations rise while the pooled separation
falls. Panel C therefore does not show that transport matters less in the full model; it shows that in the full model
the between-species separation is set by the sources, and that transport contrasts redistribute propagules between
each species' seaward and landward modes. Mechanism statements use panels A and B. Levels above a trait's calibrated
contrast are extrapolations for the comparison, not estimates for the site. Seed ranges are below 2 mm everywhere.

Table S-B1 = Table S-B of Appendix B above (trait values at c = 0, calibrated, and c = 0.75).

Table S-B2. Change in between-species separation (mean pairwise W1 of stranding elevation, calibrated species, mm)
relative to zero contrast, mean of three seeds. Columns: equal contrast c = 0.134 (the calibrated contrast of
mobility), the trait's own calibrated contrast (draft 0.552, mobility 0.134, spread 0.344), and the upper limit
c = 0.75. Slope = mm per unit CV in the boundary stratum, fitted through the origin.

| Stratum | Trait | c = 0.134 | own calibrated contrast | c = 0.75 | slope (mm per unit CV) |
|---|---|---|---|---|---|
| A. Boundary origin | Draft | +16 | +67 | +90 | 120 |
| | Mobility | +31 | +31 | +159 | 212 |
| | Lateral spread | 0 | 0 | 0 | 0 |
| B. Parent-tree origin | Draft | +6 | +26 | +36 | 48 |
| | Mobility | +10 | +10 | +60 | 80 |
| | Lateral spread | 0 | 0 | 0 | 0 |
| C. Both sources pooled | Draft | -3 | -10 | -10 | |
| | Mobility | -4 | -4 | +32 | |
| | Lateral spread | 0 | 0 | 0 | |

Zero-contrast baselines (absolute W1, m): boundary 0.023 (draft run) / 0.058 (mobility run) / 0.089 (spread run);
parent-tree 0.098 / 0.113 / 0.123; pooled 0.079 / 0.073 / 0.069. The baselines differ between traits because the other
two traits stay calibrated, which is why the figure plots increments. Source: results/equal_cv/equal_cv_w1.csv
(column dw1_calibrated).

---

## Simple methods text: drag sensitivity analysis and mobility calibration (2026-11-09)

### Vegetation drag: sensitivity analysis (`drag_sensitivity.py`, results/drag_sensitivity)

1. The calibrated model was taken as is (Tables 2 and 3). Only the retention factors F of the three zones were changed.
2. Five profiles were run, ordered from strong to no drag. In every profile zone A keeps F = 1 (open flow); zones C and
   B are weakened together so that the seaward-to-landward gradient becomes flatter step by step. A sixth, reversed
   profile (drag strongest landward, the orientation used before 2026-09-20) was run as a check and is not drawn.
3. Each profile was run for 365 days with three random seeds (42, 49, 56); the calibrated profile uses the three
   baseline runs of the final pipeline with the same seeds.
4. From every run we recorded, per species and origin (boundary / parent tree), the zone shares of stranded
   propagules and of established recruits over the whole domain and inside the survey belts, plus the 1999 zonation
   error of the calibrated species, the share of boundary propagules that cross zone C into zone B, and the zone C
   recruit count.
5. The figure plots the zone shares against the profile order; the table below gives the summary statistics.

Table S-A3. Drag profiles tried and their outcome (mean of three seeds).

| Profile | F zone C | F zone B | F zone A | 1999 zonation error | boundary propagules crossing C -> B | zone A share Avi / Cer / Rhi |
|---|---|---|---|---|---|---|
| 1 strong drag | 0.2 | 0.6 | 1.0 | 0.0143 | 2 % | 0.76 / 0.84 / 0.42 |
| 2 calibrated | 0.4 | 0.7 | 1.0 | 0.0086 | 19 % | 0.66 / 0.79 / 0.40 |
| 3 | 0.6 | 0.8 | 1.0 | 0.0176 | 26 % | 0.70 / 0.73 / 0.39 |
| 4 | 0.8 | 0.9 | 1.0 | 0.0226 | 33 % | 0.73 / 0.69 / 0.40 |
| 5 no drag | 1.0 | 1.0 | 1.0 | 0.0207 | 40 % | 0.72 / 0.68 / 0.40 |
| reversed (check only) | 1.0 | 0.6 | 0.2 | 0.0210 | 46 % | 0.77 / 0.75 / 0.33 |

Observed 1999 zone A shares: 0.54 / 0.80 / 0.37. Source: results/drag_sensitivity/drag_sensitivity.csv.
Reading: the profile acts through the boundary propagules' crossing from zone C into zone B (2 % to 46 %); it moves
the zone A shares of Avicennia and Ceriops by up to 10-16 points and leaves Rhizophora unchanged except under the
reversed profile. The calibrated profile has the lowest error, but profiles 1 and 3 are within two to three times the
seed noise of the objective, which is why the earlier calibration could not separate neighbouring profiles.

How the calibrated profile itself was chosen (`recalibrate_4T.py`, 2026-09-22; `recalibrate_T2blind.py`, 2026-09-24):
a grid of the same four gradient profiles (0.2/0.6, 0.4/0.7, 0.6/0.8, 0.8/0.9 for zones C/B, A fixed at 1.0) crossed
with the boundary-influx scale (x1, x2, x4), seed 42, 365 days, on the T2-withheld stand; objective = mean
Jensen-Shannon divergence between observed and simulated 1999 zone shares of the three calibrated species. The
0.4/0.7/1.0 profile at x1 was retained; the x2 variant and 0.2/0.6/1.0 at x4 were within seed noise of it (confirmed
with seeds 49, 56) and were not adopted.

### Mobility: calibration grid (`mobility_experiment.py`, results/calib_mobility, 2026-09-06)

1. Mobility (net landward displacement per tidal cycle) has no field measurement for the species, so it was
   calibrated. All other constants were held at their values of the time.
2. A full grid of species mobilities was run: Avicennia {30, 40, 50, 60, 70} x Ceriops {30, 40} x Rhizophora {20, 30, 40}
   m per cycle = 30 sets, plus three confirmation sets (50/45/40, 55/45/40, 55/50/40) = 33 sets, each with two seeds
   (42, 49), 365 days.
3. Every run was scored on (a) the 1999 zonation error (mean Jensen-Shannon divergence of the three species' zone
   shares; the objective), (b) the field-signature guardrails (Ceriops most landward, Rhizophora most seaward with its
   peak in zone B, Avicennia with the largest landward shift) and the zone mortality gradient, and (c) reported only:
   the 1997 error, the recruit density per zone relative to the survey, and the parent-to-recruit displacement.
4. Result: Rhizophora had to stay at 40 (30 and 20 overfilled zone C and failed the Rhizophora share); Avicennia 30, 40
   and 50 were equivalent on the 1999 error (0.0067-0.0080) and 60-70 were worse (0.009-0.012, zone A overfilled);
   Ceriops 30 and 40 were equivalent. Among the equivalent sets 55/45/40 was chosen because it keeps the ordering
   Avicennia > Ceriops > Rhizophora suggested by propagule size at no measurable cost (error 0.0085 against 0.0067 for
   the best set, a difference within seed noise).
5. Mobility was then held fixed. When the drag orientation was reversed and the stand changed (2026-09-20 to 09-22),
   five mobility sets {(55,45,40), (40,35,40), (35,30,40), (30,25,40), (30,25,30)} were re-tested inside the drag
   recalibration grid and 55/45/40 was kept.

Table S-M. Mobility sets tried (m per tidal cycle) and the 1999 zonation error (mean of two seeds, lower is better).
Rows = Avicennia mobility, columns = Ceriops / Rhizophora mobility. Blank = combination not run. All 33 sets passed
the signature and mortality guardrails; the zone C recruit density was 3x the survey at Rhi 40, 4x at Rhi 30 and 6x at
Rhi 20, in every row.

| Avicennia | Cer 30 / Rhi 20 | Cer 30 / Rhi 30 | Cer 30 / Rhi 40 | Cer 40 / Rhi 20 | Cer 40 / Rhi 30 | Cer 40 / Rhi 40 | Cer 45 / Rhi 40 | Cer 50 / Rhi 40 |
|---|---|---|---|---|---|---|---|---|
| 30 | 0.0191 | 0.0114 | 0.0077 | 0.0203 | 0.0111 | 0.0067 | | |
| 40 | 0.0203 | 0.0115 | 0.0077 | 0.0188 | 0.0111 | 0.0075 | | |
| 50 | 0.0209 | 0.0118 | 0.0079 | 0.0185 | 0.0118 | 0.0080 | 0.0074 | |
| 55 | | | | | | | **0.0085** (chosen) | 0.0094 |
| 60 | 0.0216 | 0.0136 | 0.0105 | 0.0221 | 0.0133 | 0.0092 | | |
| 70 | 0.0233 | 0.0151 | 0.0111 | 0.0240 | 0.0147 | 0.0116 | | |

Full table (33 sets, all scores): results/calib_mobility/mobility_ranking.csv.

Caveat for the text: the mobility grid was run on the model configuration of 2026-09-06 (five-transect stand, three
species, the pre-reversal drag orientation). The retained set was re-tested but not re-optimised after the drag
reversal and the change of stand; the 1999 data could not discriminate between the equivalent sets then and cannot
now.

---

## Table 1 (user's layout), corrected version for fig_equal_cv.png (2026-11-09)

Caption. Transport trait values in the one-at-a-time equal-contrast experiment. One trait at a time was varied while
the other two kept their calibrated values. At zero contrast all four species share the mean of the three calibrated
species. The between-species contrast was then raised in ten levels of the target coefficient of variation c (0, 0.045,
0.089, 0.134, 0.25, 0.344, 0.45, 0.552, 0.65, 0.75) by moving each species' value away from the common mean in
proportion to its calibrated deviation, T_i(c) = mean + (c / CV_cal)(T_i,cal - mean), so species above the mean
increase and species below it decrease, and the mean is unchanged at every level. The upper limit 0.75 is the last
level before the draft of Avicennia reaches zero (c = 0.76). The calibrated values are recovered at c = CV_cal (draft
0.55, mobility 0.13, lateral spread 0.34). Each level was run with three seeds.

| Trait varied (unit) | All species equal (c = 0) | Species | Range across the ten levels (c = 0 -> c = 0.75) | Calibrated value | Remarks |
|---|---|---|---|---|---|
| Draft (m) | 0.108 | Avicennia | 0.108 -> 0.002 | 0.030 | Upper limit of the experiment set where the smallest draft reaches zero |
| | | Ceriops | 0.108 -> 0.124 | 0.120 | |
| | | Rhizophora | 0.108 -> 0.199 | 0.175 | |
| | | Sonneratia | 0.108 -> 0.000 | 0.020 | Reaches zero at c = 0.68, held at zero above |
| Mobility (m per tidal cycle) | 46.7 | Avicennia | 46.7 -> 93.4 | 55 | Levels above c = 0.13 stretch the calibrated contrast up to 5.6x; model test, not site estimates |
| | | Ceriops | 46.7 -> 37.3 | 45 | |
| | | Rhizophora | 46.7 -> 9.3 | 40 | |
| | | Sonneratia | 46.7 -> 93.4 | 55 | = Avicennia by rule |
| Lateral spread (m per tidal cycle) | 8.3 | Avicennia | 8.3 -> 16.3 | 12 | Levels above c = 0.34 stretch the calibrated contrast up to 2.2x; same caveat |
| | | Ceriops | 8.3 -> 7.6 | 8 | |
| | | Rhizophora | 8.3 -> 1.1 | 5 | |
| | | Sonneratia | 8.3 -> 16.3 | 12 | = Avicennia by rule |

Footnote: common value and CV_cal over the three calibrated species; Sonneratia follows the same rule but its draft is
an assumption, so the separation metric of the figure uses the three calibrated species only.
Corrections to the user's draft: 17 levels -> 10; "0-75 %" -> target CV 0-0.75; draft common value 0.10 -> 0.108;
Rhizophora draft 0.17 / 0.19 -> 0.175 / 0.199; lateral-spread rows for Ceriops and Rhizophora added; Sonneratia added;
units added.

---

## Justification of the Sonneratia supply-amplified (x1000) runs (2026-11-09)

Argument: (1) the figures show normalised distributions, not amounts; (2) propagules do not interact anywhere in the model
(drag per zone, stranding by own draft and local depth, rooting per cell, burial and predation per individual with
zone-constant rates; no density dependence), so a propagule's stranding position is a draw from a distribution independent
of how many are released; both sources were multiplied by the same factor, so the boundary/parent mix is unchanged; (3)
empirical check against the 20 standard runs pooled (cascade_D1T1S1E1_seed*.csv) vs the 3 shape runs:
boundary origin n 205 vs 31,198, W1 1.0 cm, KS p 0.87, zone shares 65/34/1 vs 62/37/0 %, recruit share 0.18 vs 0.15;
parent-tree origin n 1,491 vs 225,035, W1 0.8 cm, KS p 0.18, zone shares 14/63/23 vs 15/63/22 %, recruit share 0.42 vs 0.42;
other species in the shape runs differ from the standard runs of the same seeds by < 0.1 cm W1. (4) Amplified runs are
used for distribution figures only; densities, recruit fractions and the observed-vs-simulated comparison use the standard
runs. Why not pool the 20 standard runs: 205 boundary propagules is too few for a KDE and the knockout treatments exist for
three seeds only.
Sentence for the paper: Because propagules do not interact in the model, the stranding distribution of a species is
independent of its supply; the amplified runs reproduce the Sonneratia distributions of the 20 standard runs within 1 cm
(Wasserstein-1; Kolmogorov-Smirnov p = 0.87 and 0.18 for boundary- and parent-tree-origin propagules) and leave the other
species unchanged.

---

## Methods paragraph "Species-specific parameters", revised with Sonneratia (2026-11-09)

Species-specific parameters represented propagule transport, establishment and reproduction (Table A1). For the
three focal species, propagule draft and lateral spread were assigned from published observations of propagule
morphology, buoyancy and dispersal (De Ryck et al., 2012; Van der Stocken et al., 2015, 2019) and were not adjusted
to the recruitment data. Effective landward mobility per tidal cycle was calibrated against the observed 1999
recruitment zonation of the three focal species, as comparable species-specific measurements were unavailable. A
grid of mobility combinations was evaluated on the recruitment zonation error, subject to additional constraints on
simulated transport and mortality patterns. Several combinations produced comparable agreement with observations,
indicating limited identifiability of species-specific mobility differences. Among these, we retained the combination
preserving the assumed size-based mobility ordering of A. marina > C. tagal > R. mucronata (55, 45 and 40 m per
cycle, respectively).

S. alba was added as a fourth species after this calibration and was not calibrated. Its dispersal unit is the small
buoyant seed released from the fruit, for which no draft or dispersal measurements comparable to those of the focal
species exist. Its draft (0.02 m) was therefore assumed from fruit size, and its mobility and lateral spread were set
equal to those of A. marina (55 m and 12 m per cycle), the species with the most similar propagule size. Its
establishment parameters were taken from the literature where available and otherwise assumed (Table A1). Only its
propagule supply was fitted: the per-tree fecundity and boundary influx were set by the rule used for A. marina and
multiplied by a single survival factor (3.0 x 10^-4) chosen so that the simulated S. alba juvenile density matched the
mean density recorded in the 1999 survey. The S. alba recruitment pattern is therefore reported as a model outcome
and was not used in calibration.

Notes: the mobility grid (2026-09-06) ran on the three-species model, before Sonneratia was added (2026-09-20).
The fitted density target is the mean over the 309 in-domain quadrats, 6.1e-4 per m2 (three occupied quadrats);
the fit resolves the factor to within about a factor of two and 3.0e-4 is the lowest level tested (a bound).
Citations De Ryck 2012 / Van der Stocken 2015, 2019 still to be verified by the authors.

---

## Short methods paragraph: stand reconstruction (2026-11-09)

The adult stand was reconstructed from the single point-centred quarter method (PCQM) census of the 1997 to 1999
campaigns, in which the nearest adult in each of four quarters was recorded every 10 m along each transect with its
distance, species, DBH and height. The 195 recorded trees of T1, T3, T4 and T6 were placed at their measured
distance from the tape at a random angle within their quarter; T2 was withheld as a blind test and T5 was never
used. The domain was divided into corridors bounded by the midlines between transects and into three bands by
distance from each transect's seaward fringe (0 to 60, 60 to 120 and over 120 m). For each band the recorded trees
define a reference belt whose density, species composition, size marks and spatial statistics (pair correlation and mark
correlation functions) are the target. Each corridor by band tile was then filled with the multi-trait point
pattern reconstruction algorithm of Wudel et al. (2023), which places synthetic trees with species, stem radius and
height as marks until the pair correlation function and the mark correlation functions of species, stem radius
and height of the tile match those of the reference belt, with the measured trees held fixed; trees falling on the treeless band seaward of the forest edge
were removed. This gave five replicate stands of about 2,100 trees. The representative member reproduces the band
densities within 2 %, the composition within one percentage point and the DBH and height distributions (KS p 0.28 to
1.0), and it is the parent inventory of every model run.

Limitation to state: the stand density follows the recorded nearest-tree reference (about 220 to 500 trees per ha),
two to four times lower in the landward bands than the Cottam-Curtis distance estimator from the same records, and
fine-scale spatial structure is not claimed.

Reference: Wudel, C., Schlicht, R. & Berger, U. (2023) Multi-trait point pattern reconstruction of plant ecosystems.
Methods in Ecology and Evolution 14, 2668-2679. (MTPPR v2.0.0 R code used unmodified; full account in
reconstruction/METHODS_reconstruction.md.)

---

## Short methods section: field data and stand reconstruction (2026-11-09; supersedes the reconstruction-only paragraph above)

### Field data

The field data come from the mangrove survey of Gazi Bay, Kenya, by Dahdouh-Guebas et al. (2004). Five parallel
land-to-sea transects (T1 to T4 and T6, about 73 m apart, 60 to 170 m long) and one orthogonal transect (T5) were
laid across the forest from the landward edge to the seaward Sonneratia fringe. Adults were sampled with the
point-centred quarter method (PCQM): at sampling points every 10 m along the tape, the nearest adult tree in each of
the four quarters was recorded with its distance to the point, species, DBH, height and crown diameter, giving a
single adult census of 266 trees on the five parallel transects (R. mucronata 142, S. alba 56, C. tagal 34,
A. marina 30, B. gymnorhiza 4) and 104 more on T5. Juveniles were surveyed twice, in 1997 and 1999, in a quadrat
in each of the four quarters around the same sampling points, and recorded as densities per m2 by species (the
quadrat area is not documented in the paper, so densities, not counts, are the observations). Transects T1 to T5
hold 335 juvenile quadrats, of which 309 lie inside the model domain; T6 was used for the adult stand only. The
adult records provide the parent trees of the model, the 1999 juvenile densities the calibration target, and the
1997 densities an out-of-sample comparison.

### Stand reconstruction

The parent stand of the model was reconstructed from the adult PCQM records. The 195 recorded trees of T1, T3, T4
and T6 were placed at their measured distance from the tape at a random angle within their quarter; T2 was withheld
as a blind test and T5 was never used. The domain was divided into corridors bounded by the midlines between
transects and into three bands by distance from each transect's seaward fringe (0 to 60, 60 to 120 and over
120 m). For each band the recorded trees define a reference belt whose density, species composition, size marks and
spatial statistics (pair correlation and mark correlation functions) are the target. Each corridor by band tile was then
filled with the multi-trait point pattern reconstruction algorithm of Wudel et al. (2023), which places synthetic
trees with species, stem radius and height as marks until the pair correlation function and the mark correlation
functions of species, stem radius and height of the tile match those of the reference belt, with the measured trees held fixed; trees falling on
the treeless band seaward of the forest edge were removed. This gave five replicate stands of about 2,100 trees.
The representative member reproduces the band densities within 2 %, the composition within one percentage point
and the DBH and height distributions (KS p 0.28 to 1.0), and it is the parent inventory of every model run.

Limitation to state: the stand density follows the recorded nearest-tree reference (about 220 to 500 trees per ha),
two to four times lower in the landward bands than the Cottam-Curtis distance estimator from the same records, and
fine-scale spatial structure is not claimed.

References: Dahdouh-Guebas, F. et al. (2004) Marine Ecology Progress Series (the Gazi transect study; complete the
reference from the paper). Wudel, C., Schlicht, R. & Berger, U. (2023) Multi-trait point pattern reconstruction of
plant ecosystems. Methods in Ecology and Evolution 14, 2668-2679.

---

## Notes for the authors (not manuscript text)

- The appendix drag figure is the RECRUIT version (`fig_drag_zonation.png`). The stranded version
  (`fig_drag_zonation_stranded.png`) isolates transport; if it replaces the recruit version, change 'established
  recruits' to 'stranded propagules' in Appendix A.
- The equal-contrast figure uses the standard-supply runs and the three calibrated species only; no amplified runs.
- 20 replicate runs = seeds 42 to 175 in steps of 7; three seeds = 42, 49, 56.
- The draft table the user supplied (17 levels, '0-75 %', three species, draft mean 0.10) was the pre-final version;
  the table above replaces it. `results/equal_cv/table_trait_levels.md` (regenerated today) is the long version.
- Sources for the mobility and lateral-spread values (Di Nitto et al.; Van der Stocken et al. 2015) remain unverified.

Correction 2026-11-09: the MTPPR objective is the L1 distance over the pair correlation function and the mark
correlation functions (species, stem radius, height) plus mark products (METHODS_reconstruction.md section 6);
mingling and nearest-neighbour functions were validation checks only. Earlier wording "species mingling" as a
reconstruction target in this document is superseded.
