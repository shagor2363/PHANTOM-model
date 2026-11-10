################################################################################
# One (tile, seed) MTPPR v2.0.0 reconstruction job for the Gazi Bay stand.
#
# Usage:
#   Rscript --vanilla reconstruction/mtppr_tile.R --tile=203 --seed=1 \
#       --geom=results/reconstruction/geom --out=results/reconstruction/mtppr
#       [--max-steps=400000 --min-steps=40000 --steps-per-point=200 --conv-start=5000
#        --conv-window=5000 --conv-tol=2e-6 --trace-every=500 --rmax=8 --dr=0.25
#        --bw=0.5 --pin-n=0 --lambda-target=0 --engine=PATH]
#
# The engine file (MTPPR-v.2.0.0.R) is sourced unmodified.  Two hooks are supplied
# through its own interfaces: `in_window()` is overridden (global lookup) so the
# reference core can be a union of rotated belt rectangles, and `edge_correction` is
# a function argument.  Reference for a tile = every reference point of the tile's
# band (all transects); fixed points = every reference point (any band) inside the
# tile box, plus the exclusion sectors of all transects as forbidden zones.
#
# 2026-09-20 B0 ADAPTATION "count control" (--pin-n=1 --lambda-target=<trees per m2>):
#   the engine's initial pattern is Poisson(n_core * |wnd| / |core_wnd box|) uniform in the
#   buffered window (engine lines 356-359); the declared core_wnd box enters the engine ONLY
#   there (the reference subset is selected by the overridden in_window -> belt rectangles,
#   engine line 351; the reference edge correction is edge_correction_union(r, core_wnd) ->
#   belt rectangles, engine line 327).  Declaring |box| = n_core / lambda_target therefore
#   makes the initial count lambda_target * |wnd| (the buffered window area) while the
#   reference statistics are unchanged; add/delete proposals are switched off, so the count
#   stays at its initial value.  Under count control the fixed (forbidden) region is the set
#   of reference points themselves plus the exclusion sectors, NOT the whole reference belts
#   (the belts hold only the recorded nearest-tree sample, far below lambda_target).
#   Test switches used by b0_check.py only: --core-area-factor (multiplies the declared box
#   area), --fixed-all (every point fixed: the initial pattern is the reference itself, so
#   e_initial is the engine's own reference-vs-reference energy), --dump-ref (writes the
#   reference statistics that depend on core_wnd), --ref-density-scale (candidate remedy,
#   default off: multiplies the reference per-point moments by lambda_target / lambda_belt,
#   the p-thinning-consistent normalisation, through the edge_correction argument).
################################################################################

args <- commandArgs(trailingOnly = TRUE)
opt <- list(tile = NA, seed = NA, geom = NA, out = NA, engine = NA,
            max_steps = 400000, min_steps = 40000, steps_per_point = 200,
            conv_start = 5000, conv_window = 5000, conv_tol = 2e-6,
            trace_every = 500, rmax = 8, dr = 0.25, bw = 0.5, pin_n = 0, quiet = 0, margin = 8,
            lambda_target = 0, core_area_factor = 1, fixed_all = 0, dump_ref = 0, ref_density_scale = 0,
            exact_n0 = 1, stem_suffix = "", copy = 0, fixed_mode = "belts")
# 2026-09-20 B0 SUPERPOSITION candidate (superpose.py): --copy=c (c >= 1) gives the copy an independent
# RNG stream (rseed = seed*1000 + tile + 7919*(c-1)) and the stem suffix _c{c}; --fixed-mode=belts is
# the standard run (measured trees fixed, reference belts and exclusion sectors forbidden; copy 1),
# --fixed-mode=sectors_only drops the measured trees as fixed points (they are NOT copied into the
# pattern; only the exclusion sectors stay forbidden; copies 2..ceil(k)).  The reference pattern for
# the energy is the same in every copy.
for (a in args) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  key <- gsub("-", "_", kv[1], fixed = TRUE)
  opt[[key]] <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "1"
}
num_keys <- c("tile", "seed", "max_steps", "min_steps", "steps_per_point", "conv_start",
              "conv_window", "conv_tol", "trace_every", "rmax", "dr", "bw", "pin_n", "quiet", "margin",
              "lambda_target", "core_area_factor", "fixed_all", "dump_ref", "ref_density_scale", "exact_n0", "copy")
for (k in num_keys) opt[[k]] <- as.numeric(opt[[k]])
stopifnot(!is.na(opt$tile), !is.na(opt$seed), !is.na(opt$geom), !is.na(opt$out))
tile_id <- as.integer(opt$tile)
seed <- as.integer(opt$seed)
copy_id <- as.integer(opt$copy)
stopifnot(opt$fixed_mode %in% c("belts", "sectors_only"))
if (copy_id > 0 && opt$stem_suffix == "") opt$stem_suffix <- sprintf("_c%d", copy_id)
rseed <- seed * 1000 + tile_id + (if (copy_id > 1) 7919 * (copy_id - 1) else 0)

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(sub("^--file=", "", fa[1])) else getwd()
})
engine <- if (!is.na(opt$engine)) opt$engine else file.path(
  script_dir, "..", "ChrisWudel-Multi-trait-point-pattern-reconstruction-6982c39",
  "func", "MTPPR-v.2.0.0.R")
stopifnot(file.exists(engine))
source(engine)
stopifnot(identical(environment(construct), globalenv()))

SPECIES_LEVELS <- c("Avicennia", "Bruguiera", "Ceriops", "Rhizophora", "Sonneratia")

# --- Rotated rectangles and sectors ---------------------------------------------------------

in_rects <- function(pts, R) {
  if (is.null(R) || nrow(R) == 0) return(rep(FALSE, nrow(pts)))
  n <- nrow(pts); m <- nrow(R)
  dx <- outer(pts[, 1], R$x0, "-")
  dy <- outer(pts[, 2], R$y0, "-")
  along <- sweep(dx, 2, R$ux, "*") + sweep(dy, 2, R$uy, "*")
  perp  <- sweep(dx, 2, R$uy, "*") - sweep(dy, 2, R$ux, "*")
  hit <- along >= matrix(R$s0, n, m, byrow = TRUE) & along < matrix(R$s1, n, m, byrow = TRUE) &
         abs(perp) < matrix(R$hw, n, m, byrow = TRUE)
  if (!is.null(R$cx0)) {
    cx0 <- ifelse(is.finite(R$cx0), R$cx0, -1e9); cx1 <- ifelse(is.finite(R$cx1), R$cx1, 1e9)
    cy0 <- ifelse(is.finite(R$cy0), R$cy0, -1e9); cy1 <- ifelse(is.finite(R$cy1), R$cy1, 1e9)
    hit <- hit & outer(pts[, 1], cx0, ">=") & outer(pts[, 1], cx1, "<") &
                 outer(pts[, 2], cy0, ">=") & outer(pts[, 2], cy1, "<")
  }
  rowSums(hit) > 0
}

.engine_in_window <- in_window
in_window <- function(points, wnd) {
  if (is.null(wnd$rects)) return(.engine_in_window(points, wnd))
  in_rects(points, wnd$rects)
}

# Ohser-Stoyan factor |W| / gamma_W(r) for a union of rectangles (rotation-invariant)
edge_correction_union <- function(r, wnd) {
  if (is.null(wnd$rects)) return(edge_correction_OhserStoyan(r, wnd))
  R <- wnd$rects
  A <- R$area                     # clipped pieces: rectangle of equal area and belt width
  denom <- 0
  for (i in seq_len(nrow(R))) {
    f <- edge_correction_OhserStoyan(r, list(min = c(0, 0), max = c(2 * R$hw[i], A[i] / (2 * R$hw[i]))))
    denom <- denom + A[i] / f
  }
  sum(A) / denom
}

# 2026-09-20 B0 -- REFERENCE MOMENT RESCALING (count-control adaptation; NOT a geometric edge
# correction).  The engine's energy compares per-point second-order moments f/n, which estimate
# lambda * g(r) (and lambda * k_mm(r) etc.), so they scale with the density of the pattern they
# are computed on.  Under count control the new pattern is generated at lambda_target while the
# reference is the recorded nearest-tree sample in its belts (lambda_belt << lambda_target).  The
# reference per-point moments are therefore multiplied by lambda_target / lambda_belt.  This
# treats the recorded sample AS IF it were an independent p-thinning of the stand (a p-thinning
# preserves g(r) and scales the per-point neighbour moments by p).  That is an APPROXIMATION:
# the actual PCQM procedure selects the nearest tree per quarter, so the selection depends on the
# tree's position relative to the sample point and on its neighbours, and the behaviour of the
# sample's second-order structure depends on the underlying pattern (assessed on synthetic stands
# in b0_check/synthetic).  It is applied through the engine's `edge_correction` argument (the only
# multiplicative hook on the reference moments), but it is recorded and reported separately from
# edge_correction_union.
reference_moment_rescaling <- function(wnd) {
  if (is.null(wnd$ref_moment_rescale)) 1 else wnd$ref_moment_rescale
}
edge_correction_with_rescaling <- function(r, wnd) edge_correction_union(r, wnd) * reference_moment_rescaling(wnd)

# --- Inputs ------------------------------------------------------------------------------------

tiles <- read.csv(file.path(opt$geom, "tiles.csv"))
core_rects <- read.csv(file.path(opt$geom, "core_rects.csv"))
for (cc in c("cx0", "cx1", "cy0", "cy1")) if (!is.null(core_rects[[cc]])) core_rects[[cc]] <- as.numeric(core_rects[[cc]])
reference <- read.csv(file.path(opt$geom, "reference.csv"), stringsAsFactors = FALSE)
sectors <- read.csv(file.path(opt$geom, "sectors.csv"))

tile <- tiles[tiles$tile_id == tile_id, ]
stopifnot(nrow(tile) == 1)
band <- tile$band
transect <- tile$transect
band_rects <- core_rects[core_rects$band == band, ]
box <- list(min = c(tile$xp_min, tile$yp_min), max = c(tile$xp_max, tile$yp_max))
wnd <- list(min = box$min - opt$margin, max = box$max + opt$margin)   # buffered reconstruction window
in_box <- function(pts, w) .engine_in_window(pts, w)

ref_band <- reference[reference$band == band, ]
ref_pts <- cbind(reference$xp, reference$yp)
extra <- reference[reference$band != band & in_box(ref_pts, wnd), ]
ref <- rbind(ref_band, extra)
stopifnot(nrow(ref_band) > 0, nrow(band_rects) > 0)
core_area <- tile$core_area_union

sec_pts <- cbind(sectors$px, sectors$py)
sec <- sectors[sqrt(pmin((sec_pts[, 1] - wnd$min[1])^2, (sec_pts[, 1] - wnd$max[1])^2) *
                   (sec_pts[, 1] < wnd$min[1] | sec_pts[, 1] > wnd$max[1]) +
                   pmin((sec_pts[, 2] - wnd$min[2])^2, (sec_pts[, 2] - wnd$max[2])^2) *
                   (sec_pts[, 2] < wnd$min[2] | sec_pts[, 2] > wnd$max[2])) < 25, ]
in_sector <- function(pts) {
  if (!nrow(sec)) return(rep(FALSE, nrow(pts)))
  dx <- outer(pts[, 1], sec$px, "-")
  dy <- outer(pts[, 2], sec$py, "-")
  along <- sweep(dx, 2, sec$ux, "*") + sweep(dy, 2, sec$uy, "*")
  perp  <- sweep(dx, 2, sec$uy, "*") - sweep(dy, 2, sec$ux, "*")
  hit <- sweep(perp, 2, sec$sx, "*") > 0 & sweep(along, 2, sec$sy, "*") > 0 &
    (dx^2 + dy^2) < matrix(sec$radius^2, nrow(pts), nrow(sec), byrow = TRUE)
  rowSums(hit) > 0
}
# 2026-09-20 B0: count control = pinned count at a target density (see header)
count_control <- (opt$pin_n > 0 && opt$lambda_target > 0)
ref_key <- paste(ref$xp, ref$yp)
is_ref <- function(pts) !is.na(match(paste(pts[, 1], pts[, 2]), ref_key))
n_core <- sum(in_rects(cbind(ref_band$xp, ref_band$yp), band_rects))   # reference points the engine will select
fixed_fn_belts <- function(p) {                 # original rule: whole reference belts forbidden
  pts <- p$points
  inside <- in_box(pts, wnd)
  out <- inside & in_rects(pts, core_rects)
  if (nrow(sec) > 0) out <- out | (inside & in_sector(pts))
  out
}
fixed_fn_sectors <- function(p) {               # count control: reference points + sectors only
  pts <- p$points
  inside <- in_box(pts, wnd)
  out <- inside & is_ref(pts)
  if (nrow(sec) > 0) out <- out | (inside & in_sector(pts))
  out
}
fixed_fn_sectors_only <- function(p) {           # superposition copies >= 2: measured trees NOT fixed, sectors forbidden
  pts <- p$points
  inside <- in_box(pts, wnd)
  out <- if (nrow(sec) > 0) inside & in_sector(pts) else rep(FALSE, nrow(pts))
  out & !is_ref(pts)   # a reference point lying inside another quarter's sector (placement conflict) must not be copied
}
fixed_fn <- if (opt$fixed_all > 0) function(p) rep(TRUE, nrow(p$points)) else
            if (opt$fixed_mode == "sectors_only") fixed_fn_sectors_only else
            if (count_control) fixed_fn_sectors else fixed_fn_belts

set.seed(rseed)
mc <- cbind(runif(40000, wnd$min[1], wnd$max[1]), runif(40000, wnd$min[2], wnd$max[2]))
mc_fixed <- (if (opt$fixed_all > 0) fixed_fn_belts else fixed_fn)(list(points = mc))
mc_tile <- in_box(mc, box)
free_frac <- 1 - mean(mc_fixed)
wnd_area <- prod(wnd$max - wnd$min)
tile_free_area <- wnd_area * mean(mc_tile & !mc_fixed)          # tile area outside the forbidden sectors
ref_pts_all <- cbind(ref$xp, ref$yp)
n_fixed_tile <- sum(fixed_fn(list(points = ref_pts_all)) & in_box(ref_pts_all, box))   # recorded trees fixed inside the tile
if (count_control) {
  # 2026-09-20 B0 (after the first check: explicit initial count, tile-targeted).  New points are
  # uniform over the buffered window at density rho; in the tile they add to the fixed recorded
  # trees, so rho = (lambda_target * |tile| - n_fixed_tile) / |tile free area| and the initial
  # count over the whole buffered window is N0 = rho * |wnd|.  The declared box area is then
  # n_core * |wnd| / N0 (engine lines 356-357) and, with --exact-n0=1 (default), the engine's
  # rpois(1, N0) is replaced by round(N0) through a global override of rpois (engine untouched:
  # construct() lives in the global environment, so the override is found first; rpois is used
  # nowhere else in the engine), removing the Poisson noise of the pinned count.
  rho <- (opt$lambda_target * tile$area - n_fixed_tile) / tile_free_area
  n0_target <- rho * wnd_area
  virtual_area <- n_core * wnd_area / n0_target
} else if (opt$pin_n > 0) {
  virtual_area <- core_area / max(free_frac, 1e-6)   # legacy pin-n path (unused)
} else {
  virtual_area <- core_area
}
virtual_area <- virtual_area * opt$core_area_factor    # test switch only (b0_check), 1 otherwise
core_wnd <- list(min = c(0, 0), max = rep(sqrt(virtual_area), 2), rects = band_rects)
lambda_belt <- n_core / core_area
rescale_applied <- (count_control && opt$ref_density_scale > 0)
rescale_factor <- if (rescale_applied) opt$lambda_target / lambda_belt else 1
if (rescale_applied) core_wnd$ref_moment_rescale <- rescale_factor   # reference moment rescaling (see above)
n0_expected <- n_core * wnd_area / virtual_area
if (count_control && opt$exact_n0 > 0) {
  .engine_rpois <- stats::rpois
  rpois <- function(n, lambda) rep.int(as.integer(round(lambda)), n)   # exact initial count (see above)
}

marked_pattern <- list(
  points = cbind(xp = ref$xp, yp = ref$yp),
  marks = list(r_stem = ref$r_stem_s, h_stem = ref$h_stem_s,
               species = factor(ref$species, levels = SPECIES_LEVELS)),
  core_wnd = core_wnd, wnd = wnd)

probabilities <- c(move_point = 0.3, exchange_coords = 0.1, exchange_mark = 0.2,
                   draw_mark = 0.2, delete_point = 0.1, add_point = 0.1)
if (opt$pin_n > 0) {
  probabilities[c("delete_point", "add_point")] <- 0
  probabilities <- probabilities / sum(probabilities)
}

if (opt$dump_ref > 0) {
  # reference statistics that could depend on the declared core_wnd (b0_check (ii))
  rr <- seq(opt$dr, opt$rmax, by = opt$dr)
  refstats <- list(tile_id = tile_id, seed = seed, virtual_area = virtual_area, box_side = sqrt(virtual_area),
                   core_area_factor = opt$core_area_factor, n_core = n_core,
                   n_in_window_core = sum(in_window(marked_pattern$points, core_wnd)),
                   wedge_core = edge_correction_union(rr, core_wnd), wedge_wnd = edge_correction_union(rr, wnd),
                   ref_moment_rescale = rescale_factor,
                   n_fixed_pattern = sum(fixed_fn(marked_pattern)), wnd_area = wnd_area,
                   n0_expected = n0_expected, lambda_target = opt$lambda_target, lambda_belt = lambda_belt,
                   free_frac = free_frac)
  dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
  writeLines(jsonlite::toJSON(refstats, auto_unbox = TRUE, pretty = TRUE, digits = NA),
             file.path(opt$out, sprintf("tile%d_seed%d%s_refstats.json", tile_id, seed, opt$stem_suffix)))
}

# --- Trace / stopping closure -------------------------------------------------------------------

tr <- new.env()
tr$step <- integer(0); tr$e <- numeric(0); tr$n <- integer(0)
tr$e_prev <- NA; tr$reason <- "running"; tr$last <- NULL; tr$n0 <- NA
stop_fn <- function(p, step, stat, ...) {
  n <- nrow(p$points)
  if (step == 0) tr$n0 <- n                    # initial count (buffered window), b0_check (ii)
  tr$last <- list(step = step, e = stat$e, e0 = stat$e_initial, n = n,
                  attempt = stat$attempt, success = stat$success)
  if (step %% opt$trace_every == 0) {
    tr$step <- c(tr$step, step); tr$e <- c(tr$e, stat$e); tr$n <- c(tr$n, n)
    if (opt$quiet == 0 && step %% (opt$trace_every * 20) == 0) {
      cat(sprintf("tile %d seed %d step %d energy %.5g n %d\n", tile_id, seed, step, stat$e, n))
      flush.console()
    }
  }
  if (step >= opt$max_steps) { tr$reason <- "max_steps"; return(TRUE) }
  if (step >= max(opt$min_steps, opt$steps_per_point * n)) { tr$reason <- "budget"; return(TRUE) }
  if (step >= opt$conv_start && step %% opt$conv_window == 0) {
    if (!is.na(tr$e_prev) && (tr$e_prev - stat$e) / opt$conv_window <= opt$conv_tol * stat$e) {
      tr$reason <- "converged"; return(TRUE)
    }
    tr$e_prev <- stat$e
  }
  FALSE
}

# --- Run ------------------------------------------------------------------------------------------

set.seed(rseed)
t0 <- proc.time()[["elapsed"]]
res <- tryCatch(construct_pattern(
  marked_pattern,
  fixed = fixed_fn,
  r = seq(opt$dr, opt$rmax, by = opt$dr),
  kernel = select_kernel("epanechnikov", bw = opt$bw, divisor = "r"),
  edge_correction = edge_correction_with_rescaling,   # = edge_correction_union x reference moment rescaling (1 unless applied)
  probabilities = probabilities,
  progress_stop = stop_fn
), error = function(e) e)
elapsed <- proc.time()[["elapsed"]] - t0

dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
stem <- file.path(opt$out, sprintf("tile%d_seed%d%s", tile_id, seed, opt$stem_suffix))

write_json <- function(x, path) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"), path)
  } else {
    writeLines(paste0("{", paste(sprintf('"%s": %s', names(x),
      vapply(x, function(v) if (is.character(v)) sprintf('"%s"', v) else as.character(v), "")),
      collapse = ", "), "}"), path)
  }
}

if (inherits(res, "error")) {
  write_json(list(tile_id = tile_id, seed = seed, error = conditionMessage(res)), paste0(stem, "_summary.json"))
  stop(sprintf("tile %d seed %d failed: %s", tile_id, seed, conditionMessage(res)))
}

pts <- res$points
# 2026-09-20 B0: under count control generated points may lie inside the belts; fixed = reference points
is_fixed <- if (count_control || opt$fixed_all > 0 || opt$fixed_mode == "sectors_only") is_ref(pts) else in_rects(pts, core_rects)
in_tile <- .engine_in_window(pts, box)
out <- data.frame(xp = pts[, 1], yp = pts[, 2], r_stem_s = res$marks$r_stem, h_stem_s = res$marks$h_stem,
                  species = as.character(res$marks$species), is_fixed = as.integer(is_fixed),
                  in_tile = as.integer(in_tile))
write.csv(out, paste0(stem, "_points.csv"), row.names = FALSE)
write.csv(data.frame(step = tr$step, energy = tr$e, n = tr$n), paste0(stem, "_trace.csv"), row.names = FALSE)

last <- tr$last
summary <- list(
  tile_id = tile_id, seed = seed, rseed = rseed, copy = copy_id, fixed_mode = opt$fixed_mode, band = band, transect = transect,
  n_ref = nrow(ref_band), n_extra_fixed = nrow(extra), n_sectors = nrow(sec),
  n_fixed = sum(is_fixed & in_tile), n_final = sum(in_tile), n_final_buffered = nrow(pts),
  margin = opt$margin, expected_n = tile$expected_n,
  expected_n_eligible = if (!is.null(tile$expected_n_eligible)) tile$expected_n_eligible else NA,
  eligible_area = if (!is.null(tile$eligible_area)) tile$eligible_area else NA,
  tile_area = tile$area, core_area = core_area, free_frac = free_frac, virtual_area = virtual_area,
  # 2026-09-20 B0 count control record
  count_control = count_control, lambda_target = opt$lambda_target, n_core = n_core, lambda_belt = lambda_belt,
  n0_expected = n0_expected, n_initial = tr$n0, wnd_area = wnd_area, tile_free_area = tile_free_area,
  n_fixed_tile_ref = n_fixed_tile, exact_n0 = opt$exact_n0,
  core_area_factor = opt$core_area_factor, fixed_all = opt$fixed_all, ref_density_scale = opt$ref_density_scale,
  reference_moment_rescaling = list(
    applied = rescale_applied, factor = rescale_factor, lambda_target = opt$lambda_target, lambda_belt = lambda_belt,
    description = "reference moment rescaling (count-control adaptation): reference per-point second-order moments x lambda_target / lambda_belt; independent-thinning APPROXIMATION of the nearest-tree sample, not a geometric edge correction"),
  e_initial = last$e0, e_final = last$e, steps = last$step, stop_reason = tr$reason,
  seconds = elapsed, pin_n = opt$pin_n, rmax = opt$rmax, dr = opt$dr, bw = opt$bw,
  max_steps = opt$max_steps, min_steps = opt$min_steps, steps_per_point = opt$steps_per_point,
  conv_tol = opt$conv_tol, conv_window = opt$conv_window)
for (a in names(last$attempt)) {
  summary[[paste0("attempt_", a)]] <- unname(last$attempt[a])
  summary[[paste0("success_", a)]] <- unname(last$success[a])
}
write_json(summary, paste0(stem, "_summary.json"))
cat(sprintf("tile %d seed %d done: n %d in tile (expected %.0f), energy %.4g -> %.4g, %d steps (%s), %.0f s\n",
            tile_id, seed, sum(in_tile), tile$expected_n, last$e0, last$e, last$step, tr$reason, elapsed))
