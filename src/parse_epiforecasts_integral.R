#!/usr/bin/env Rscript
#
# Parse epiforecasts/BVDOutbreakSize INTEGRAL-model release estimates into
# hubverse submission files for the BVBD Modeling Hub (model epiforecasts-integral).
#
# The upstream repository (https://github.com/epiforecasts/BVDOutbreakSize)
# published a series of tagged `results-v*` releases. The early "integral" model
# vintages are parsed here; the later "renewal" model is handled by a separate
# script (src/parse_epiforecasts_renewal.R).
#
# Each release attaches a `posterior_draws.csv` (the posterior sample for that
# vintage) and an `observations.toml` (recording the data cut-off in
# `as_of_date`). The integral model's headline cumulative outbreak-size
# posterior is the `cumulative_cases` column of the draws file.
#
# Target semantics. The hub target ("cumulative cases") is the cumulative
# number of *symptomatic* cases, explicitly NOT underlying infections and NOT
# confirmed-only counts (see hub-config/tasks.json). The integral model reports
# outbreak size as a single coarse exponential-growth quantity that the authors
# label `cumulative_cases`; this is the value intended to map onto the hub
# target. (The integral model does not separately resolve infections vs symptom
# onsets the way the later renewal model does.)
#
# This script downloads the relevant release assets, computes the quantiles,
# median and mean required by the hub directly from the posterior draws, and
# writes one hubverse model-output CSV per integral vintage.
#
# Usage (from the hub root):
#   Rscript src/parse_epiforecasts_integral.R
#
# Requires the `gh` CLI to be installed and authenticated, plus the base R
# packages `utils` and `stats` (no extra package dependencies).

suppressWarnings(suppressMessages({
  library(utils)
  library(stats)
}))

## ---------------------------------------------------------------------------
## Configuration
## ---------------------------------------------------------------------------

UPSTREAM_REPO <- "epiforecasts/BVDOutbreakSize"

## Hub-required output. The `quantile` output type must carry exactly these
## probability levels; `median` and `mean` are point estimates with an "NA"
## output_type_id. See hub-config/tasks.json.
QUANTILE_LEVELS <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
TARGET <- "cumulative cases"
LOCATION <- "CD"
TEAM_ABBR <- "epiforecasts"
MODEL_ABBR <- "integral"
MODEL_ID <- paste0(TEAM_ABBR, "-", MODEL_ABBR)

## The integral posterior-draws column carrying the headline outbreak-size
## posterior (cumulative symptomatic cases, per the target semantics above).
DRAWS_COL <- "cumulative_cases"

## The integral-model releases to parse. One row per submission file.
## For 2026-05-18 the upstream repo published two integral vintages
## (results-v1.0.0 and results-v1.1.0); we use the later revision v1.1.0.
RELEASE_TAGS <- c("results-v1.1.0", "results-v1.2.0", "results-v1.3.0")

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

## Download a single named asset from an upstream release into `dir`,
## returning its path. Stops with an informative error on failure.
fetch_asset <- function(tag, file, dir) {
  dest <- file.path(dir, file)
  status <- system2(
    "gh",
    c("release", "download", tag, "-R", UPSTREAM_REPO,
      "-p", file, "-O", dest, "--clobber"),
    stdout = FALSE, stderr = FALSE
  )
  if (status != 0 || !file.exists(dest)) {
    stop(sprintf("failed to download asset '%s' from release '%s'", file, tag))
  }
  dest
}

## Read the data cut-off date from an observations.toml file.
read_as_of_date <- function(toml_path) {
  lines <- readLines(toml_path, warn = FALSE)
  hit <- grep("^\\s*as_of_date\\s*=", lines, value = TRUE)
  if (length(hit) == 0) {
    stop(sprintf("no as_of_date found in %s", toml_path))
  }
  date_str <- sub('.*=\\s*"?([0-9]{4}-[0-9]{2}-[0-9]{2})"?.*', "\\1", hit[1])
  if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date_str)) {
    stop(sprintf("could not parse as_of_date from line: %s", hit[1]))
  }
  date_str
}

## Build the long-format hubverse table for one release.
build_submission <- function(reference_date, draws) {
  ## value column rounded to whole cases (the target counts people).
  q_values <- round(as.numeric(quantile(draws, probs = QUANTILE_LEVELS,
                                         names = FALSE, type = 7)))

  quantile_rows <- data.frame(
    reference_date = reference_date,
    target = TARGET,
    location = LOCATION,
    output_type = "quantile",
    output_type_id = as.character(QUANTILE_LEVELS),
    value = q_values,
    stringsAsFactors = FALSE
  )

  point_rows <- data.frame(
    reference_date = reference_date,
    target = TARGET,
    location = LOCATION,
    output_type = c("median", "mean"),
    output_type_id = c("NA", "NA"),
    value = c(round(median(draws)), round(mean(draws))),
    stringsAsFactors = FALSE
  )

  rbind(quantile_rows, point_rows)
}

## ---------------------------------------------------------------------------
## Main
## ---------------------------------------------------------------------------

## Resolve the hub root relative to this script so it can be run from anywhere.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else NA
hub_root <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), ".."))
} else {
  normalizePath(getwd())
}

tmp <- tempfile("bvd_integral_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

for (tag in RELEASE_TAGS) {
  message(sprintf("Processing %s (%s) ...", tag, MODEL_ID))
  rel_dir <- file.path(tmp, tag)
  dir.create(rel_dir, showWarnings = FALSE)

  obs_path <- fetch_asset(tag, "observations.toml", rel_dir)
  draws_path <- fetch_asset(tag, "posterior_draws.csv", rel_dir)

  reference_date <- read_as_of_date(obs_path)

  draws_df <- utils::read.csv(draws_path, check.names = FALSE)
  if (!DRAWS_COL %in% names(draws_df)) {
    stop(sprintf("column '%s' not found in posterior_draws.csv for %s",
                 DRAWS_COL, tag))
  }
  draws <- draws_df[[DRAWS_COL]]

  submission <- build_submission(reference_date, draws)

  out_dir <- file.path(hub_root, "model-output", MODEL_ID)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir,
                        sprintf("%s-%s.csv", reference_date, MODEL_ID))
  utils::write.csv(submission, out_file, row.names = FALSE, quote = FALSE)

  message(sprintf("  -> wrote %s (%d draws, median %s, quantity: symptomatic cases)",
                  out_file, length(draws), round(median(draws))))
}

message("Done.")
