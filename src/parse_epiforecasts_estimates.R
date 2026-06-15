#!/usr/bin/env Rscript
#
# Parse epiforecasts/BVDOutbreakSize release estimates into hubverse
# submission files for the BVBD Modeling Hub.
#
# The upstream repository (https://github.com/epiforecasts/BVDOutbreakSize)
# publishes a series of tagged `results-v*` releases. Each release attaches a
# `posterior_draws.csv` (the posterior sample for that vintage) and an
# `observations.toml` (recording the data cut-off in `as_of_date`). The
# headline cumulative outbreak-size posterior lives in one column of the draws
# file: `cumulative_cases` for the early "integral" model, and `C_T` for the
# later discrete-time "renewal" model.
#
# IMPORTANT - target semantics. The hub target ("cumulative cases") is the
# cumulative number of *symptomatic* cases, explicitly NOT underlying
# infections and NOT confirmed-only counts (see hub-config/tasks.json).
#   * integral `cumulative_cases` IS cumulative symptomatic cases: the
#     integral posterior_summary.csv reports `cumulative_cases` and
#     `cumulative_infections` as distinct quantities, with cases < infections.
#     This maps cleanly onto the hub target.
#   * renewal `C_T` is `cumsum(infections)[n]`, i.e. cumulative *infections*
#     (see src/models/priors.jl in the upstream repo). The renewal release
#     assets publish no symptomatic-case quantity, so the renewal submission
#     reports infections as a documented stand-in. This deviation from the
#     symptomatic-case target is recorded in the renewal model metadata.
#
# This script downloads the relevant release assets, computes the quantiles,
# median and mean required by the hub directly from the posterior draws, and
# writes one hubverse model-output CSV per (reference_date, model).
#
# Two hub models are produced:
#   * epiforecasts-integral  (vintages prior to the renewal model)
#   * epiforecasts-renewal   (the renewal-model vintage)
#
# Usage (from the hub root):
#   Rscript src/parse_epiforecasts_estimates.R
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

## The releases to parse. One row per submission file.
##   tag         - upstream release tag carrying the posterior_draws.csv asset
##   model       - hub model name (directory suffix); also selects the draws column
##   draws_col   - column of posterior_draws.csv holding the headline posterior
##   quantity    - what draws_col represents, relative to the hub target
##
## For 2026-05-18 the upstream repo published two integral vintages
## (results-v1.0.0 and results-v1.1.0); we use the later revision v1.1.0.
RELEASES <- data.frame(
  tag = c("results-v1.1.0", "results-v1.2.0", "results-v1.3.0", "results-v1.4.0"),
  model = c("integral", "integral", "integral", "renewal"),
  draws_col = c("cumulative_cases", "cumulative_cases", "cumulative_cases", "C_T"),
  quantity = c("symptomatic cases", "symptomatic cases", "symptomatic cases",
               "infections (stand-in for symptomatic cases)"),
  stringsAsFactors = FALSE
)

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

tmp <- tempfile("bvd_releases_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

for (i in seq_len(nrow(RELEASES))) {
  tag <- RELEASES$tag[i]
  model <- RELEASES$model[i]
  draws_col <- RELEASES$draws_col[i]
  model_id <- paste0(TEAM_ABBR, "-", model)

  message(sprintf("Processing %s (%s) ...", tag, model_id))
  rel_dir <- file.path(tmp, tag)
  dir.create(rel_dir, showWarnings = FALSE)

  obs_path <- fetch_asset(tag, "observations.toml", rel_dir)
  draws_path <- fetch_asset(tag, "posterior_draws.csv", rel_dir)

  reference_date <- read_as_of_date(obs_path)

  draws_df <- utils::read.csv(draws_path, check.names = FALSE)
  if (!draws_col %in% names(draws_df)) {
    stop(sprintf("column '%s' not found in posterior_draws.csv for %s",
                 draws_col, tag))
  }
  draws <- draws_df[[draws_col]]

  submission <- build_submission(reference_date, draws)

  out_dir <- file.path(hub_root, "model-output", model_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir,
                        sprintf("%s-%s.csv", reference_date, model_id))
  utils::write.csv(submission, out_file, row.names = FALSE, quote = FALSE)

  message(sprintf("  -> wrote %s (%d draws, median %s, quantity: %s)",
                  out_file, length(draws), round(median(draws)),
                  RELEASES$quantity[i]))
}

message("Done.")
