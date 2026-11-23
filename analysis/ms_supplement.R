# The manuscript's supplement: which figures are in it, and in what order.
#
# figures/supplement/ holds every supporting figure the analysis
# scripts draw, most of them for checking the results. This script copies the
# ones chosen for the manuscript's supplement into
# figures/ms_supplement/, numbered in the order they appear (S01,
# S02, ...). The table below is the list; edit it to add, drop, or reorder a
# figure, and rerun this script after the scripts that draw the figures.
#
# The folder is rebuilt from scratch each time, so it holds only the figures
# listed here.

library(dplyr)
library(readr)
library(tibble)
library(here)

# ---------------------------------------------------------------------------
# 1. The list, in the order of the manuscript
# ---------------------------------------------------------------------------

supp_dir <- here("figures", "supplement")
ms_dir   <- here("figures", "ms_supplement")

ms_figures <- tribble(
  ~file,                             ~from,    ~shows,
  # the synchronous counterparts of main-text Figs 2, 4, and 6 first, then
  # the rest in the order the Results cite them
  "synchronous_diversity.pdf",       supp_dir, "diversity with synchronous patches: dispersal replaces dormancy only where patches are asynchronous",
  "synchronous_stability.pdf",       supp_dir, "stability at four levels with synchronous patches",
  "synchronous_strategy_layers.pdf", supp_dir, "cued against bet-hedging with synchronous patches, active and dormant",
  "pathways_seedbank.pdf",           supp_dir, "each seed bank against none: metacommunity stability split into population, among species, and among patches",
  "pathways_strategy.pdf",           supp_dir, "cued against bet-hedging: the same split, active and dormant",
  "noise_strategy_layers.pdf",       supp_dir, "cued against bet-hedging under red, white, and blue noise, active and dormant",
  "sensitivity.pdf",                 supp_dir, "metacommunity stability as one parameter is varied: (A) each seed bank against none, (B) cued against bet-hedging",
  "combined_stability.pdf",          supp_dir, "active and dormant combined: stability at four levels for every strategy and no dormancy",
  "dormant_diversity.pdf",           supp_dir, "diversity of the dormant community: alpha, spatial beta, gamma"
)

# ---------------------------------------------------------------------------
# 2. Copy them, numbered
# ---------------------------------------------------------------------------

ms_figures <- ms_figures |>
  mutate(number = sprintf("S%02d", row_number()),
         source = file.path(from, file),
         target = file.path(ms_dir, paste0(number, "_", file)),
         found  = file.exists(source))

if (!all(ms_figures$found)) {
  stop("Not found: ", paste(ms_figures$file[!ms_figures$found], collapse = ", "))
}

dir.create(ms_dir, recursive = TRUE, showWarnings = FALSE)
unlink(Sys.glob(file.path(ms_dir, "*.pdf")))
invisible(file.copy(ms_figures$source, ms_figures$target))

# The list, with numbers, for the figure captions
ms_figures |>
  transmute(number, file, from = sub(paste0(here(), "/"), "", from), shows) |>
  write_csv(here("docs", "ms_supplement_figures.csv"))

message("Copied ", nrow(ms_figures), " figures to ", ms_dir)
