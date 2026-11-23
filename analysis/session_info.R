# Which R and which packages did the code use, and at which versions?
#
# Finds every package the simulation and analysis scripts load, by library()
# or require(), or call with pkg::, loads each one, and records its version.
# This is the list behind the manuscript's software section (Supporting
# Information): run it on the machine that ran the simulations.
#
# Writes
#   docs/session_info.csv       one row per package: name, version, base or not
#   docs/session_info.txt       sessionInfo() with every one of those packages
#                               loaded, for the full record (platform, BLAS, ...)
#   docs/package_citations.bib  citation() for R and each package besides base
#                               R, keyed by package name ("R" for R itself;
#                               a second citation of a package gets ".2", ...)

library(here)

# ---------------------------------------------------------------------------
# 1. The packages the code uses
# ---------------------------------------------------------------------------

scripts <- list.files(here(c("simulations", "analysis")), pattern = "\\.R$",
                      full.names = TRUE)
code <- unlist(lapply(scripts, readLines, warn = FALSE))
code <- sub("#.*$", "", code)   # ignore comments

loaded <- regmatches(code, regexpr("(library|require|requireNamespace)\\(\\s*[\"']?[A-Za-z][A-Za-z0-9.]*",
                                   code))
loaded <- sub(".*\\(\\s*[\"']?", "", loaded)
called <- unlist(regmatches(code, gregexpr("\\b[A-Za-z][A-Za-z0-9.]*(?=::)", code,
                                           perl = TRUE)))
packages <- sort(unique(c(loaded, called)))

# ---------------------------------------------------------------------------
# 2. Their versions
# ---------------------------------------------------------------------------

base_packages <- rownames(installed.packages(priority = "base"))
info <- data.frame(
  package = packages,
  version = vapply(packages, function(p) as.character(packageVersion(p)), ""),
  base    = packages %in% base_packages,
  row.names = NULL
)

dir.create(here("docs"), showWarnings = FALSE)
write.csv(info, here("docs", "session_info.csv"), row.names = FALSE)

invisible(lapply(packages, requireNamespace, quietly = TRUE))
writeLines(capture.output(sessionInfo()), here("docs", "session_info.txt"))

message(R.version.string, "; ", sum(!info$base), " packages besides base R:")
print(info[!info$base, c("package", "version")], row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Their citations, for the supplement
# ---------------------------------------------------------------------------

# The BibTeX entries of one citation() result, each keyed by the package name
bibtex_for <- function(package) {
  entries <- if (package == "R") citation() else citation(package)
  keys <- package
  if (length(entries) > 1) keys <- c(keys, paste0(package, ".", 2:length(entries)))
  vapply(seq_along(entries), function(i) {
    lines <- toBibtex(entries[i])
    lines[1] <- sub("\\{[^,]*,", paste0("{", keys[i], ","), lines[1])
    paste(lines, collapse = "\n")
  }, "")
}

cited <- c("R", info$package[!info$base])
bib <- unlist(lapply(cited, bibtex_for))
writeLines(paste(bib, collapse = "\n\n"), here("docs", "package_citations.bib"),
           useBytes = TRUE)
first_lines <- vapply(strsplit(bib, "\n"), `[`, "", 1)
message("Citation keys: ",
        paste(sub("^@[A-Za-z]+\\{([^,]*),.*$", "\\1", first_lines), collapse = ", "))
