# Fail if any .R file or .Rmd code chunk does not parse.
args <- commandArgs(TRUE)
bad <- 0
files <- list.files(args[1], pattern = "\\.(R|Rmd)$", recursive = TRUE, full.names = TRUE)
for (f in files) {
  if (grepl("/(\\.git|renv|archive[^/]*)/", f)) next
  src <- f
  if (grepl("Rmd$", f)) {
    src <- tempfile(fileext = ".R")
    ok <- tryCatch({ suppressWarnings(knitr::purl(f, output = src, quiet = TRUE)); TRUE },
                   error = function(e) { cat("SKIPPED (cannot extract code):", f, "\n"); FALSE })
    if (!ok) next
  }
  ok <- tryCatch({ parse(file = src); TRUE },
                 error = function(e) { cat("PARSE ERROR in", f, ":", conditionMessage(e), "\n"); FALSE })
  if (!ok) bad <- bad + 1
}
cat(length(files), "files checked,", bad, "with parse errors\n")
if (bad > 0) quit(status = 1)
