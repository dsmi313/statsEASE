# Smoke test for the estimator-histories learnr tutorial.
# Run from the repository root:  Rscript tests/smoke_test_tutorial.R
# Starts from a clean R session, verifies packages, structure, a full render,
# and the statistical invariants the tutorial relies on. Stops non-zero on any
# failure so it can gate CI.

rmd  <- "inst/tutorials/estimator-histories/estimator-histories.Rmd"
help <- "R/tutorial_helpers.R"
fail <- 0L
ok   <- function(cond, label) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) fail <<- fail + 1L
}

# 1. required packages ------------------------------------------------------
for (p in c("learnr", "shiny", "rmarkdown", "knitr"))
  ok(requireNamespace(p, quietly = TRUE), paste("package available:", p))
ok(file.exists(rmd),  "tutorial Rmd present")
ok(file.exists(help), "tutorial_helpers.R present")

# 2. chunk structure --------------------------------------------------------
lines <- readLines(rmd)
hdr   <- grep("^```\\{r", lines, value = TRUE)
labels <- sub("^```\\{r ([A-Za-z0-9_.-]+).*", "\\1", hdr)
labels <- labels[labels != hdr]                       # drop headerless chunks
ok(!any(duplicated(labels)),
   paste("no duplicate chunk labels",
         if (any(duplicated(labels)))
           paste("(dupes:", paste(unique(labels[duplicated(labels)]), collapse = ", "), ")")
         else ""))

ex_hdr <- grep("exercise=TRUE", hdr, value = TRUE)
ex_lab <- sub("^```\\{r ([A-Za-z0-9_.-]+).*", "\\1", ex_hdr)
for (e in ex_lab) {
  ok(paste0(e, "-solution") %in% labels, paste("solution exists for", e))
  ok(paste0(e, "-hint")     %in% labels, paste("hint exists for", e))
}
setup_refs <- regmatches(ex_hdr, regexpr('exercise.setup *= *"[^"]+"', ex_hdr))
setup_refs <- sub('.*"([^"]+)".*', "\\1", setup_refs)
for (s in setup_refs)
  ok(s %in% labels, paste("exercise.setup chunk exists:", s))

# 3. full render (detects chunk errors, runs non-exercise chunks) -----------
render_ok <- tryCatch({
  out <- tempfile(fileext = ".html")
  rmarkdown::render(rmd, output_file = out, quiet = TRUE,
                    envir = new.env())
  file.exists(out)
}, error = function(e) { cat("render error:", conditionMessage(e), "\n"); FALSE })
ok(render_ok, "tutorial renders without chunk errors")

# 4. statistical invariants -------------------------------------------------
source(help)
d <- simulate_smolt_season(seed = 6L, B = 300L)
ok(all(is.finite(unlist(d$truth))), "smolt truth values finite")
ok(all(d$ge_draws > 0 & d$ge_draws <= 1), "GE draws in (0, 1]")
ok(all(d$fish$SR > 0 & d$fish$SR <= 1),   "SR in (0, 1]")

# additivity of a full compound row: stockA + stockB == WildSmolts, exactly
one_row_additivity <- local({
  b <- 1L; w <- d$weeks
  pb <- w$rate * d$ge_draws[, b]
  cs <- rbinom(6, round(w$Tally / pb), pb)
  ps <- tapply(cs / pb, w$stratum, sum)[as.character(1:3)]
  Tb <- w$rate[d$rear$week] * d$ge_draws[cbind(d$rear$week, b)]
  wild <- 0; stockA <- 0
  for (h in 1:3) {
    ri  <- which(d$rear$stratum == h); idx <- sample(ri, length(ri), replace = TRUE)
    pW  <- sum((1 / Tb[idx]) * (d$rear$rear[idx] == "W")) / sum(1 / Tb[idx])
    wild_h <- pW * ps[h]; wild <- wild + wild_h
    fi  <- which(d$fish$stratum == h); fx <- sample(fi, length(fi), replace = TRUE)
    pa  <- sum((1 / d$fish$SR[fx]) * (d$gsi_draws[fx, b] == "A")) / sum(1 / d$fish$SR[fx])
    stockA <- stockA + pa * wild_h
  }
  abs((stockA + (wild - stockA)) - wild)
})
ok(one_row_additivity < 1e-6, "full-compound row additivity stockA+stockB==WildSmolts")

# adult truths use the correct estimands: escapement = ascensions*(1 - p_fall)
a <- simulate_adult_season(seed = 6L, B = 100L)
ok(all(is.finite(unlist(a$truth))), "adult truth values finite")
ok(a$truth$true_escapement < a$truth$true_ascensions,
   "true_escapement < true_ascensions (fallback removes fish)")
ok(abs(a$truth$true_escapement -
       round(a$truth$true_ascensions * (1 - a$truth$p_fall))) <= 1,
   "true_escapement == ascensions * (1 - p_fall)")
ok(abs((a$truth$escA + a$truth$escB) - a$truth$true_escapement) <= 1,
   "escA + escB == true_escapement")
ok(all(a$gsi_draws %in% c("A", "B")), "adult GSI draws are valid stock labels")

# 5. summary ----------------------------------------------------------------
if (fail > 0) stop(sprintf("smoke test FAILED: %d check(s) failed", fail))
cat("\nAll smoke-test checks passed.\n")
