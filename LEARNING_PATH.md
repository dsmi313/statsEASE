# Learning path

This repository has two tracks. They cover the same two production estimators
(EASE/`escapeLGD` for adults, SCRAPI2/`smoltEASE` for smolts) at different
depths.

## Track 1: the short talk track (`R/session*.R`)

Nine self-contained R scripts, each one idea, under a minute to run, one figure.
They exist to make a 45-minute talk explainable out loud. Start here for the
big picture. See `PLAN.md` and `README.md`; run each from the repository root:

```r
source("R/session01_count_expansion.R")
```

Each prints its recovery, writes a figure to `figs/`, and writes a plain-English
explanation to `docs/sessionNN_explain.md`. `docs/session02_shared_skeleton.md`
is the component-by-component comparison table and is the spine of the talk.

## Track 2: the interactive course (`inst/tutorials/estimator-histories/`)

One cumulative `learnr` tutorial that rebuilds the compound bootstrap by hand
against a single known-truth simulated fishery, following the estimators from
SCOBI and SCRAPI to EASE and SCRAPI2. Unlike the talk track, it is interactive:
you edit and run the estimator code yourself, predict before revealing, and
answer graded conceptual questions.

### How to launch it

The talk-track scripts in `R/` execute on load, so this project is not intended
to be installed as a package. Launch the tutorial directly from the repository
root:

```r
rmarkdown::run("inst/tutorials/estimator-histories/estimator-histories.Rmd")
```

or, equivalently, with learnr:

```r
learnr::run_tutorial(
  "inst/tutorials/estimator-histories",
  package = NULL,
  shiny_args = list(launch.browser = TRUE)
)
```

`Rscript run_estimator_histories.R` from the repository root does the same. If
the repository is later restructured into an installed package, the tutorial
would then also launch with
`learnr::run_tutorial("estimator-histories", package = "statsEASE")`.

Requires `learnr`, `shiny`, `rmarkdown`, and `knitr`. It uses only simulated
data with fixed seeds and known truth, and writes nothing to any production-data
directory.

### The five modules, in order

1. **The shared estimator.** Count, inclusion/detection probability, expanded
   abundance, composition, uncertainty propagation, in one place.
2. **Build a bootstrap from scratch.** Plug-in, then the parametric passage-count
   bootstrap, nonparametric fish resampling, one coherent GE posterior column
   per replicate, one coherent GSI column per replicate, up to the full compound
   SCRAPI2 bootstrap.
3. **SCRAPI to SCRAPI2.** Fixed vs uncertain GE and GSI; original SCRAPI's
   SR-weighted resampling vs SCRAPI2's uniform within-stratum resampling; the
   bootstrap-row mean vs the raw point estimate; percentile intervals and
   additivity.
4. **SCOBI to EASE.** Old SCOBI `wc_prop = 1`, EASE `wc_prop = 5/6`, nighttime
   passage, composition, GSI, fallback and reascension, the final interval.
5. **Model diagnostics and failure modes.** Low GE and right-skewed abundance,
   `E(1/GE)` vs `1/E(GE)`, poor PIT sample sizes, incorrect fish-resampling
   weights, biased GE, mixing GE columns across weeks, and a mean outside a
   percentile interval.

## Questions to ask at every step

- What is random in this step, and what am I conditioning on as observed?
- What statistical target does this operation estimate, and is the interval
  centered on it or on something else?
- If I widen the interval, does the problem go away? (Sampling noise: yes. A
  biased center or a cancelled weight: no.)
- Does this quantity add up to its parent total within every replicate?
- Which production function does this, and does it match what I just built?

## Suggested experiments

- Push the trap sample rate or the low-GE stratum harder and watch the right
  tail grow.
- Change the two stock sampling rates so they are equal, and confirm the
  SR-weighted and uniform resamples then agree (the cancellation only bites when
  sampling rate correlates with composition).
- Increase the PIT sample sizes on the adult side and watch the nighttime and
  fallback intervals tighten without moving their center.
- Bias the GE posterior and confirm no interval width restores coverage.

## Where the production claims come from

Every "Locate this in production" note and every row of `docs/bootstrap_map.md`
is traceable to the current `smoltEASE`, `escapeLGD`, or `EASE-code` source, not
to general knowledge. Key files: `smoltEASE/R/SCRAPI.R`,
`smoltEASE/R/SCRAPI2.R`, `smoltEASE/R/generate_ge_draws.R`,
`escapeLGD/R/night_fall_reascend_wc_binom.R`,
`escapeLGD/R/fallback_reascend_likelihood.R`,
`escapeLGD/R/apply_fallback_rates_summarize.R`, and the `EASE-code` "Step 3"
annual scripts for `wc_prop`.
