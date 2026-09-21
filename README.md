# stats

Working repository for a 45-minute talk to IDFG colleagues: *"uncertainty in
escapement estimation at Lower Granite: adults and smolts, one framework."*

## What this is

Material that reverse engineers the two production escapement tools at Lower
Granite Dam — `escapeLGD` (EASE, adults) and `smoltEASE` (SCRAPI2, smolts) —
treating them as one estimator structure applied at two life stages. There are
two tracks.

### Track 1 — the short talk track (`R/session*.R`)

Nine self-contained R sessions, each one idea, under a minute to run, one
figure. They exist to make a 45-minute talk explainable out loud. Most sessions
write a toy version by hand, then point at the exact production file and function
that does the same job; session 7 instead runs production `SCRAPI2()` on the real
MY2025 inputs and checks it by hand.

### Track 2 — the interactive course (`inst/tutorials/estimator-histories/`)

One cumulative `learnr` tutorial that rebuilds the compound bootstrap by hand
against a single known-truth simulated fishery, following the estimators from
SCOBI and SCRAPI to EASE and SCRAPI2. It is interactive: you edit and run the
estimator code, predict before revealing, and answer graded questions. Launch it
from the repository root:

```r
rmarkdown::run("inst/tutorials/estimator-histories/estimator-histories.Rmd")
```

`LEARNING_PATH.md` explains both tracks and the five course modules;
`docs/bootstrap_map.md` is the precise production reference table for what each
estimator resamples and where the draw enters.

## Where to start

1. Read `PLAN.md`. It states the one idea the talk is built on, the five-part
   session format, and the nine sessions in order.
2. Read `docs/session02_shared_skeleton.md` for the component-by-component
   comparison of the adult and smolt estimators. That table is the spine.
3. Work the sessions in order in `R/`. Run each from the repository root
   (`source("R/session01_count_expansion.R")`); each prints its recovery, writes
   a figure to `figs/`, and writes a plain-English explanation to
   `docs/sessionNN_explain.md`. Session 2 also writes the skeleton table to
   `docs/session02_shared_skeleton.md`. Base R throughout, with two exceptions:
   session 5 uses `glmmTMB` (with a `glm` fallback), and session 7 requires the
   `smoltEASE` package because it runs production `SCRAPI2()` on the real data.

## What is not here

Everything cut from the earlier construction-oriented plan is in `BACKLOG.md`,
each item with one line saying why it was cut. Nothing was deleted.

## Ground rules

These apply to the **talk track** only. Each session teaches one idea and one
idea only: simulate, build the estimator by hand, recover against truth, one
figure, Locate. A second figure means a second session. Single-use helpers are
inlined; no defensive input checking, since these are teaching scripts, not
package code. No session over a minute to run. No banner comments. No placeholder
headings. No dates — those get added by hand.

The interactive course is cumulative by design and does not follow the
one-idea/one-figure/one-minute rules; it preserves one simulated fishery across
all five modules.
