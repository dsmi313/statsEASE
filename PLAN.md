# Reverse engineering EASE and SCRAPI2

This repository has one purpose: to get me able to give a 45-minute talk to
IDFG colleagues, working title *"uncertainty in escapement estimation at Lower
Granite: adults and smolts, one framework."* Every session here exists to make
some part of that talk explainable out loud. If a session does not serve the
talk, it is not here; it is in `BACKLOG.md` with one line saying why.

## The one idea the talk is built on

EASE for adults (`escapeLGD`) and SCRAPI2 for smolts (`smoltEASE`) are not two
topics. They are one estimator structure applied at two life stages. The whole
talk is the claim that the same skeleton runs at both stages, that the pieces
line up almost row for row, and that the few places they diverge are the
interesting places. The generic skeleton, in the order the estimators apply it:

1. A total count at the dam (window count for adults, trap count for smolts).
2. A sampled fraction that the count is expanded by (counting-window open
   fraction for adults, trap sample rate for smolts).
3. A detection/observation expansion for fish the count structurally misses
   (nighttime passage for adults, guidance efficiency into the bypass for
   smolts).
4. Composition proportions that split the expanded total into groups (PBT and
   genetic stock, rear type, age).
5. A stage-specific adjustment (fallback and reascension for adults; wild
   fraction split for smolts).
6. Uncertainty propagated to an interval by resampling and posterior draws.

Both production tools are already built and in production. This plan is not a
path toward building them. It is reverse engineering: every session starts from
a toy version I write by hand, then points at the exact production file and
function that does the same job, and says how the two production
implementations handle it and where they part ways.

## Where every description of EASE comes from

Every description of EASE in this repository is derived from the `escapeLGD`
source, never from general knowledge and never from the package README. The
same holds for SCRAPI2 and the `smoltEASE` source. The "Locate" part of each
session names files and functions by their real names in those two repositories.

## Session format

Every session is one self-contained R script in `R/`, written MIT 18.05 style:
a stated objective, a worked example, and a short exercise. Each follows the
Kery and Schaub simulate-build-function loop. Five parts, in this order:

1. **Objective.** One sentence naming what I should be able to explain out loud
   afterward. Not what I should be able to code.
2. **Simulate.** Generate data from known truth. Truth is a variable at the top
   of the script.
3. **Build.** Write the estimator by hand. The likelihood is written out and
   optimized with `optim` or `nlminb`. No package does the estimating work.
4. **Recover.** Show the estimate against known truth, show the interval, plot
   it.
5. **Locate.** Name the exact file and function in **both** `escapeLGD` and
   `smoltEASE` that does this same job in production. In two sentences per side,
   say how each production version differs from the toy, and whether the two
   production implementations handle it the same way.

Each script also writes `docs/sessionNN_explain.md`: how I would explain that
one concept to a colleague in three minutes, in plain English, with no notation.

## The nine sessions

The sequence is centered on likelihood profiling. Sessions 1 and 2 lay the
generic skeleton; sessions 3 and 4 are the profiling core, with session 3 also
reading the four interval recipes off one rate; session 5 fits guidance
efficiency two ways and session 6 asks whether the composed interval covers what
it claims; session 7 grounds the whole machine in the real MY2025 run; session 8
is the sharpest adult/smolt divergence, composition by likelihood versus
accounting; session 9 assembles the talk.

### Session 1 — One count, expanded, by likelihood

**Objective:** explain why a count at a dam divided by the fraction of fish it
saw is a maximum likelihood estimate, and why that single idea is the base of
both the adult and the smolt estimator.

Binomial window count, likelihood written by hand, MLE recovered. This session
already existed under the old plan and is rewritten to the new format. The toy
estimator is written generically enough that both production estimators map
onto it: quantities are named generically (total count, sampled fraction,
composition proportions, expansion, uncertainty sources), not with
smolt-specific names.

**Locate:** `escapeLGD` `expand_wc_binom_night()` in
`R/night_fall_reascend_wc_binom.R`; `smoltEASE` `thetahat()` inside `SCRAPI2()`
in `R/SCRAPI2.R` (`dailypass <- Tally / Ptrue`).

### Session 2 — The shared skeleton

**Objective:** explain, component by component, how the adult and smolt
estimators are the same machine, and name the two or three places they are
genuinely not.

No simulation here. The deliverable is one table, derived from reading both
sources, with one row per estimator component. Columns: the generic component;
how `escapeLGD` implements it for adults; how `smoltEASE` implements it for
smolts; and whether the two differ structurally or only in parameterization.
Every uncertainty source each tool propagates is listed with its mechanism
(bootstrap, posterior draws, delta method, or not propagated at all), and any
source one propagates while the other does not is flagged. Written to
`docs/session02_shared_skeleton.md`. This table is the spine of the talk.

### Session 3 — Profile likelihood, and four intervals off one curve

**Objective:** explain how an interval falls out of the shape of the likelihood
curve for one parameter, without any normal approximation, and what each of the
four interval recipes means when all four are read off that same rate.

Plot the log-likelihood curve for one parameter, drop the line where twice the
log-likelihood falls by the chi-squared cutoff, and read the interval off the
curve. On the same rate, compute the other three intervals too — delta method,
bootstrap percentile, Bayesian credible — and show they nearly coincide here and
where they would part. This session absorbs the old standalone four-intervals
session; the profile curve stays as its figure.

**Locate:** the single-parameter binomial rate inside `escapeLGD` `nightFall()`
(nighttime passage, `p_night`) and the fallback rate `p_fa`; the per-rate
resampling in `smoltEASE` `SCRAPI2()`. The bootstrap percentile is the interval
both production tools actually ship (`quantile()` on the bootstrap matrix in
`SCRAPI2()` and in `apply_fallback_rates()`); the profile and delta intervals are
what the by-hand likelihood work implies; the Bayesian credible interval is what
`smoltEASE` `fit_ge_model()` produces for GE, the interval session 5 fits.

### Session 4 — Profile versus marginalize a nuisance parameter

**Objective:** explain the difference between profiling out a nuisance parameter
and integrating it out, and show they give nearly the same interval here and why.

Fit the two-parameter joint likelihood surface, then on one figure show the
profile trace for the parameter of interest against the marginalized curve.
The production anchor is exact: `escapeLGD`'s fallback likelihood is literally a
two-parameter binomial likelihood in P(fallback) and the nuisance P(reascend |
fallback). This is the one place the profiling core touches escapeLGD's own
two-parameter machinery, and it is one of the three genuine adult/smolt
divergences, since smolts have no reascension likelihood at all.

**Locate:** `escapeLGD` `fallback_log_likelihood()` and
`gradient_fallback_log_likelihood()` in `R/fallback_reascend_likelihood.R`,
optimized in `nightFall()`; `smoltEASE` has no two-parameter reascension
likelihood, which is itself a point the talk makes.

### Session 5 — Guidance efficiency fit both ways

**Objective:** explain what guidance efficiency is, why it needs its own model,
and what fitting it by maximum likelihood versus Bayesian buys and costs.

Fit the GE relationship both ways on simulated data with known truth: maximum
likelihood in `glmmTMB` and Bayesian by hand. Simulate-only, so the truth stays
known and the two fits can be judged against it; the truth values are drawn to
look like the real MY2025 GE-versus-spill shape (see `data/`) without the script
depending on the data to run. The exercise pushes the fit to a separation case
where maximum likelihood diverges and only the prior keeps it finite.

**Locate:** `smoltEASE` `fit_ge_model()` (JAGS multistate mark-recapture for the
route-selection probability) and `prep_ge_data()` in `R/fit_ge_model.R` and
`R/prep_ge_data.R`; `escapeLGD` has no GE model, because nighttime passage plays
the structurally analogous role and is estimated as a plain binomial rate in
`nightFall()`.

### Session 6 — What the composed interval actually claims

**Objective:** explain what the SCRAPI2 interval claims to cover and whether a
simulation says it delivers that coverage, and how the adult interval compares.

The SCRAPI2 interval is composed of a nonparametric bootstrap stacked on
posterior draws of GE and GSI. State what that composed interval claims, then
run a simulation study to check whether it has the coverage the claim implies.
Do the same for the adult interval from `escapeLGD` and compare. The Monte Carlo
identity — that drawing a posterior value and pushing it through the calculation
is integrating the nuisance out — is used here rather than given its own session.

**Locate:** `smoltEASE` `SCRAPI2()` CI construction (`quantile(theta.b, ...)`
over a bootstrap that folds in GE and GSI draws, where `ge_day_mat[, b]` enters
`thetahat()`); `escapeLGD` `apply_fallback_rates()` CI construction (`quantile()`
over composition bootstrap times fallback bootstrap) and `HNC_expand_unkGSI()`
GSI posterior draw columns.

### Session 7 — The real run, checked by hand

**Objective:** explain that the whole smolt machine is the session-1 move applied
to real inputs, by expanding one real day by hand and watching production
reproduce it and the run total.

No simulation. Read the real MY2025 steelhead files in `data/` and run production
`smoltEASE::SCRAPI2()` on them with fixed guidance efficiency and a bootstrap size
cut well below the 5000 default for runtime. Pick one day, expand it by hand as
Tally / (SampleRate × GuidanceEfficiency), and show that hand number reproduces
what SCRAPI2 reports for that day, and that the daily expansions sum to the total
SCRAPI2 prints. Same move as session 1, real inputs, checkable with a calculator.

**Locate:** `smoltEASE` `SCRAPI2()` in `R/SCRAPI2.R`, where `pass$estimated` is
`SampleCount / (SampleRate * GuidanceEfficiency)` summed by week and stratum;
`escapeLGD` `expand_wc_binom_night()`, `round(wc / wc_prop)` summed by week.

### Session 8 — Composition two ways, accounting versus likelihood

**Objective:** explain why splitting an expanded total into origin groups can be
done by bookkeeping or by a likelihood, why the two agree on the point estimate,
and why only the likelihood tells you honestly what an imperfect PBT tag rate
costs your certainty. This is the sharpest adult/smolt divergence and nothing
else in the talk demonstrates it.

Simulate a stratum with known origin composition and known PBT tag rates. Estimate
the proportions both ways — the accounting expansion and the multinomial MLE — show
the point estimates coincide because accounting is the interior MLE, then profile
the likelihood over the wild fraction and show its interval widens as the tag rate
worsens, an honesty about uncertainty the single accounting number cannot express.

**Locate:** `escapeLGD` `HNC_expand_unkGSI()` in `R/wrappers_HNC_expand.R`, which
takes `method = c("Account", "MLE")` and branches to `HNC_expand_one_strat()` or
`HNC_expand_one_strat_MLE()`; the MLE path runs `PBT_expand_calc_MLE()` in
`R/composition_estimation_utils.R` (multinomial over softmax proportions, analytic
gradient, optim BFGS) and the accounting path `PBT_expand_calc()`. `smoltEASE` has
no likelihood path; `thetahat()` in `R/SCRAPI2.R` does inverse-sample-rate
weighting and `prop.table` only.

### Session 9 — Talk assembly

**Objective:** deliver the talk from two or three figures.

Assemble the two or three figures that carry the whole thing: the shared
skeleton, the profiling core, and the coverage result. No new estimator work;
this session selects and finishes figures already produced.

## Rules this repository holds itself to

These rules govern the **talk track** (`R/session*.R`). The interactive `learnr`
course under `inst/tutorials/estimator-histories/` is cumulative by design and
is deliberately exempt from the one-idea, one-figure, and one-minute rules; see
`LEARNING_PATH.md`.

- Each session teaches one idea and one idea only: simulate from a known truth,
  build the estimator by hand, recover it against that truth, show it in one
  figure, and Locate it in the two production sources. If a session needs a second
  figure to make its point, it is two sessions.
- If a helper function is used once, inline it. No defensive input checking —
  these are teaching scripts, not package code.
- No script takes more than a minute to run.
- No new session is added without deleting one. Tangents go to `BACKLOG.md`.
- No banner comments (no `# =====`, no `# -----`). Plain comments only.
- Every heading in every markdown file has real text. No empty or placeholder
  headings.
- No dates are set anywhere in this repository. I add dates myself.

## Repository layout

- `PLAN.md` — this file.
- `README.md` — one-screen orientation.
- `BACKLOG.md` — everything cut from the earlier 22-section plan, each with one
  line saying why it was cut.
- `R/sessionNN_*.R` — one self-contained script per session.
- `docs/sessionNN_explain.md` — the three-minute plain-English explanation each
  session writes.
- `docs/session02_shared_skeleton.md` — the component-by-component comparison
  table.
- `data/` — real MY2025 steelhead inputs, both the reference shapes the
  simulations imitate and the actual inputs session 7 runs `SCRAPI2()` on; see
  `data/README.md`.
