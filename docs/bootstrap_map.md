# Bootstrap map: what each estimator resamples and where the draw enters

Read out of the `smoltEASE` and `escapeLGD` source. Every row names the real
production function that implements it. "Draw type" is one of: parametric
bootstrap (resample a fitted distribution), nonparametric bootstrap (resample
observed units), posterior draw (a fitted Bayesian posterior), or plug-in (a
fixed constant, not resampled).

## Smolts: SCRAPI2 (`smoltEASE`)

| Quantity | What is resampled / drawn | Unit | Held fixed | Where the draw enters the estimator | Draw type | Production function |
| --- | --- | --- | --- | --- | --- | --- |
| Passage counts | daily trap count | one day | strata, SR, sample rate | `est_daily_b <- round(Tally / ptrue_b)`, then `cntstar[i] <- rbinom(1, est_daily_b[i], ptrue_b[i])`; expanded as `Tally / ptrue_b` | parametric bootstrap | `SCRAPI2()` bootstrap loop, `R/SCRAPI2.R` |
| Fish (composition) | genotyped fish, **uniform within collapsed stratum** | one fish | SR (`trap_rate * genotype_rate`) | resampled fish enter `thetahat()`, weighted by `1/Fish$SR` in `Primarystrata` | nonparametric bootstrap | `SCRAPI2()`: `sample.int(nrow(jw), replace = TRUE)` (no `prob`) |
| Rearing units (wild split) | rearing records, uniform within stratum | one rearing record | — (True rebuilt per replicate, see GE) | `1/RearData$True` in `thetahat()` gives `PWild` | nonparametric bootstrap | `SCRAPI2()`: `sample.int(nrow(jw), replace = TRUE)` |
| Guidance efficiency (GE) | one coherent posterior column per replicate | whole season (all days share the column) | GE excluded from SR | `ptrue_b <- rate * ge_day_mat[, b]` in passage expansion, and `True_b <- ptrue_b[...]` in the rearing weight | posterior draw | columns from `generate_ge_draws()`; consumed in `SCRAPI2()` |
| Genetic stock (GSI) | one per-fish assignment column per replicate | whole fish sample | — | `ap_b[[FISHpndx]] <- gsiDraws[[gsi_idx_boot[b] + 1L]][im]` replaces the stock column before `thetahat()` | posterior draw | `gsiDraws` argument to `SCRAPI2()`; test generator `sim_gsi_draws()` |
| Trap/genotype sampling rate (SR) | not resampled | — | fixed across replicates | `1/SR` weight inside `thetahat()` (`Primarystrata`, `getAvgProp`, `Freqs`) | plug-in | computed in `SCRAPI2()` as `trap_rate * subrate` |

Point estimate: `pointVec <- colMeans(theta.b)` (mean of the B rows). Interval:
`quantile(theta.b[, j], c(alph/2, 1 - alph/2))` (percentile). GE and GSI columns
are drawn with independent indices (`ge_idx_boot`, `gsi_idx_boot`) so two
unrelated MCMC chains are not paired by shared column position. Composition
proportions sum to 1 within every replicate, so stock and age/sex cells add to
their parent total row by row; the mean preserves that additivity, which is why
`pointEst = "mean"` is enforced.

## Original SCRAPI vs current SCRAPI2

The one estimator-design difference that changes the target, not just the
spread:

- **Original SCRAPI** (`smoltEASE/R/SCRAPI.R`) resamples fish **weighted by SR**:
  `sample.int(nwk, replace = TRUE, prob = unlist(justwk$SR))` (and rearing units
  weighted by `True`), while `thetahat` still divides by `SR`. Resampling
  proportional to SR and then dividing by SR cancels, so the bootstrap drifts
  toward the **unweighted** sample composition, discarding the inverse-
  probability correction.
- **Current SCRAPI2** (`smoltEASE/R/SCRAPI2.R`) resamples fish and rearing units
  **uniformly** within each stratum (`sample.int(..., replace = TRUE)`, no
  `prob`) and keeps `1/SR` and `1/True` only inside `thetahat()`. This is the
  behavior the interactive tutorial's coverage exercise reproduces: uniform
  resampling covers the true composition near the nominal rate; SR-weighted
  resampling collapses coverage because it estimates the wrong quantity.

## Adults: EASE (`escapeLGD`)

| Quantity | What is resampled / drawn | Unit | Held fixed | Where the draw enters the estimator | Draw type | Production function |
| --- | --- | --- | --- | --- | --- | --- |
| Window count | expanded daytime count | one statistical week | `wc_prop` (5/6 for EASE, 1 for old SCOBI) | `rbinom(boots, round(wc / wc_prop), wc_prop) / wc_prop` | parametric bootstrap | `expand_wc_binom_night()`, `R/night_fall_reascend_wc_binom.R` |
| Counting-window fraction (`wc_prop`) | not resampled | — | fixed constant | divides the count: `wc / wc_prop` | plug-in | `expand_wc_binom_night()` argument; set in `EASE-code` Step 3 |
| Nighttime passage | binomial resample of the PIT night rate | PIT tags per stratum | totalPass | `p_night` bootstrap; expansion `wc / (1 - p_night)` | nonparametric bootstrap (binomial-equivalent) | `nightFall()` (`nightPassage_rates`), `expand_wc_binom_night()` |
| Composition | trap fish per stratum (PBT + GSI accounting or multinomial MLE) | one fish | tag rates | resampled counts feed `PBT_breakdown()` / `subGroup_breakdown()` | nonparametric bootstrap (+ MLE) | `HNC_expand()` / `ascension_composition()`; `composition_estimation_utils.R` |
| Genetic stock (GSI) | posterior draw columns, one per iteration | fish sample | — | per-iteration stock assignment in the composition breakdown | posterior draw | `HNC_expand_unkGSI()`, `R/wrappers_HNC_expand.R` |
| Fallback | binomial resample of reascension + spillway counts, re-fit each iteration | PIT tags per stratum/stockGroup | totalPass, totalFall | two-parameter likelihood re-optimized per replicate; escapement scaled by fallback | nonparametric bootstrap + MLE | `nightFall()` with `optimllh()`/`fallback_log_likelihood()` (`fallback_reascend_likelihood.R`) |
| Reascension | folded into the fallback likelihood | PIT tags | — | `P(reascend | fallback)`; when no spillway data, assumed 1 (no fallback without reascension) | nonparametric bootstrap + MLE | `nightFall()` |

Interval: percentile `quantile()` over the composition bootstrap multiplied by
the fallback bootstrap, assembled in `apply_fallback_rates()`
(`R/apply_fallback_rates_summarize.R`).

## The one structural difference between the two detection terms

Smolt GE is a fitted Bayesian **posterior** (`fit_ge_model()` +
`generate_ge_draws()`), so it can be **miscentered** by the model: a biased GE
posterior biases escapement and no wider interval fixes it. Adult nighttime
passage is a bootstrapped **binomial rate** from PIT tags, so it carries only
sampling noise (poor PIT sample sizes widen it) and cannot be miscentered by a
model. Everything else in the two skeletons lines up piece for piece.
