# Session 2: the shared skeleton
# Objective: explain, component by component, how the adult and smolt estimators
# are the same machine, and name the two or three places they genuinely are not.
# No simulation here; the deliverable is one table, read out of the two sources,
# written to docs/session02_shared_skeleton.md. Columns: generic component;
# escapeLGD (adults); smoltEASE (smolts); same or different, and how.

tbl <- c(
"| Component | escapeLGD (adults) | smoltEASE (smolts) | Same or different |",
"| --- | --- | --- | --- |",
"| Total count | window count wc (expand_wc_binom_night) | daily trap count Tally (thetahat in SCRAPI2) | same: a raw count at the dam |",
"| Sampled fraction | wc_prop, counting-window open fraction | SampleRate, trap sample rate | parameterization: adults one term, smolts one of two |",
"| Detection expansion | nighttime passage, 1/(1-p_night) (nightFall) | guidance efficiency, Ptrue = SampleRate * GE | STRUCTURAL: unseen night passage vs route selection into the bypass |",
"| Count uncertainty | binomial bootstrap rbinom(boots, round(wc/wc_prop), wc_prop)/wc_prop | binomial bootstrap rbinom(1, est_daily, Ptrue) | same mechanism: parametric binomial bootstrap on the expanded count |",
"| Composition estimator | PBT multinomial MLE via softmax optim, or accounting expansion | accounting only: inverse-SR weighting, prop.table (thetahat) | STRUCTURAL: adults have a likelihood path, smolts do not |",
"| Composition uncertainty | nonparametric resample of trap fish per stratum (sample_n) | uniform resample within stratum; 1/SR kept only inside the estimator (thetahat) | same: both nonparametric within stratum. Original SCRAPI weighted the resample by SR; SCRAPI2 corrected that to uniform |",
"| Genetic stock (GSI) uncertainty | posterior draw columns, one per iter (HNC_expand_unkGSI) | gsiDraws, one column per iter; point = colMeans(theta.b) over the B rows (SCRAPI2; n_point deprecated) | SAME: smoltEASE deliberately copies the adult GSI format |",
"| Detection-rate uncertainty | night passage: binomial bootstrap only, NO posterior | GE: Bayesian posterior draws (fit_ge_model, generate_ge_draws) | FLAG: smolts propagate a Bayesian posterior here, adults do not |",
"| Stage adjustment | fallback and reascension: two-parameter binomial likelihood MLE + bootstrap | none (smolts do not reascend); wild-fraction split instead | FLAG: a source adults propagate and smolts have no analogue for |",
"| Origin / wild split | H / HNC / W composition (HNC_expand_one_strat) | explicit PWild by rear-strata inverse weighting (thetahat) | parameterization of the same split |",
"| Interval | bootstrap percentile, quantile() (apply_fallback_rates) | bootstrap percentile, quantile(theta.b) (SCRAPI2) | same: percentile of the propagated bootstrap |")

notes <- c(
"## Uncertainty sources and how each is propagated",
"",
"- Count: both propagate it by parametric binomial bootstrap.",
"- Composition and genetic stock: both by nonparametric bootstrap plus GSI",
"  posterior draw columns; both resample uniformly within stratum. SCRAPI2 keeps",
"  the 1/SR inverse-rate weight only inside the estimator; the original SCRAPI",
"  weighted the resample by SR as well, which cancelled that correction.",
"- Detection rate: smolts push GE posterior draws through the bootstrap; adults",
"  bootstrap nighttime passage as a plain binomial rate, no posterior. This is",
"  the one source one tool models and the other does not.",
"- Fallback and reascension: adults only, a two-parameter likelihood with its",
"  own bootstrap. Smolts have no reascension, so no analogue exists.",
"- Every interval on both sides is a bootstrap percentile; neither production",
"  tool ships a delta-method interval.")

md <- c("# The shared skeleton: one estimator, two life stages", "",
        "Read out of the escapeLGD and smoltEASE source. Adults (EASE) and smolts",
        "(SCRAPI2) run the same six-part machine; the FLAG rows are where they",
        "genuinely diverge.", "", tbl, "", notes)
writeLines(md, "docs/session02_shared_skeleton.md")
cat(paste(tbl, collapse = "\n"), "\n")

# Count, from the table itself, how many components line up and how many diverge.
comp_rows <- tbl[-(1:2)]
n_comp <- length(comp_rows)
n_flag <- sum(grepl("STRUCTURAL|FLAG", comp_rows))
cat(sprintf("components compared: %d | rows carrying a divergence marker: %d\n",
            n_comp, n_flag))

# Locate. This whole session is Locate: every cell names the real file and
# function. Adults: night_fall_reascend_wc_binom.R, fallback_reascend_likelihood.R,
# composition_estimation_utils.R, HNC_expand_one_strat.R, wrappers_HNC_expand.R,
# apply_fallback_rates_summarize.R. Smolts: SCRAPI2.R, fit_ge_model.R,
# generate_ge_draws.R, sim_gsi_draws.R.

writeLines(c(
"How I would explain session 2 in three minutes",
"",
"Adults and smolts look like separate problems with separate software, but they",
"are the same assembly line. Count the fish you can, divide by the fraction you",
"saw, expand for the ones the count structurally misses, split the total into",
sprintf("groups, and shake the whole thing to get a range. Lining the two tools up gives"),
sprintf("%d components; %d of those match cleanly and only %d carry a divergence marker.",
        n_comp, n_comp - n_flag, n_flag),
"",
"Those flagged rows collapse to three talk-level differences. Smolts give guidance",
"efficiency a full statistical model and carry its uncertainty; adults treat the",
"matching quantity, night passage, as a simple rate. Adults carry fallback and",
"reascension, which smolts do not have. And adults can estimate composition by",
"likelihood, while smolts only account for it. Everything else is the same",
"machine wearing two life stages."
), "docs/session02_explain.md")
