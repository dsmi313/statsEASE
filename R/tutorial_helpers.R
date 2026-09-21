# Helpers for the estimator-histories learnr tutorial.
# These only generate the simulated fishery and format output. The estimator
# itself, expansion, resampling, and uncertainty propagation, is written by the
# learner inside the tutorial, never hidden here.

# One simulated smolt season, preserved across the whole tutorial.
#
# What is random vs fixed: truth (total smolts, per-week passage, wild fraction
# by stratum, stock proportion, guidance efficiency by stratum) is fixed by the
# seed and treated as known. The observed trap counts, the rear-type sample, and
# the genotyped wild sample are the one realized dataset the learner conditions
# on. GE and GSI posterior draws stand in for the JAGS / genetics-lab posteriors
# that generate_ge_draws() and the GSI pipeline supply in production.
simulate_smolt_season <- function(seed = 6L, B = 1000L) {
  set.seed(seed)

  # Six statistical weeks collapsed into three strata (two weeks each), the same
  # collapse SCRAPI2 does through its Collapse column.
  week    <- 1:6
  stratum <- c(1L, 1L, 2L, 2L, 3L, 3L)

  # Truth. Total smolts per week (wild + hatchery), the stratum-level guidance
  # efficiency, and the wild fraction by stratum. Stratum 3 has deliberately low
  # GE (0.08) so the 1/GE expansion develops a right tail the learner can see.
  N_week     <- c(1200L, 2000L, 2400L, 1800L, 900L, 700L)
  GE_stratum <- c(0.22, 0.14, 0.08)
  GE_week    <- GE_stratum[stratum]
  rate_week  <- rep(0.12, 6)                    # trap sample rate
  Ptrue_week <- rate_week * GE_week             # inclusion prob of the trap count
  pWild_str  <- c(0.75, 0.60, 0.85)             # wild fraction by stratum
  pWild_week <- pWild_str[stratum]

  # Observed trap count: how many of the true passers (all rear types) were
  # counted at the trap.
  Tally <- rbinom(6, N_week, Ptrue_week)

  weeks <- data.frame(week = week, stratum = stratum, N_week = N_week,
                      rate = rate_week, GE_true = GE_week, Ptrue = Ptrue_week,
                      pWild = pWild_week, Tally = Tally)

  # Rear-type sample (wild W vs hatchery H). A smolt enters this sample at rate
  # True = rate * GE (trap rate times guidance efficiency), so GE enters the
  # wild-proportion estimate here. Each rear fish carries its week, so the
  # estimator can rebuild True_b = rate * GE_b from that replicate's GE column.
  # The wild proportion is estimated by 1/True weighting within stratum.
  rear <- do.call(rbind, lapply(1:6, function(i) {
    n_samp <- max(rbinom(1, N_week[i], Ptrue_week[i]), 1L)
    r      <- ifelse(runif(n_samp) < pWild_week[i], "W", "H")
    data.frame(week = i, stratum = stratum[i], rear = r,
               stringsAsFactors = FALSE)
  }))

  # Genotyped wild sample for stock composition. Among wild smolts, stock A vs B
  # with true proportion pA. Teaching device: a fish's genetic sampling rate SR
  # depends on its stock (stock A is genotyped at a lower rate). In production SR
  # is a per-day trap_rate * genotyping_rate that EXCLUDES GE; the real analogue
  # is stock mix correlating with the day-to-day sampling rate. Either way, 1/SR
  # is a Horvitz-Thompson correction that is only honest if the resample does not
  # also weight by SR (which would cancel it). SR carries no GE, so GE never
  # directly reweights composition.
  pA  <- 0.35
  srA <- 0.10
  srB <- 0.30
  fish <- do.call(rbind, lapply(1:6, function(i) {
    n_wild <- round(N_week[i] * pWild_week[i])
    stock  <- ifelse(runif(n_wild) < pA, "A", "B")
    sr     <- ifelse(stock == "A", srA, srB)    # SR = rate * genotype_rate, no GE
    keep   <- runif(n_wild) < sr
    if (!any(keep)) keep[1] <- TRUE
    data.frame(week = i, stratum = stratum[i], stock = stock[keep],
               SR = sr[keep],           # trap_rate * genotype_rate, GE excluded
               stringsAsFactors = FALSE)
  }))
  fish$fishID <- seq_len(nrow(fish))
  n_fish <- nrow(fish)

  # GE posterior draws: weeks x B. Each COLUMN is one coherent season-wide draw
  # (all weeks share the stratum draw), exactly how SCRAPI2 consumes one
  # ge_day_mat column per bootstrap replicate.
  ge_draws <- matrix(0, nrow = 6, ncol = B)
  for (b in seq_len(B)) {
    g_s <- rbeta(3, GE_stratum * 250, (1 - GE_stratum) * 250)
    ge_draws[, b] <- g_s[stratum]
  }

  # GSI posterior draws: n_fish x B. Each column is one coherent per-fish stock
  # assignment draw; a fish keeps its observed stock with probability 0.9.
  gsi_draws <- matrix(fish$stock, nrow = n_fish, ncol = B)
  flip  <- matrix(runif(n_fish * B) > 0.9, nrow = n_fish, ncol = B)
  other <- ifelse(fish$stock == "A", "B", "A")
  gsi_draws[flip] <- other[row(gsi_draws)[flip]]

  WildSmolts <- sum(N_week * pWild_week)

  list(weeks = weeks, rear = rear, fish = fish,
       ge_draws = ge_draws, gsi_draws = gsi_draws,
       truth = list(N_total = sum(N_week), WildSmolts = round(WildSmolts),
                    pA = pA, stockA = round(WildSmolts * pA),
                    stockB = round(WildSmolts * (1 - pA))),
       params = list(srA = srA, srB = srB, pA = pA, B = B))
}

# One simulated adult season for the SCOBI-to-EASE module. Window counts by
# statistical week, PIT-based nighttime-passage and fallback/reascension data,
# and a trap composition sample with GSI draws. Two truths are stored because
# fallback separates them: ascensions (fish that climbed the ladder) and
# escapement (ascensions that did not fall back without reascending).
simulate_adult_season <- function(seed = 6L, B = 1000L) {
  set.seed(seed)
  sWeek   <- 1:6
  stratum <- c(1L, 1L, 2L, 2L, 3L, 3L)

  N_week  <- c(1200L, 2600L, 3100L, 2400L, 1500L, 800L)  # true ascensions/week
  wc_prop <- 5 / 6                              # counting window open fraction
  p_night <- 0.12                               # fraction passing at night
  p_fall  <- 0.05                               # net fallback without reascension
  pA      <- 0.40                               # stock A share of escapement

  # Observed daytime window count: day fish passing while the window is open.
  seen_day <- rbinom(6, N_week, wc_prop * (1 - p_night))
  wc <- data.frame(sWeek = sWeek, stratum = stratum, wc = seen_day)

  night_tot  <- c(300L, 320L, 290L)
  night_pass <- rbinom(3, night_tot, p_night)
  night <- data.frame(stratum = 1:3, nightPass = night_pass, totalPass = night_tot)

  fall_tot   <- c(280L, 300L, 275L)
  fall_reasc <- rbinom(3, fall_tot, p_fall)
  fallback <- data.frame(stratum = 1:3, numReascend = fall_reasc, totalPass = fall_tot)

  # Trap composition sample: two stock groups, proportion pA true, by stratum.
  n_comp <- 900L
  comp <- data.frame(
    stratum = sample(1:3, n_comp, replace = TRUE),
    stock   = ifelse(runif(n_comp) < pA, "A", "B"),
    stringsAsFactors = FALSE)
  comp$fishID <- seq_len(n_comp)

  # GSI posterior draws for the composition sample: n_comp x B coherent columns.
  gsi_draws <- matrix(comp$stock, nrow = n_comp, ncol = B)
  flip  <- matrix(runif(n_comp * B) > 0.9, nrow = n_comp, ncol = B)
  other <- ifelse(comp$stock == "A", "B", "A")
  gsi_draws[flip] <- other[row(gsi_draws)[flip]]

  true_ascensions <- sum(N_week)
  true_escapement <- sum(N_week * (1 - p_fall))

  list(wc = wc, night = night, fallback = fallback, comp = comp,
       gsi_draws = gsi_draws, wc_prop = wc_prop,
       truth = list(true_ascensions = true_ascensions,
                    true_escapement = round(true_escapement),
                    pA = pA,
                    escA = round(true_escapement * pA),
                    escB = round(true_escapement * (1 - pA)),
                    p_night = p_night, p_fall = p_fall))
}

# Sample skewness, base R. Positive means a right tail.
skewness <- function(x) {
  x <- x[is.finite(x)]
  mean((x - mean(x))^3) / stats::sd(x)^3
}

# One-line summary of a bootstrap column against a known truth.
boot_summary <- function(theta, truth, alpha = 0.10) {
  ci <- stats::quantile(theta, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  data.frame(
    truth    = truth,
    mean     = mean(theta),
    median   = stats::median(theta),
    bias     = mean(theta) - truth,
    rel_bias = (mean(theta) - truth) / truth,
    lci      = ci[1],
    uci      = ci[2],
    width    = ci[2] - ci[1],
    skew     = skewness(theta))
}

# One row of a cumulative "watch the estimator evolve" table.
stage_row <- function(stage, sources, theta, truth, alpha = 0.10,
                      coverage = NA_real_) {
  ci <- stats::quantile(theta, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  data.frame(
    stage    = stage,
    sources  = sources,
    estimate = round(mean(theta)),
    truth    = round(truth),
    bias     = round(mean(theta) - truth),
    rel_bias = round((mean(theta) - truth) / truth, 3),
    lci      = round(ci[1]),
    uci      = round(ci[2]),
    width    = round(ci[2] - ci[1]),
    skew     = round(skewness(theta), 2),
    coverage = coverage,
    stringsAsFactors = FALSE)
}

# Histogram of a bootstrap distribution with truth and the point estimate drawn on.
plot_boot_hist <- function(theta, truth = NULL, point = NULL, main = "",
                           xlab = "estimate") {
  hist(theta, breaks = 40, col = "grey85", border = "white",
       main = main, xlab = xlab)
  if (!is.null(truth)) abline(v = truth, col = "firebrick", lwd = 3)
  if (!is.null(point)) abline(v = point, col = "steelblue", lwd = 3, lty = 2)
}
