# Helpers for the estimator-histories learnr tutorial.
# These only generate the simulated fishery and format output. The estimator
# itself, expansion, resampling, and uncertainty propagation, is written by the
# learner inside the tutorial, never hidden here.

# One simulated smolt season, preserved across the whole tutorial.
#
# What is random vs fixed: truth (total wild smolts, per-week passage, stock
# proportion, guidance efficiency by stratum) is fixed by the seed and treated
# as known. The observed trap counts and the genotyped sample are the one
# realized dataset the learner conditions on. GE and GSI posterior draws stand
# in for the JAGS / genetics-lab posteriors that generate_ge_draws() and the
# GSI pipeline supply in production.
simulate_smolt_season <- function(seed = 6L, B = 1000L) {
  set.seed(seed)

  # Six statistical weeks collapsed into three strata (two weeks each), the
  # same collapse SCRAPI2 does through its Collapse column.
  week    <- 1:6
  stratum <- c(1L, 1L, 2L, 2L, 3L, 3L)

  # Truth. True wild passage per week and the stratum-level guidance efficiency.
  # Stratum 3 has deliberately low GE (0.08) so the 1/GE expansion develops a
  # right tail the learner can see.
  N_week      <- c(900L, 1500L, 1800L, 1400L, 700L, 500L)
  GE_stratum  <- c(0.22, 0.14, 0.08)
  GE_week     <- GE_stratum[stratum]
  rate_week   <- rep(0.12, 6)                 # trap sample rate
  Ptrue_week  <- rate_week * GE_week          # inclusion prob of the trap count
  N_total     <- sum(N_week)

  # Observed trap count: how many of the true passers were counted at the trap.
  Tally <- rbinom(6, N_week, Ptrue_week)

  weeks <- data.frame(week = week, stratum = stratum, N_week = N_week,
                      rate = rate_week, GE_true = GE_week,
                      Ptrue = Ptrue_week, Tally = Tally)

  # Composition. Among wild smolts, stock A vs B with true proportion pA.
  # Teaching device: a fish's genetic sampling rate SR depends on its stock
  # (stock A is genotyped at a lower rate than stock B). In production SR is a
  # per-day trap-rate * genotyping-rate and does not depend on stock; the real
  # analogue is stock mix correlating with the day-to-day sampling rate. The
  # statistical point is identical either way: 1/SR inside the estimator is a
  # Horvitz-Thompson correction, and it is only honest if the resample does not
  # also weight by SR (which would cancel it).
  pA  <- 0.35
  srA <- 0.10                                 # stock A genotyped less often
  srB <- 0.30
  n_wild        <- 4000L                      # wild smolts eligible for genetics
  stock_true    <- ifelse(runif(n_wild) < pA, "A", "B")
  sr_fish       <- ifelse(stock_true == "A", srA, srB)
  sampled       <- runif(n_wild) < sr_fish    # genotyped inclusion
  fish <- data.frame(
    fishID  = seq_len(sum(sampled)),
    stratum = sample(1:3, sum(sampled), replace = TRUE),
    stock   = stock_true[sampled],
    SR      = sr_fish[sampled],
    stringsAsFactors = FALSE)
  n_fish <- nrow(fish)

  # GE posterior draws: weeks x B. Each COLUMN is one coherent season-wide draw
  # (all weeks share the stratum draw), exactly how SCRAPI2 consumes one
  # ge_day_mat column per bootstrap replicate. Concentration 250 gives a
  # realistic posterior spread.
  ge_draws <- matrix(0, nrow = 6, ncol = B)
  for (b in seq_len(B)) {
    g_s <- rbeta(3, GE_stratum * 250, (1 - GE_stratum) * 250)
    ge_draws[, b] <- g_s[stratum]
  }

  # GSI posterior draws: n_fish x B. Each column is one coherent per-fish stock
  # assignment draw. A fish keeps its observed stock with probability 0.9 and
  # flips otherwise, standing in for genotyping-assignment uncertainty.
  gsi_draws <- matrix(fish$stock, nrow = n_fish, ncol = B)
  flip <- matrix(runif(n_fish * B) > 0.9, nrow = n_fish, ncol = B)
  other <- ifelse(fish$stock == "A", "B", "A")
  gsi_draws[flip] <- other[row(gsi_draws)[flip]]

  list(weeks = weeks, fish = fish,
       ge_draws = ge_draws, gsi_draws = gsi_draws,
       truth = list(N_total = N_total, pA = pA,
                    stockA = round(N_total * pA),
                    stockB = round(N_total * (1 - pA))),
       params = list(srA = srA, srB = srB, pA = pA, B = B))
}

# One simulated adult season for the SCOBI-to-EASE module. Window counts by
# statistical week, PIT-based nighttime-passage and fallback/reascension data,
# and a trap composition sample. Truth is the true escapement the pieces should
# reconstruct.
simulate_adult_season <- function(seed = 6L) {
  set.seed(seed)
  sWeek   <- 1:6
  stratum <- c(1L, 1L, 2L, 2L, 3L, 3L)

  # Truth: true adults passing per week, the daytime counting-window fraction,
  # the true nighttime-passage rate, and the true net fallback rate.
  N_week    <- c(1200L, 2600L, 3100L, 2400L, 1500L, 800L)
  wc_prop   <- 5 / 6                          # counting window open fraction
  p_night   <- 0.12                           # fraction passing at night
  p_fall    <- 0.05                           # net fallback without reascension

  # Observed daytime window count: day fish that pass during the open window.
  seen_day  <- rbinom(6, N_week, wc_prop * (1 - p_night))
  wc        <- data.frame(sWeek = sWeek, stratum = stratum, wc = seen_day)

  # Nighttime-passage PIT data by stratum: totalPass tags, nightPass at night.
  night_tot   <- c(300L, 320L, 290L)
  night_pass  <- rbinom(3, night_tot, p_night)
  night <- data.frame(stratum = 1:3, nightPass = night_pass, totalPass = night_tot)

  # Fallback / reascension PIT data by stratum: of totalPass ascending tags,
  # numReascend were reascensions (net fallback signal).
  fall_tot    <- c(280L, 300L, 275L)
  fall_reasc  <- rbinom(3, fall_tot, p_fall)
  fallback <- data.frame(stratum = 1:3, numReascend = fall_reasc, totalPass = fall_tot)

  # Trap composition sample: two stock groups, proportion pA true.
  pA        <- 0.40
  n_comp    <- 900L
  comp <- data.frame(
    stratum = sample(1:3, n_comp, replace = TRUE),
    stock   = ifelse(runif(n_comp) < pA, "A", "B"),
    stringsAsFactors = FALSE)

  list(wc = wc, night = night, fallback = fallback, comp = comp,
       wc_prop = wc_prop,
       truth = list(N_total = sum(N_week), pA = pA,
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
    truth     = truth,
    mean      = mean(theta),
    median    = stats::median(theta),
    bias      = mean(theta) - truth,
    rel_bias  = (mean(theta) - truth) / truth,
    lci       = ci[1],
    uci       = ci[2],
    width     = ci[2] - ci[1],
    skew      = skewness(theta))
}

# Histogram of a bootstrap distribution with truth and the point estimate drawn on.
plot_boot_hist <- function(theta, truth = NULL, point = NULL, main = "",
                           xlab = "estimate") {
  hist(theta, breaks = 40, col = "grey85", border = "white",
       main = main, xlab = xlab)
  if (!is.null(truth))
    abline(v = truth, col = "firebrick", lwd = 3)
  if (!is.null(point))
    abline(v = point, col = "steelblue", lwd = 3, lty = 2)
}
