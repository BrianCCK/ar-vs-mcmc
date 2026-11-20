###############################################################################
# MATH 3428 Statistical Computing -- Project code
# "A Comparative Study of Acceptance-Rejection Sampling and MCMC Algorithms"
# Author: CHAN Chi Kong Brian
#
# This file contains all the code for the report:
#   Part 0  Shared settings and functions: effective sample size, Gelman-Rubin
#           statistic, burn-in rule, table formatting, plot style
#   Part 1  Experiment 1: one-dimensional targets of different shapes  (Section 3.1)
#   Part 2  Experiment 2: dimension and correlation, with MALA          (Section 3.2, App. D)
#   Part 3  Real-data example: the Challenger O-ring data               (Appendix E)
#
# Each part writes its results, tables and figures to its own folder. 
# Every setting has its own fixed random seed, so a part gives the same results whether or not the other parts are run.
#
# Requirements: base R only (the package coda is optional and used only as a check).
###############################################################################

RUN_EXP1 <- TRUE        # Part 1
RUN_EXP2 <- TRUE        # Part 2
RUN_REAL <- TRUE        # Part 3

###############################################################################
# PART 0 -- Shared settings and functions
###############################################################################

## ========================== 0.1 Common settings ==============================
B        <- 100000      # budget: evaluations of the unnormalised density f~ per run
N_REP    <- 50          # independent replicates per setting
N_PILOT  <- B/2         # length of each pilot chain used to choose the burn-in
RHAT_TOL <- 1.1         # Gelman-Rubin threshold
SEED     <- 3428        # base seed; each setting derives its own seed from it

## ========================== 0.2 Effective sample size ========================
# Geyer's (1992) initial positive sequence estimator of tau; ESS = n / tau.
ess_geyer <- function(z) {
  n <- length(z)
  if (n < 4 || var(z) == 0) return(NA_real_)
  zc  <- z - mean(z)
  f   <- fft(c(zc, numeric(n)))                               
  acv <- Re(fft(Mod(f)^2, inverse = TRUE))[1:(n - 1)]
  rho <- acv / acv[1]                                         
  m     <- floor(length(rho) / 2)
  Gamma <- rho[seq(1, 2 * m, by = 2)] + rho[seq(2, 2 * m, by = 2)]   
  neg   <- which(Gamma <= 0)
  M <- if (length(neg) > 0) neg[1] - 1 else length(Gamma)     
  n / (-1 + 2 * sum(Gamma[seq_len(M)]))
}

## ========================== 0.3 Gelman-Rubin statistic and burn-in ===========
# Course formula, R-hat = V-hat / W ,
# applied to the draws in the second half of the first t iterations of each chain.
rhat_second_half <- function(X, t) {          # X: chains in rows, draws in columns
  Z <- X[, (t / 2 + 1):t, drop = FALSE]; n <- ncol(Z)
  W <- mean(apply(Z, 1, var))                 # within-chain variance
  if (W == 0) return(Inf)
  Bv <- n * var(rowMeans(Z))                  # between-chain variance
  (W * (n - 1) / n + Bv / n) / W
}

# Burn-in rule used in every part. X holds the pilot chains (from dispersed starting
# points) in its rows. R-hat is computed for t = 200, 400, ...; if t* is the first t
# after which R-hat stays below RHAT_TOL, the burn-in is b = t*/2, the discarded half.
# If R-hat never settles, b = max_b and the setting is reported as not converged.
burnin_from_pilots <- function(X, max_b) {
  grid  <- seq(200, ncol(X), by = 200)
  rhat  <- sapply(grid, function(t) rhat_second_half(X, t))
  above <- which(rhat >= RHAT_TOL)
  b <- if (length(above) == 0) grid[1] / 2
       else if (max(above) == length(grid)) max_b             
       else min(max_b, grid[max(above) + 1] / 2)
  c(b = b, rhat_end = rhat[length(grid)])
}

## ========================== 0.4 Table formatting =============================
# 3 significant figures, thousands separators, powers of ten for very large or small
fmt <- function(x, sig = 3, small = 1e-3) sapply(x, function(v) {
  if (is.na(v) || !is.finite(v)) return("--")
  if (v != 0 && (abs(v) >= 1e5 || abs(v) < small)) {
    e <- floor(log10(abs(v))); m <- round(v / 10^e, 1)
    if (abs(m) >= 10) { m <- m / 10; e <- e + 1 }               
    return(sprintf("$%s\\times10^{%d}$", formatC(m, digits = 1, format = "f"), e))
  }
  s <- sub("\\.$", "", formatC(signif(v, sig), digits = sig, format = "fg", flag = "#", big.mark = ","))
  if (v < 0) paste0("$", s, "$") else s                         
})
fx <- function(v, k = 3) { r <- round(v, k); if (r == 0) r <- 0; sprintf("$%s$", formatC(r, digits = k, format = "f")) }
fmt_mse <- function(v) {
  e <- floor(log10(v)); m <- round(v / 10^e, 1)
  if (m >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("$%.1f\\mathrm{e}{%d}$", m, e)
}
DG <- "$^{\\dagger}$"; DD <- "$^{\\ddagger}$"; SEC <- "$^{\\S}$"
mark <- function(flag, notconv)
  if (flag && notconv) "$^{\\dagger\\ddagger}$" else if (flag) DG else if (notconv) DD else ""
NOTE_BOTH <- "$\\dagger$: the bias exceeds two standard errors; $\\ddagger$: the pilot chains did not converge, so $b = B/2$."
NOTE_CONV <- "$\\ddagger$: the pilot chains did not converge, so $b = B/2$."
foot <- function(text, size = "\\footnotesize") c("\\par\\vspace{3pt}", sprintf("{%s %s\\par}", size, text))
write_tex <- function(lines, dir, file) writeLines(lines, file.path(dir, file))

## ========================== 0.5 Plot style ===================================

STYLE <- data.frame(sampler = c("AR", "IMH", "RWM", "Gibbs", "MALA"),
                    col = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"),
                    pch = c(16, 17, 15, 18, 8), lty = c(1, 2, 4, 5, 3), stringsAsFactors = FALSE)
INK <- "#333333"; GRID <- "#e6e6e6"


###############################################################################
# PART 1 -- Experiment 1: One-Dimensional Targets of Different Shapes (Section 3.1)
#
# Targets  : Gamma(2.5, 1)  (skewed), Student t_3  (heavy tails),
#            0.5 N(-Delta, 1) + 0.5 N(Delta, 1)  (two modes), Delta = 1, 2, 4, 6
# Samplers : AR, IMH (same proposal g as AR), RWM tuned to an acceptance rate of 0.44
# Estimands: the mean E_f[X] and the tail probability P_f(X > q_0.95) = 0.05
# Output   : CSV files, LaTeX tables and one PDF figure in Exp1_output/
###############################################################################
run_experiment1 <- function() {

  ## ---------------------------- 1.0 Settings ---------------------------------
  ACC_TARGET  <- 0.44                # acceptance rate RWM is tuned to (lectures, d = 1)
  DELTAS      <- c(1, 2, 4, 6)       # half-distance between the two modes
  LAMBDAS     <- c(0.2, 0.4, 0.8)    # rates of the exponential proposal for the Gamma target
  SEED_OFFSET <- 1e7                 
  OUT_DIR     <- "Exp1_output"
  dir.create(OUT_DIR, showWarnings = FALSE)

  ## ---------------------------- 1.1 Targets ----------------------------------
  make_gamma <- function(a) list(
    key = "gamma", label = sprintf("Gamma(%g, 1)", a), tex = sprintf("Gamma$(%g, 1)$", a),
    lf = function(x) { out <- rep(-Inf, length(x)); p <- x > 0
                       out[p] <- (a - 1) * log(x[p]) - x[p]; out },
    logZ = lgamma(a), rexact = function(n) rgamma(n, a),
    cdf = function(x) pgamma(x, a), q = function(p) qgamma(p, a), mean = a, var = a)

  make_t <- function(nu) list(
    key = "t", label = sprintf("t_%g", nu), tex = sprintf("$t_{%g}$", nu),
    lf = function(x) -(nu + 1) / 2 * log1p(x^2 / nu),
    logZ = 0.5 * log(nu) + lbeta(0.5, nu / 2), rexact = function(n) rt(n, nu),
    cdf = function(x) pt(x, nu), q = function(p) qt(p, nu), mean = 0, var = nu / (nu - 2))

  make_mixture <- function(Delta) {
    cdf <- function(x) 0.5 * pnorm(x - Delta) + 0.5 * pnorm(x + Delta)
    list(key = paste0("mix", Delta), label = sprintf("Mixture, Delta = %g", Delta),
         tex = sprintf("Mixture, $\\Delta = %g$", Delta),
         lf = function(x) -(abs(x) - Delta)^2 / 2 + log1p(exp(-2 * Delta * abs(x))),
         logZ = log(2) + 0.5 * log(2 * pi),
         rexact = function(n) rnorm(n, sample(c(-Delta, Delta), n, replace = TRUE)),
         cdf = cdf,
         q = function(p) sapply(p, function(pp)
           uniroot(function(x) cdf(x) - pp, c(-Delta - 10, Delta + 10), tol = 1e-12)$root),
         mean = 0, var = 1 + Delta^2, Delta = Delta)
  }

  ## ---------------------------- 1.2 Proposals for AR and IMH ------------------
  log_sup_ratio <- function(lf, lg, lower, upper) {
    x  <- seq(lower, upper, length.out = 2e5 + 1)
    i  <- which.max(lf(x) - lg(x)); h <- x[2] - x[1]
    optimize(function(z) lf(z) - lg(z), c(x[i] - h, x[i] + h), maximum = TRUE, tol = 1e-12)$objective
  }
  prop_exp <- function(lam, a) {                 # g = Exp(lam) for Gamma(a, 1); needs lam < 1
    xs <- (a - 1) / (1 - lam)                    # maximiser of f~/g
    list(label = sprintf("Exp(%g)", lam), tex = sprintf("Exp$(%g)$", lam),
         r = function(n) rexp(n, lam), lg = function(x) log(lam) - lam * x,
         logM = (a - 1) * log(xs) - (1 - lam) * xs - log(lam), range = c(1e-8, 200))
  }
  prop_cauchy <- function(nu) {                  # g = Cauchy(0, 1) for t_nu; f~/g is largest at x = 1
    list(label = "Cauchy(0, 1)", tex = "Cauchy$(0, 1)$",
         r = function(n) rcauchy(n), lg = function(x) dcauchy(x, log = TRUE),
         logM = log(pi) + log(2) - (nu + 1) / 2 * log1p(1 / nu), range = c(-50, 50))
  }
  prop_normal <- function(tg) {                  # g = N(0, sd^2) with the variance of the target
    s  <- sqrt(tg$var)
    lg <- function(x) dnorm(x, 0, s, log = TRUE)
    list(label = sprintf("N(0, %g)", tg$var), tex = sprintf("$N(0, %g)$", tg$var),
         r = function(n) rnorm(n, 0, s), lg = lg,
         logM = log_sup_ratio(tg$lf, lg, 0, tg$Delta + 10) + 1e-10, range = c(0, tg$Delta + 10))
  }

  ## ---------------------------- 1.3 Samplers ---------------------------------
  # Acceptance-rejection (vectorised): exactly B proposals = B evaluations
  ar_sampler <- function(tg, pr, B) {
    y    <- pr$r(B)
    lw   <- tg$lf(y) - pr$lg(y) - pr$logM          # log{ f~(y) / (M g(y)) }, must be <= 0
    keep <- log(runif(B)) <= lw
    list(x = y[keep], acc = mean(keep), max_ratio = exp(max(lw)))
  }
  # Independence MH with the same proposal as AR
  imh_sampler <- function(tg, pr, n_iter, x0) {
    y <- pr$r(n_iter); lwy <- tg$lf(y) - pr$lg(y); logu <- log(runif(n_iter))
    x <- numeric(n_iter); cur <- x0; cur_lw <- tg$lf(x0) - pr$lg(x0); n_acc <- 0
    for (t in seq_len(n_iter)) {
      if (logu[t] <= lwy[t] - cur_lw) { cur <- y[t]; cur_lw <- lwy[t]; n_acc <- n_acc + 1 }
      x[t] <- cur
    }
    list(x = x, acc = n_acc / n_iter)
  }
  # Random-walk Metropolis with proposal N(x, sigma^2)
  rwm_sampler <- function(tg, n_iter, x0, sigma) {
    z <- sigma * rnorm(n_iter); logu <- log(runif(n_iter))
    x <- numeric(n_iter); cur <- x0; cur_lf <- tg$lf(cur); n_acc <- 0
    for (t in seq_len(n_iter)) {
      y <- cur + z[t]; ly <- tg$lf(y)
      if (logu[t] <= ly - cur_lf) { cur <- y; cur_lf <- ly; n_acc <- n_acc + 1 }
      x[t] <- cur
    }
    list(x = x, acc = n_acc / n_iter)
  }
  run_mcmc <- function(sampler, tg, pr, n_iter, x0, step = NA) {
    if (sampler == "IMH") imh_sampler(tg, pr, n_iter, x0) else rwm_sampler(tg, n_iter, x0, step)
  }

  ## ---------------------------- 1.4 Setup: step size and burn-in --------------
  # RWM: sigma is found by bisection (on the log scale) so that the acceptance rate is
  # ACC_TARGET; the acceptance rate decreases as sigma grows.
  choose_step <- function(sampler, tg, n_iter = 5000, n_bisect = 16) {
    if (sampler != "RWM") return(NA_real_)
    lo <- sqrt(tg$var) / 100; hi <- sqrt(tg$var) * 100
    for (i in seq_len(n_bisect)) {
      mid <- sqrt(lo * hi)
      acc <- rwm_sampler(tg, n_iter, tg$rexact(1), mid)$acc
      if (acc > ACC_TARGET) lo <- mid else hi <- mid
    }
    sqrt(lo * hi)
  }
  # Pilot chains start at the 0.1%, 10%, 90% and 99.9% quantiles of the target
  choose_burnin <- function(sampler, tg, pr, step) {
    if (sampler == "AR") return(c(b = 0, rhat_end = NA))
    starts <- tg$q(c(0.001, 0.1, 0.9, 0.999))
    X <- t(sapply(starts, function(x0) run_mcmc(sampler, tg, pr, N_PILOT, x0, step)$x))
    burnin_from_pilots(X, max_b = B / 2)
  }

  ## ---------------------------- 1.5 One setting ------------------------------
  # Runs one sampler on one target: setup (timed), then N_REP replicates at budget B.
  # Chains start at the 99.9% quantile of the target, a deliberately poor point.
  run_setting <- function(tg, pr, sampler, setting_id) {
    x0 <- tg$q(0.999); q95 <- tg$q(0.95)
    set.seed(SEED + SEED_OFFSET + 1e5 * setting_id)
    t_setup <- system.time({
      step <- choose_step(sampler, tg)
      bi   <- choose_burnin(sampler, tg, pr, step)
    })[["elapsed"]]
    b <- bi[["b"]]
    reps <- do.call(rbind, lapply(seq_len(N_REP), function(r) {
      set.seed(SEED + SEED_OFFSET + 1e5 * setting_id + r)
      secs <- system.time(
        out <- if (sampler == "AR") ar_sampler(tg, pr, B) else run_mcmc(sampler, tg, pr, B, x0, step)
      )[["elapsed"]]
      z <- if (sampler == "AR") out$x else out$x[(b + 1):B]          # discard burn-in
      tail_ind <- as.numeric(z > q95)
      data.frame(rep = r, n_kept = length(z), acc = out$acc, secs = secs,
                 est_mean = mean(z), est_tail = mean(tail_ind),
                 ess_mean = if (sampler == "AR") length(z) else ess_geyer(z),
                 ess_tail = if (sampler == "AR") length(z) else ess_geyer(tail_ind),
                 ess_coda = if (r == 1 && sampler != "AR" && requireNamespace("coda", quietly = TRUE))
                   as.numeric(coda::effectiveSize(z)) else NA,
                 max_ratio = if (sampler == "AR") out$max_ratio else NA)
    }))
    cbind(target = tg$key, proposal = if (is.null(pr)) NA else pr$label, sampler = sampler,
          c_theory = if (is.null(pr)) NA else exp(pr$logM - tg$logZ),
          step = step, burnin = b, rhat_pilot_end = bi[["rhat_end"]], setup_secs = t_setup, reps)
  }

  ## ---------------------------- 1.6 Summaries --------------------------------
  summarise_setting <- function(df, tg) {
    truth <- list(mean = c(theta = tg$mean, var_h = tg$var),
                  tail = c(theta = 0.05,    var_h = 0.05 * 0.95))
    do.call(rbind, lapply(names(truth), function(h) {
      est <- df[[paste0("est_", h)]]; ess <- df[[paste0("ess_", h)]]
      th  <- truth[[h]][["theta"]]
      ess[is.na(ess)] <- 0                                    # chain never moved: ESS = 0
      mse <- mean((est - th)^2)
      data.frame(
        target = df$target[1], proposal = df$proposal[1], sampler = df$sampler[1], estimand = h,
        bias = mean(est) - th, se_bias = sd(est) / sqrt(length(est)), mse = mse,
        ess = mean(ess),
        k_ess = B / mean(ess),                                # evaluations per effective sample
        k_mse = mse * B / truth[[h]][["var_h"]],              # the same quantity implied by the MSE
        acc = mean(df$acc), c_theory = df$c_theory[1], step = df$step[1], burnin = df$burnin[1],
        converged = df$burnin[1] < B / 2,                     # FALSE: pilot chains never reached R-hat < 1.1
        secs_per_run = mean(df$secs), ess_per_sec = mean(ess) / mean(df$secs),
        setup_secs = df$setup_secs[1], max_ratio = max(df$max_ratio),
        ess_own_rep1 = if (h == "mean") df$ess_mean[1] else NA, ess_coda_rep1 = if (h == "mean") df$ess_coda[1] else NA)
    }))
  }

  ## ---------------------------- 1.7 Self-checks ------------------------------
  # For every target/proposal pair: (a) the numerical maximisation of f~/g against the
  # closed form of M (where one exists), (b) c = M/Z with Z computed by numerical
  # integration against the value used, (c) the AR acceptance rate against 1/c,
  # (d) that the envelope is never exceeded, (e) a Kolmogorov-Smirnov test of the AR
  # output against the exact distribution. (c)-(e) are repeated N_CHECK times.
  N_CHECK <- 20
  run_checks <- function(pairs) {
    set.seed(SEED + SEED_OFFSET)
    cat(sprintf("---- Experiment 1 self-checks (%d repetitions each) ----\n", N_CHECK))
    for (p in pairs) {
      tg <- p$tg; pr <- p$pr
      logM_num <- log_sup_ratio(tg$lf, pr$lg, pr$range[1], pr$range[2])
      Z_num    <- integrate(function(x) exp(tg$lf(x)), if (tg$key == "gamma") 0 else -Inf, Inf, rel.tol = 1e-10)$value
      c_used   <- exp(pr$logM - tg$logZ)
      out <- replicate(N_CHECK, { o <- ar_sampler(tg, pr, 1e5)
        # (R's uniform generator has 32-bit resolution, so a few of the ~50,000 draws coincide
        #  exactly; the warning about ties is harmless and is suppressed)
        c(acc = o$acc, max_ratio = o$max_ratio, ks_p = suppressWarnings(ks.test(o$x, tg$cdf))$p.value) })
      cat(sprintf("%-20s g = %-13s log M: used %.6f, numerical %.6f | c: used %.4f, from numerical M and Z %.4f | 1/c = %.5f, acceptance %.5f (se %.5f) | largest ratio %.6f | KS p < 0.05: %d of %d\n",
                  tg$label, pr$label, pr$logM, logM_num, c_used, exp(logM_num) / Z_num, 1 / c_used,
                  mean(out["acc", ]), sd(out["acc", ]) / sqrt(N_CHECK), max(out["max_ratio", ]),
                  sum(out["ks_p", ] < 0.05), N_CHECK))
    }
    cat("---------------------\n")
  }

  ## ---------------------------- 1.8 Main -------------------------------------
  tg_gamma <- make_gamma(2.5); tg_t <- make_t(3); tg_mix <- lapply(DELTAS, make_mixture)
  targets  <- c(list(tg_gamma, tg_t), tg_mix); names(targets) <- sapply(targets, `[[`, "key")

  # target/proposal pairs for AR and IMH
  pairs <- c(lapply(LAMBDAS, function(l) list(tg = tg_gamma, pr = prop_exp(l, 2.5))),
             list(list(tg = tg_t, pr = prop_cauchy(3))),
             lapply(tg_mix, function(tg) list(tg = tg, pr = prop_normal(tg))))

  # All settings. The position of a setting in this list fixes its random seed, so new
  # settings must be added at the end.
  settings <- c(
    unlist(lapply(pairs, function(p) lapply(c("AR", "IMH"), function(s) c(p, sampler = s))), recursive = FALSE),
    lapply(unname(targets), function(tg) list(tg = tg, pr = NULL, sampler = "RWM")))

  sink(file.path(OUT_DIR, "exp1_checks.txt"), split = TRUE); run_checks(pairs); sink()

  t_start <- Sys.time()
  raw <- vector("list", length(settings)); res <- raw
  for (i in seq_along(settings)) {
    s <- settings[[i]]
    raw[[i]] <- run_setting(s$tg, s$pr, s$sampler, setting_id = i)
    res[[i]] <- summarise_setting(raw[[i]], s$tg)
    cat(sprintf("[Exp 1: %2d/%d] %-20s %-13s %-5s  burn-in=%6d  acc=%.3f  k(mean)=%9.2f  (%.0f s elapsed)\n",
                i, length(settings), s$tg$label, if (is.null(s$pr)) "" else s$pr$label, s$sampler,
                raw[[i]]$burnin[1], mean(raw[[i]]$acc), res[[i]]$k_ess[1],
                as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  }
  replicates <- do.call(rbind, raw); results <- do.call(rbind, res); rownames(results) <- NULL
  write.csv(replicates, file.path(OUT_DIR, "exp1_replicates.csv"), row.names = FALSE)
  write.csv(results,    file.path(OUT_DIR, "exp1_results.csv"),    row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "exp1_sessionInfo.txt"))

  ## ---------------------------- 1.9 LaTeX tables -----------------------------
  SAMPLERS <- c("AR", "IMH", "RWM")
  # one row of results: RWM does not use a proposal g
  pick <- function(tg, pr, s, h) {
    r <- results[results$target == tg$key & results$sampler == s & results$estimand == h, ]
    if (s %in% c("AR", "IMH")) r <- r[r$proposal == pr$label, ]
    r
  }
  flagged <- function(r) abs(r$bias / r$se_bias) > 2
  # Table 1 (cost k, main text) and Table 8 (MSE, appendix): one row per proposal
  two_block_table <- function(cell, file, note = NOTE_BOTH) {
    rows <- sapply(seq_along(pairs), function(i) {
      tg <- pairs[[i]]$tg; pr <- pairs[[i]]$pr
      first <- i == 1 || pairs[[i - 1]]$tg$key != tg$key      # first row of this target
      cells <- unlist(lapply(c("mean", "tail"), function(h) sapply(SAMPLERS, function(s)
        if (s == "RWM" && !first) "" else cell(pick(tg, pr, s, h)))))
      paste0(if (first) tg$tex else "", " & ", pr$tex, " & ", sprintf("%.2f", exp(pr$logM - tg$logZ)), " & ",
             paste(cells, collapse = " & "), " \\\\")
    })
    write_tex(c(
      "\\begin{tabular}{@{}ll r rrr rrr@{}}", "\\toprule",
      " & & & \\multicolumn{3}{c}{Mean $\\E_f[X]$} & \\multicolumn{3}{c}{Tail probability} \\\\",
      "\\cmidrule(lr){4-6}\\cmidrule(l){7-9}",
      paste0("Target & Proposal $g$ & $c$ & ", paste(rep(SAMPLERS, 2), collapse = " & "), " \\\\"),
      "\\midrule", rows, "\\bottomrule", "\\end{tabular}", foot(note)), OUT_DIR, file)
  }
  # Section mark: RWM is not geometrically ergodic on the polynomial tails of t_3, so its
  # ESS-based cost for the mean is not reliable (k from the MSE is much larger).
  unreliable <- function(r) r$sampler == "RWM" && r$target == "t" && r$estimand == "mean"
  two_block_table(function(r) paste0(fmt(r$k_ess), mark(flagged(r), !r$converged), if (unreliable(r)) SEC else ""),
                  "tab_exp1_cost.tex",
                  paste0(sub("\\.$", "", NOTE_BOTH), "; $\\S$: not reliable, see Section~\\ref{sec:shape-analysis}."))
  two_block_table(function(r) paste0(fmt_mse(r$mse), mark(flagged(r), !r$converged)), "tab_exp1_accuracy.tex")

  # Tables 4, 9 and 10 (appendix): one row per sampler
  ord <- unlist(lapply(seq_along(targets), function(j) {
    tg <- targets[[j]]
    c(unlist(lapply(pairs[sapply(pairs, function(p) p$tg$key == tg$key)], function(p)
      lapply(c("AR", "IMH"), function(s) list(tg = tg, pr = p$pr, s = s))), recursive = FALSE),
      list(list(tg = tg, pr = NULL, s = "RWM")))
  }), recursive = FALSE)
  long_table <- function(header, align, row_fun, note, file) {
    rows <- c(); last <- ""
    for (o in ord) {
      if (o$tg$key != last && last != "") rows <- c(rows, "\\addlinespace[2pt]")
      lab <- if (o$tg$key != last) o$tg$tex else ""; last <- o$tg$key
      rows <- c(rows, paste0(lab, " & ", if (is.null(o$pr)) "--" else o$pr$tex, " & ", row_fun(o), " \\\\"))
    }
    write_tex(c(sprintf("\\begin{tabular}{@{}%s@{}}", align), "\\toprule", header, "\\midrule", rows,
                "\\bottomrule", "\\end{tabular}", foot(note, "")), OUT_DIR, file)
  }
  for (h in c("mean", "tail"))
    long_table("Target & Proposal $g$ & Sampler & Bias & SE & Bias/SE & MSE & ESS & $k$ (ESS) & $k$ (MSE) \\\\",
               "l l l r r r r r r r",
               function(o) { r <- pick(o$tg, o$pr, o$s, h)
                 paste(paste0(o$s, if (!r$converged) DD else ""), fmt(r$bias, small = 0.01), fmt(r$se_bias, small = 0.01),
                       paste0(fx(r$bias / r$se_bias, 1), if (flagged(r)) DG else ""),
                       fmt(r$mse, small = 0.01), fmt(r$ess), fmt(r$k_ess), fmt(r$k_mse), sep = " & ") },
               NOTE_BOTH, sprintf("tab_exp1_full_%s.tex", h))
  long_table("Target & Proposal $g$ & Sampler & $c$ & Step $\\sigma$ & Acceptance & Burn-in $b$ & Setup (s) & Run (s) & ESS/s \\\\",
             "l l l r r r r r r r",
             function(o) { r <- pick(o$tg, o$pr, o$s, "mean")
               paste(o$s, if (is.na(r$c_theory)) "--" else sprintf("%.3f", r$c_theory),
                     if (is.na(r$step)) "--" else sprintf("%.2f", r$step), sprintf("%.3f", r$acc),
                     if (o$s == "AR") "--" else paste0(format(r$burnin, big.mark = ","), if (!r$converged) DD else ""),
                     sprintf("%.2f", r$setup_secs), sprintf("%.3f", r$secs_per_run), fmt(r$ess_per_sec), sep = " & ") },
             NOTE_CONV, "tab_exp1_setup.tex")

  ## ---------------------------- 1.10 Figure (Figure 3) -----------------------
  # Cost k for the mean on the two-mode target as the modes move apart
  st  <- STYLE[match(SAMPLERS, STYLE$sampler), ]
  mix <- subset(results, estimand == "mean" & grepl("^mix", target))
  mix$Delta <- as.numeric(sub("mix", "", mix$target))
  pdf(file.path(OUT_DIR, "fig_exp1_mixture.pdf"), width = 4.6, height = 3.5)
  par(mar = c(3.6, 4.4, 1.0, 1.2), mgp = c(2.4, 0.6, 0), tcl = -0.25)
  plot(NA, xlim = range(DELTAS), ylim = range(mix$k_ess), log = "y", axes = FALSE, xlab = "", ylab = "")
  yt <- axTicks(2); abline(h = yt, col = GRID, lwd = 0.8)
  axis(1, at = DELTAS, col = NA, col.ticks = INK, col.axis = INK, cex.axis = 0.85)
  axis(2, at = yt, labels = yt, las = 1, col = NA, col.ticks = INK, col.axis = INK, cex.axis = 0.85)
  title(xlab = "Half-distance between the modes (Delta)", line = 2.3, col.lab = INK, cex.lab = 0.9)
  title(ylab = "k  (lower is better)", line = 3.0, col.lab = INK, cex.lab = 0.9)
  for (j in seq_len(nrow(st))) {
    ds <- mix[mix$sampler == st$sampler[j], ]; ds <- ds[order(ds$Delta), ]
    lines(ds$Delta, ds$k_ess, col = st$col[j], lty = st$lty[j], lwd = 2)
    points(ds$Delta, ds$k_ess, col = st$col[j], pch = st$pch[j], cex = 1.1)
  }
  legend("topleft", legend = st$sampler, col = st$col, pch = st$pch, lty = st$lty, lwd = 2,
         bty = "n", cex = 0.8, text.col = INK, seg.len = 2.6)
  dev.off()

  ## ---------------------------- 1.11 Console report --------------------------
  show <- function(h) {
    cat(sprintf("\n==== Experiment 1 results: %s ====\n", if (h == "mean") "mean of X" else "tail probability P(X > q_0.95)"))
    print(format(subset(results, estimand == h)[, c("target", "proposal", "sampler", "bias", "se_bias", "mse", "ess",
          "k_ess", "k_mse", "acc", "c_theory", "step", "burnin", "converged", "ess_per_sec")], digits = 3), row.names = FALSE)
  }
  show("mean"); show("tail")
  cat(sprintf("\nLargest value of f~/(M g) over all AR runs: %.6f (must not exceed 1)\n", max(results$max_ratio, na.rm = TRUE)))
  cat(sprintf("Experiment 1 run time: %.1f minutes. Output written to '%s/'.\n\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins")), OUT_DIR))
}
if (RUN_EXP1) run_experiment1()


###############################################################################
# PART 2 -- Experiment 2: Dimension and Correlation (Section 3.2 and Appendix D)
#
# Target   : N_d(0, Sigma_rho), (Sigma_rho)_ij = rho^|i-j|
#            dimension sweep d = 1, 2, 5, 10, 20, 50 at rho = 0;
#            correlation sweep rho = 0, 0.5, 0.9, 0.99 at d = 5
# Samplers : AR, IMH (same proposal as AR), RWM, Gibbs; MALA (Appendix D)
# Estimands: the mean E_f[X_1] = 0 and the tail probability P_f(X_1 > 1.645) = 0.05
# Output   : CSV files, LaTeX tables and PDF figures in Exp2_output/
###############################################################################
run_experiment2 <- function() {

  ## ---------------------------- 2.0 Settings ---------------------------------
  DIMS         <- c(1, 2, 5, 10, 20, 50)    # dimension sweep (rho = 0)
  RHOS         <- c(0, 0.5, 0.9, 0.99)      # correlation sweep
  D_RHO        <- 5                         #   ... at this dimension
  INCLUDE_MALA <- TRUE                      # FALSE runs the four course samplers only
  Q95          <- qnorm(0.95)               # tail threshold: P(X_1 > Q95) = 0.05
  N_CORES      <- 1                         # > 1 runs replicates in parallel (Mac/Linux)
  OUT_DIR      <- "Exp2_output"
  dir.create(OUT_DIR, showWarnings = FALSE)

  ## ---------------------------- 2.1 Target -----------------------------------
  # Sigma_rho has a tridiagonal precision matrix Q, so x'Qx costs O(d).
  make_target <- function(d, rho) {
    if (d == 1) {
      qd <- 1; qo <- numeric(0); lam_max <- 1
    } else {
      qd <- c(1, rep(1 + rho^2, d - 2), 1) / (1 - rho^2)     # diagonal of Q
      qo <- rep(-rho / (1 - rho^2), d - 1)                   # off-diagonal of Q
      Sigma   <- rho^abs(outer(1:d, 1:d, "-"))
      lam_max <- max(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
    }
    s2 <- 1.44 * lam_max                                     # proposal variance (AR, IMH)
    list(d = d, rho = rho, qd = qd, qo = qo, s = sqrt(s2),
         c = exp((d / 2) * log(s2) - 0.5 * (d - 1) * log(1 - rho^2)))  # c = s^d (1-rho^2)^{-(d-1)/2}
  }
  # x'Qx for one point (vector) and for many points (rows of a matrix)
  quad_vec <- function(x, tg) {
    q <- sum(tg$qd * x * x)
    if (tg$d > 1) q <- q + 2 * sum(tg$qo * x[-1] * x[-tg$d])
    q
  }
  quad_mat <- function(Y, tg) {
    q <- drop((Y * Y) %*% tg$qd)
    if (tg$d > 1)
      q <- q + 2 * drop((Y[, -1, drop = FALSE] * Y[, -tg$d, drop = FALSE]) %*% tg$qo)
    q
  }
  log_ftilde      <- function(x, tg) -0.5 * quad_vec(x, tg)     # log f~(x), one evaluation
  grad_log_ftilde <- function(x, tg) {                          # gradient = -Qx
    g <- tg$qd * x
    if (tg$d > 1) g <- g + c(0, tg$qo * x[-tg$d]) + c(tg$qo * x[-1], 0)
    -g
  }
  # log{ f~(y) / (M g(y)) } for proposal g = N(0, s^2 I): always <= 0 when s^2 >= lambda_max
  log_w_mat <- function(Y, tg) -0.5 * (quad_mat(Y, tg) - rowSums(Y * Y) / tg$s^2)
  # one exact draw from the target (used only to start tuning runs in equilibrium)
  r_target <- function(tg) {
    x <- rnorm(tg$d)
    if (tg$d > 1) for (j in 2:tg$d) x[j] <- tg$rho * x[j - 1] + sqrt(1 - tg$rho^2) * x[j]
    x
  }

  ## ---------------------------- 2.2 Samplers ---------------------------------
  # Every sampler returns the first coordinate of its draws (x1) and its acceptance rate.
  # Acceptance-rejection (vectorised): exactly B proposals = B evaluations
  ar_sampler <- function(tg, B, chunk = 2e4) {
    x1 <- numeric(0); max_lw <- -Inf; done <- 0
    while (done < B) {
      m    <- min(chunk, B - done)
      Y    <- matrix(rnorm(m * tg$d, 0, tg$s), m, tg$d)
      lw   <- log_w_mat(Y, tg)
      keep <- log(runif(m)) <= lw
      x1   <- c(x1, Y[keep, 1])
      max_lw <- max(max_lw, lw); done <- done + m
    }
    list(x1 = x1, acc = length(x1) / B, max_ratio = exp(max_lw))  # max_ratio must be <= 1
  }
  # Independence MH with the same proposal as AR; proposals do not depend on the
  # current state, so they are generated in bulk.
  imh_sampler <- function(tg, n_iter, x0, chunk = 2e4) {
    y1 <- numeric(n_iter); lwy <- numeric(n_iter); done <- 0
    while (done < n_iter) {
      m   <- min(chunk, n_iter - done)
      Y   <- matrix(rnorm(m * tg$d, 0, tg$s), m, tg$d)
      idx <- done + seq_len(m)
      y1[idx] <- Y[, 1]; lwy[idx] <- log_w_mat(Y, tg); done <- done + m
    }
    logu <- log(runif(n_iter))
    x1 <- numeric(n_iter); cur1 <- x0[1]; n_acc <- 0
    cur_lw <- log_w_mat(matrix(x0, 1), tg)
    for (t in seq_len(n_iter)) {
      if (logu[t] <= lwy[t] - cur_lw) { cur1 <- y1[t]; cur_lw <- lwy[t]; n_acc <- n_acc + 1 }
      x1[t] <- cur1
    }
    list(x1 = x1, acc = n_acc / n_iter)
  }
  # Random-walk Metropolis
  rwm_sampler <- function(tg, n_iter, x0, sigma) {
    d <- tg$d; x <- x0; lx <- log_ftilde(x, tg)
    x1 <- numeric(n_iter); n_acc <- 0
    for (t in seq_len(n_iter)) {
      y  <- x + sigma * rnorm(d)
      ly <- log_ftilde(y, tg)
      if (log(runif(1)) <= ly - lx) { x <- y; lx <- ly; n_acc <- n_acc + 1 }
      x1[t] <- x[1]
    }
    list(x1 = x1, acc = n_acc / n_iter)
  }
  # Gibbs sampler with exact normal full conditionals (systematic scan)
  gibbs_sampler <- function(tg, n_sweep, x0) {
    d <- tg$d; rho <- tg$rho; stopifnot(d >= 2)
    a     <- rho / (1 + rho^2)                    # interior conditional mean coefficient
    s_in  <- sqrt((1 - rho^2) / (1 + rho^2))      # interior conditional sd
    s_end <- sqrt(1 - rho^2)                      # end-point conditional sd
    x <- x0; x1 <- numeric(n_sweep)
    for (t in seq_len(n_sweep)) {
      z    <- rnorm(d)
      x[1] <- rho * x[2] + s_end * z[1]
      if (d > 2) for (j in 2:(d - 1)) x[j] <- a * (x[j - 1] + x[j + 1]) + s_in * z[j]
      x[d] <- rho * x[d - 1] + s_end * z[d]
      x1[t] <- x[1]
    }
    list(x1 = x1, acc = 1)
  }
  # MALA (Appendix D)
  mala_sampler <- function(tg, n_iter, x0, h) {
    d <- tg$d; x <- x0; lx <- log_ftilde(x, tg)
    mx <- x + 0.5 * h * grad_log_ftilde(x, tg)
    x1 <- numeric(n_iter); n_acc <- 0
    for (t in seq_len(n_iter)) {
      y  <- mx + sqrt(h) * rnorm(d)
      ly <- log_ftilde(y, tg)
      my <- y + 0.5 * h * grad_log_ftilde(y, tg)
      log_alpha <- ly - lx + (sum((y - mx)^2) - sum((x - my)^2)) / (2 * h)
      if (log(runif(1)) <= log_alpha) { x <- y; lx <- ly; mx <- my; n_acc <- n_acc + 1 }
      x1[t] <- x[1]
    }
    list(x1 = x1, acc = n_acc / n_iter)
  }
  # Budget units charged per draw (a Gibbs sweep costs about one evaluation of f~ here;
  # MALA: one for f~ and one for its gradient)
  UNITS_PER_DRAW <- c(AR = 1, IMH = 1, RWM = 1, Gibbs = 1, MALA = 2)
  run_mcmc <- function(sampler, tg, n_draws, x0, step = NA) {
    switch(sampler,
           IMH   = imh_sampler(tg, n_draws, x0),
           RWM   = rwm_sampler(tg, n_draws, x0, step),
           Gibbs = gibbs_sampler(tg, n_draws, x0),
           MALA  = mala_sampler(tg, n_draws, x0, step))
  }

  ## ---------------------------- 2.3 Setup: step size and burn-in --------------
  # Step size. At rho = 0 the theoretically optimal value is used (2.38/sqrt(d) for RWM,
  # 1.65^2 d^(-1/3) for MALA). For rho > 0 it is found by bisection (on the log scale)
  # so that the acceptance rate hits 0.234 (RWM) or 0.574 (MALA).
  choose_step <- function(sampler, tg, n_iter = 5000, n_bisect = 12) {
    if (!sampler %in% c("RWM", "MALA")) return(NA_real_)
    theory <- if (sampler == "RWM") 2.38 / sqrt(tg$d) else 1.65^2 * tg$d^(-1/3)
    if (tg$rho == 0) return(theory)
    target <- if (sampler == "RWM") 0.234 else 0.574
    lo <- theory / 1e4; hi <- 2 * theory
    for (i in seq_len(n_bisect)) {
      mid <- sqrt(lo * hi)
      acc <- run_mcmc(sampler, tg, n_iter, r_target(tg), mid)$acc
      if (acc > target) lo <- mid else hi <- mid
    }
    sqrt(lo * hi)
  }
  # Pilot chains start at +-(3, ..., 3) and +-(3, -3, 3, ...); R-hat uses the draws of X_1
  choose_burnin <- function(sampler, tg, step, n_pilot, max_b) {
    if (sampler == "AR") return(c(b = 0, rhat_end = NA))
    d <- tg$d; alt <- rep(c(1, -1), length.out = d)
    starts <- list(rep(3, d), rep(-3, d), 3 * alt, -3 * alt)
    X <- t(sapply(starts, function(x0) run_mcmc(sampler, tg, n_pilot, x0, step)$x1))
    burnin_from_pilots(X, max_b)
  }

  ## ---------------------------- 2.4 One setting ------------------------------
  # Runs one sampler on one target: setup (timed), then N_REP replicates at budget B.
  run_setting <- function(sweep, d, rho, sampler, setting_id) {
    tg <- make_target(d, rho); x0 <- rep(3, d)
    n_draws <- B / UNITS_PER_DRAW[[sampler]]                 # draws affordable with budget B
    set.seed(SEED + 1e5 * setting_id)
    t_setup <- system.time({
      step <- choose_step(sampler, tg)
      bi   <- choose_burnin(sampler, tg, step, n_pilot = min(N_PILOT, n_draws / 2),
                            max_b = n_draws / 2)
    })[["elapsed"]]
    b <- bi[["b"]]
    one_rep <- function(r) {
      set.seed(SEED + 1e5 * setting_id + r)
      secs <- system.time(
        out <- if (sampler == "AR") ar_sampler(tg, B) else run_mcmc(sampler, tg, n_draws, x0, step)
      )[["elapsed"]]
      z <- if (sampler == "AR") out$x1 else out$x1[(b + 1):n_draws]   # discard burn-in
      n <- length(z); tail_ind <- as.numeric(z > Q95)
      data.frame(rep = r, n_kept = n, acc = out$acc, secs = secs,
                 est_mean = if (n > 0) mean(z) else NA,
                 est_tail = if (n > 0) mean(tail_ind) else NA,
                 ess_mean = if (sampler == "AR") n else ess_geyer(z),
                 ess_tail = if (sampler == "AR") n else ess_geyer(tail_ind),
                 ess_coda = if (r == 1 && sampler != "AR" && requireNamespace("coda", quietly = TRUE))
                   as.numeric(coda::effectiveSize(z)) else NA,
                 max_ratio = if (sampler == "AR") out$max_ratio else NA)
    }
    reps <- if (N_CORES > 1) parallel::mclapply(seq_len(N_REP), one_rep, mc.cores = N_CORES)
            else lapply(seq_len(N_REP), one_rep)
    reps <- do.call(rbind, reps)
    cbind(sweep = sweep, d = d, rho = rho, sampler = sampler, c_theory = tg$c,
          step = step, burnin = b, rhat_pilot_end = bi[["rhat_end"]], setup_secs = t_setup, reps)
  }

  ## ---------------------------- 2.5 Summaries --------------------------------
  # True values and variances of the two estimands (X_1 ~ N(0,1) for every d and rho)
  TRUTH <- list(mean = c(theta = 0,    var_h = 1),
                tail = c(theta = 0.05, var_h = 0.05 * 0.95))
  summarise_setting <- function(df) {
    out <- lapply(names(TRUTH), function(h) {
      est <- df[[paste0("est_", h)]]; ess <- df[[paste0("ess_", h)]]
      th  <- TRUTH[[h]][["theta"]]; ok <- !is.na(est)
      ess[is.na(ess)] <- 0                                    # no usable draws: ESS = 0
      mse <- mean((est[ok] - th)^2)                           # AR at rho = 0.99: runs with a draw only
      data.frame(
        sweep = df$sweep[1], d = df$d[1], rho = df$rho[1], sampler = df$sampler[1],
        estimand = h, n_fail = sum(!ok),
        bias = mean(est[ok]) - th, se_bias = sd(est[ok]) / sqrt(sum(ok)), mse = mse,
        ess = mean(ess),
        k_ess = B / mean(ess),                                # evaluations per effective sample
        k_mse = mse * B / TRUTH[[h]][["var_h"]],              # the same quantity implied by the MSE
        acc = mean(df$acc), c_theory = df$c_theory[1], step = df$step[1], burnin = df$burnin[1],
        secs_per_run = mean(df$secs), ess_per_sec = mean(ess) / mean(df$secs),
        setup_secs = df$setup_secs[1],
        ess_own_rep1 = if (h == "mean") df$ess_mean[1] else NA, ess_coda_rep1 = if (h == "mean") df$ess_coda[1] else NA)
    })
    do.call(rbind, out)
  }

  ## ---------------------------- 2.6 Self-checks ------------------------------
  # Each check is repeated N_CHECK times with different random numbers, so that it
  # does not depend on a single lucky or unlucky run.
  N_CHECK <- 20
  run_checks <- function() {
    set.seed(SEED)
    ms <- function(v) sprintf("%.3f (sd %.3f)", mean(v), sd(v))
    cat(sprintf("---- Experiment 2 self-checks (%d repetitions each) ----\n", N_CHECK))
    # (i) ESS estimator on AR(1) series with known tau = (1+phi)/(1-phi)
    for (phi in c(0, 0.5, 0.9)) {
      tau <- t(replicate(N_CHECK, {
        z <- as.numeric(stats::filter(rnorm(1e5), phi, method = "recursive"))   # AR(1) series
        c(own = 1e5 / ess_geyer(z),
          coda = if (requireNamespace("coda", quietly = TRUE)) 1e5 / as.numeric(coda::effectiveSize(z)) else NA)
      }))
      cat(sprintf("ESS check, AR(1) phi=%.1f: tau true = %6.3f | own estimator %s | coda %s\n",
                  phi, (1 + phi) / (1 - phi), ms(tau[, "own"]), ms(tau[, "coda"])))
    }
    # (ii) Gibbs at d = 2 against tau = (1+rho^2)/(1-rho^2)
    for (rho in c(0.5, 0.9)) {
      tau <- replicate(N_CHECK, 1e5 / ess_geyer(gibbs_sampler(make_target(2, rho), 1e5, c(0, 0))$x1))
      cat(sprintf("Gibbs check, d=2 rho=%.1f: tau true = %6.3f | estimated %s\n",
                  rho, (1 + rho^2) / (1 - rho^2), ms(tau)))
    }
    # (iii) AR: acceptance rate against 1/c, envelope never exceeded, output is N(0,1)
    for (d in c(1, 5)) {
      tg  <- make_target(d, 0.5 * (d > 1))
      out <- replicate(N_CHECK, { o <- ar_sampler(tg, 1e5)
        c(acc = o$acc, max_ratio = o$max_ratio, ks_p = ks.test(o$x1, "pnorm")$p.value) })
      cat(sprintf("AR check, d=%d rho=%.1f: 1/c = %.5f | acceptance %.5f (se %.5f) | largest ratio %.3f | KS p-values below 0.05: %d of %d\n",
                  d, tg$rho, 1 / tg$c, mean(out["acc", ]), sd(out["acc", ]) / sqrt(N_CHECK),
                  max(out["max_ratio", ]), sum(out["ks_p", ] < 0.05), N_CHECK))
    }
    # (iv) Gelman-Rubin statistic for four chains of independent N(0,1) draws, which need
    #      no burn-in: the course version (on running means) against the version used here
    gr_course <- function(psi) {
      n <- ncol(psi); W <- mean(apply(psi, 1, var)); Bv <- n * var(rowMeans(psi))
      (W * (n - 1) / n + Bv / n) / W
    }
    rh <- replicate(N_CHECK, {
      X   <- matrix(rnorm(4 * 5e4), 4)
      psi <- t(apply(X, 1, cumsum)) / matrix(seq_len(5e4), 4, 5e4, byrow = TRUE)
      c(running = gr_course(psi), draws = rhat_second_half(X, 200))
    })
    cat(sprintf("R-hat check, independent draws: on running means after 50,000 iterations %s, below %.1f in %d of %d | on the draws after 200 iterations %s, below %.1f in %d of %d\n",
                ms(rh["running", ]), RHAT_TOL, sum(rh["running", ] < RHAT_TOL), N_CHECK,
                ms(rh["draws", ]), RHAT_TOL, sum(rh["draws", ] < RHAT_TOL), N_CHECK))
    cat("---------------------\n")
  }

  ## ---------------------------- 2.7 Main -------------------------------------
  samplers <- c("AR", "IMH", "RWM", "Gibbs", if (INCLUDE_MALA) "MALA")
  make_grid <- function(s) rbind(
    expand.grid(sweep = "dimension",   d = DIMS,  rho = 0,    sampler = s, stringsAsFactors = FALSE),
    expand.grid(sweep = "correlation", d = D_RHO, rho = RHOS, sampler = s, stringsAsFactors = FALSE))
  grid <- make_grid(c("AR", "IMH", "RWM", "Gibbs"))
  grid <- grid[!(grid$sampler == "Gibbs" & grid$d == 1), ]      # Gibbs needs d >= 2
  # MALA settings are appended at the end. The row number of a setting fixes its random
  # seed, so the results of the four course samplers do not change when MALA is added.
  if (INCLUDE_MALA) grid <- rbind(grid, make_grid("MALA"))

  sink(file.path(OUT_DIR, "exp2_checks.txt"), split = TRUE); run_checks(); sink()

  t_start <- Sys.time()
  raw <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    raw[[i]] <- run_setting(g$sweep, g$d, g$rho, g$sampler, setting_id = i)
    cat(sprintf("[Exp 2: %2d/%d] %-11s d=%2d rho=%.2f %-5s  burn-in=%6d  acc=%.3f  ESS(mean)=%9.1f  (%.0f s elapsed)\n",
                i, nrow(grid), g$sweep, g$d, g$rho, g$sampler, raw[[i]]$burnin[1], mean(raw[[i]]$acc),
                mean(raw[[i]]$ess_mean, na.rm = TRUE), as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  }
  replicates <- do.call(rbind, raw)
  results    <- do.call(rbind, lapply(raw, summarise_setting))
  rownames(results) <- NULL

  # Growth of cost with dimension (mean estimand): log k = a + beta log d  and  log k = a + d log gamma
  growth <- do.call(rbind, lapply(samplers, function(s) {
    df <- subset(results, sweep == "dimension" & estimand == "mean" & sampler == s & is.finite(k_ess))
    fp <- lm(log(k_ess) ~ log(d), data = df); fe <- lm(log(k_ess) ~ d, data = df)
    data.frame(sampler = s,
               beta = coef(fp)[[2]], beta_lo = confint(fp)[2, 1], beta_hi = confint(fp)[2, 2],
               r2_polynomial = summary(fp)$r.squared,
               log_gamma = coef(fe)[[2]], log_gamma_lo = confint(fe)[2, 1], log_gamma_hi = confint(fe)[2, 2],
               r2_exponential = summary(fe)$r.squared)
  }))

  write.csv(replicates, file.path(OUT_DIR, "exp2_replicates.csv"), row.names = FALSE)
  write.csv(results,    file.path(OUT_DIR, "exp2_results.csv"),    row.names = FALSE)
  write.csv(growth,     file.path(OUT_DIR, "exp2_growth.csv"),     row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "exp2_sessionInfo.txt"))

  ## ---------------------------- 2.8 Figures (Figures 1 and 2) ----------------
  short_num <- function(y)                      # 1500 -> "1.5k", 2e6 -> "2M"
    ifelse(y >= 1e6, paste0(y / 1e6, "M"), ifelse(y >= 1e3, paste0(y / 1e3, "k"), as.character(y)))
  panel <- function(df, xvar, yvar, ylab, main, x_at, x_lab, logx, theory = NULL) {
    df <- df[is.finite(df[[yvar]]) & df[[yvar]] > 0, ]
    yr <- range(c(df[[yvar]], if (!is.null(theory)) theory$y[is.finite(theory$y)]))
    plot(NA, xlim = range(x_at), ylim = yr, log = if (logx) "xy" else "y", axes = FALSE,
         xlab = "", ylab = "", main = "")
    yt <- axTicks(2); abline(h = yt, col = GRID, lwd = 0.8)
    axis(1, at = x_at, labels = x_lab, col = NA, col.ticks = INK, col.axis = INK, cex.axis = 0.85)
    axis(2, at = yt, labels = short_num(yt), las = 1, col = NA, col.ticks = INK, col.axis = INK, cex.axis = 0.85)
    title(main = main, adj = 0, font.main = 1, cex.main = 0.95, col.main = INK)
    title(xlab = if (xvar == "d") "Dimension d" else expression("Correlation " * rho), line = 2.3, col.lab = INK, cex.lab = 0.9)
    title(ylab = ylab, line = 3.0, col.lab = INK, cex.lab = 0.9)
    if (!is.null(theory)) lines(theory$x, theory$y, col = "#8a8a8a", lty = 3, lwd = 1.2)
    for (s in unique(df$sampler)) {
      st <- STYLE[STYLE$sampler == s, ]; ds <- df[df$sampler == s, ]; ds <- ds[order(ds[[xvar]]), ]
      xs <- if (xvar == "d") ds$d else match(ds$rho, RHOS)    # rho values evenly spaced, as categories
      lines(xs, ds[[yvar]], col = st$col, lty = st$lty, lwd = 2)
      points(xs, ds[[yvar]], col = st$col, pch = st$pch, cex = if (st$pch == 18) 1.5 else 1.1)
    }
  }
  add_legend <- function(present, theory_label = NULL, pos = "topleft") {
    st <- STYLE[STYLE$sampler %in% present, ]
    legend(pos, legend = c(st$sampler, theory_label), col = c(st$col, if (!is.null(theory_label)) "#8a8a8a"),
           pch = c(st$pch, if (!is.null(theory_label)) NA), lty = c(st$lty, if (!is.null(theory_label)) 3),
           lwd = c(rep(2, nrow(st)), if (!is.null(theory_label)) 1.2), bty = "n", cex = 0.8, text.col = INK, seg.len = 2.6)
  }
  # cost k for the mean: (a) dimension sweep, (b) correlation sweep. The figure for the
  # main text shows the four course samplers; fig_exp2_cost_mala.pdf adds MALA.
  res_mean <- subset(results, estimand == "mean")
  COURSE   <- c("AR", "IMH", "RWM", "Gibbs")
  cost_figure <- function(file, who) {
    pdf(file.path(OUT_DIR, file), width = 8.4, height = 3.7)
    par(mfrow = c(1, 2), mar = c(3.6, 4.4, 2.2, 1.2), mgp = c(2.4, 0.6, 0), tcl = -0.25)
    dfd <- subset(res_mean, sweep == "dimension" & sampler %in% who)
    dfr <- subset(res_mean, sweep == "correlation" & sampler %in% who)
    panel(dfd, "d", "k_ess", "k  (lower is better)", "(a) Increasing dimension", DIMS, DIMS, TRUE,
          list(x = DIMS, y = 1.2^DIMS))
    add_legend(unique(dfd$sampler), "theory for AR: k = c", "topleft")
    panel(dfr, "rho", "k_ess", "k  (lower is better)", "(b) Increasing correlation", seq_along(RHOS), RHOS, FALSE,
          list(x = seq_along(RHOS), y = sapply(RHOS, function(r) make_target(D_RHO, r)$c)))
    dev.off()
  }
  cost_figure("fig_exp2_cost.pdf", COURSE)
  if (INCLUDE_MALA) cost_figure("fig_exp2_cost_mala.pdf", samplers)

  ## ---------------------------- 2.9 LaTeX tables -----------------------------
  r  <- results[results$sampler %in% COURSE, ]
  r$z    <- r$bias / r$se_bias
  r$flag <- is.finite(r$z) & abs(r$z) > 2 & r$n_fail == 0          # AR at rho = 0.99 is left aside
  r$nc   <- r$sampler != "AR" & r$burnin >= B / 2
  ra <- results; ra$z <- ra$bias / ra$se_bias
  ordr <- function(df) df[order(match(df$sweep, c("dimension", "correlation")), df$d, df$rho, match(df$sampler, COURSE)), ]

  # Table 2 (main text): MSE of both estimands
  acc_rows <- function(sweep, key, values, label) sapply(values, function(v) {
    cells <- c()
    for (h in c("mean", "tail")) for (s in COURSE) {
      x <- r[r$sweep == sweep & r[[key]] == v & r$estimand == h & r$sampler == s, ]
      cells <- c(cells, if (nrow(x) == 0) "--" else paste0(fmt_mse(x$mse), mark(x$flag, x$nc), if (x$n_fail > 0) SEC else ""))
    }
    paste0("$", label, " = ", v, "$ & ", paste(cells, collapse = " & "), " \\\\")
  })
  write_tex(c("\\begin{tabular}{@{}l rrrr rrrr@{}}", "\\toprule",
    " & \\multicolumn{4}{c}{Mean $\\E_f[X_1]$} & \\multicolumn{4}{c}{Tail probability $P_f(X_1 > 1.645)$} \\\\",
    "\\cmidrule(lr){2-5}\\cmidrule(l){6-9}", paste0(" & ", paste(rep(COURSE, 2), collapse = " & "), " \\\\"), "\\midrule",
    "\\multicolumn{9}{@{}l}{\\emph{Increasing Dimension} ($\\rho = 0$)} \\\\", acc_rows("dimension", "d", DIMS, "d"),
    "\\addlinespace",
    "\\multicolumn{9}{@{}l}{\\emph{Increasing Correlation} ($d = 5$)} \\\\", acc_rows("correlation", "rho", RHOS, "\\rho"),
    "\\bottomrule", "\\end{tabular}",
    foot(sprintf("$\\dagger$: the bias exceeds two standard errors; $\\ddagger$: the pilot chains did not converge; $\\S$: over the %d of %d runs that returned a draw.",
                 N_REP - max(r$n_fail), N_REP))),
    OUT_DIR, "tab_exp2_accuracy.tex")

  # Table 3 (main text): growth of k with dimension
  theory <- c(AR = "$\\log\\gamma = \\log 1.2 = 0.182$", IMH = "--", RWM = "$\\beta = 1$", Gibbs = "$\\beta = 0$")
  gr <- growth[growth$sampler %in% COURSE, ]
  write_tex(c(
    "\\begin{tabular}{@{}l c c c c l@{}}", "\\toprule",
    " & \\multicolumn{2}{c}{Polynomial fit, $k \\propto d^{\\beta}$} & \\multicolumn{2}{c}{Exponential fit, $k \\propto \\gamma^{d}$} & \\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
    "Sampler & $\\hat\\beta$ (95\\% CI) & $R^2$ & $\\log\\hat\\gamma$ (95\\% CI) & $R^2$ & Theory \\\\", "\\midrule",
    sapply(seq_len(nrow(gr)), function(i) with(gr[i, ], sprintf(
      "%s & %s (%s, %s) & %.3f & %s (%s, %s) & %.3f & %s \\\\",
      sampler, fx(beta, 2), fx(beta_lo, 2), fx(beta_hi, 2), r2_polynomial,
      fx(log_gamma), fx(log_gamma_lo), fx(log_gamma_hi), r2_exponential, theory[[sampler]]))),
    "\\bottomrule", "\\end{tabular}"), OUT_DIR, "tab_exp2_growth.tex")

  # Tables 11 and 12 (appendix): full results for each estimand
  full_table <- function(h) {
    df <- ordr(r[r$estimand == h, ]); rows <- c(); last <- ""
    for (i in seq_len(nrow(df))) { x <- df[i, ]; grp <- paste(x$sweep, x$d, x$rho)
      if (grp != last && last != "") rows <- c(rows, "\\addlinespace[2pt]")
      lab <- if (grp != last) sprintf("%d & %s", x$d, format(x$rho)) else " & "; last <- grp
      rows <- c(rows, sprintf("%s & %s%s & %s & %s & %s%s & %s & %s & %s & %s & %s \\\\", lab, x$sampler, if (x$nc) DD else "",
        fmt(x$bias, small = 0.01), fmt(x$se_bias, small = 0.01), if (is.finite(x$z)) fx(x$z, 1) else "--", if (x$flag) DG else "",
        fmt(x$mse, small = 0.01), fmt(x$ess), fmt(x$k_ess), fmt(x$k_mse), if (x$n_fail > 0) sprintf("%d", x$n_fail) else "")) }
    c("\\begin{tabular}{@{}r r l r r r r r r r r@{}}", "\\toprule",
      "$d$ & $\\rho$ & Sampler & Bias & SE & Bias/SE & MSE & ESS & $k$ (ESS) & $k$ (MSE) & Failed \\\\", "\\midrule",
      rows, "\\bottomrule", "\\end{tabular}", foot(NOTE_BOTH, ""))
  }
  write_tex(full_table("mean"), OUT_DIR, "tab_exp2_full_mean.tex")
  write_tex(full_table("tail"), OUT_DIR, "tab_exp2_full_tail.tex")

  # Table 5 (appendix): setup and run time of every setting
  st <- ordr(r[r$estimand == "mean", ]); rows <- c(); last <- ""
  for (i in seq_len(nrow(st))) { x <- st[i, ]; grp <- paste(x$sweep, x$d, x$rho)
    if (grp != last && last != "") rows <- c(rows, "\\addlinespace[2pt]")
    lab <- if (grp != last) sprintf("%d & %s & %s", x$d, format(x$rho), fmt(x$c_theory)) else " & & "; last <- grp
    rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", lab, x$sampler,
      if (is.na(x$step)) "--" else sprintf("%.3f", x$step), fmt(x$acc),
      if (x$sampler == "AR") "--" else paste0(format(x$burnin, big.mark = ","), if (x$nc) DD else ""), sprintf("%.2f", x$setup_secs),
      sprintf("%.3f", x$secs_per_run), fmt(x$ess_per_sec))) }
  write_tex(c("\\begin{tabular}{@{}r r r l r r r r r r@{}}", "\\toprule",
    "$d$ & $\\rho$ & $c$ & Sampler & Step size & Acceptance & Burn-in $b$ & Setup (s) & Run (s) & ESS/s \\\\", "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}", foot(NOTE_CONV, "")), OUT_DIR, "tab_exp2_setup.tex")

  # Table 6 (appendix D): MALA against RWM
  if (INCLUDE_MALA) {
    get <- function(sweep, key, v, s, h, col) ra[ra$sweep == sweep & ra[[key]] == v & ra$sampler == s & ra$estimand == h, col]
    mala_rows <- function(sweep, key, values, label) sapply(values, function(v) {
      mse <- sapply(c("mean", "tail"), function(h) paste0(fmt_mse(get(sweep, key, v, "MALA", h, "mse")),
                                                          if (abs(get(sweep, key, v, "MALA", h, "z")) > 2) DG else ""))
      rw  <- sapply(c("mean", "tail"), function(h) paste0(fmt(get(sweep, key, v, "RWM", h, "k_ess")),
                                                          if (get(sweep, key, v, "RWM", h, "burnin") >= B / 2) DD else ""))
      sprintf("$%s = %s$ & %.3f & %.3f & %s & %s & %s & %s & %s & %s & %s \\\\", label, v,
              get(sweep, key, v, "MALA", "mean", "step"), get(sweep, key, v, "MALA", "mean", "acc"),
              paste0(format(get(sweep, key, v, "MALA", "mean", "burnin"), big.mark = ","),   # MALA's cap is n_draws / 2 = B / 4
                     if (get(sweep, key, v, "MALA", "mean", "burnin") >= B / (2 * UNITS_PER_DRAW[["MALA"]])) DD else ""),
              rw[1], fmt(get(sweep, key, v, "MALA", "mean", "k_ess")), rw[2], fmt(get(sweep, key, v, "MALA", "tail", "k_ess")),
              mse[1], mse[2])
    })
    write_tex(c("\\begin{tabular}{@{}l rrr rr rr rr@{}}", "\\toprule",
      " & \\multicolumn{3}{c}{MALA setup} & \\multicolumn{2}{c}{$k$, mean} & \\multicolumn{2}{c}{$k$, tail probability} & \\multicolumn{2}{c}{$\\widehat{\\MSE}$ of MALA} \\\\",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(l){9-10}",
      " & Step $h$ & Acceptance & Burn-in $b$ & RWM & MALA & RWM & MALA & Mean & Tail \\\\", "\\midrule",
      "\\multicolumn{10}{@{}l}{\\emph{Increasing Dimension} ($\\rho = 0$)} \\\\", mala_rows("dimension", "d", DIMS, "d"),
      "\\addlinespace",
      "\\multicolumn{10}{@{}l}{\\emph{Increasing Correlation} ($d = 5$)} \\\\", mala_rows("correlation", "rho", RHOS, "\\rho"),
      "\\bottomrule", "\\end{tabular}",
      foot("$\\dagger$: the bias exceeds two standard errors; $\\ddagger$: the pilot chains did not converge.")),
      OUT_DIR, "tab_exp2_mala.tex")
  }

  ## ---------------------------- 2.10 Console report --------------------------
  cat("\n==== Experiment 2 results: mean of X_1 ====\n")
  print(format(res_mean[, c("sweep", "d", "rho", "sampler", "n_fail", "bias", "se_bias", "mse", "ess",
                            "k_ess", "k_mse", "acc", "burnin", "ess_per_sec")], digits = 3), row.names = FALSE)
  cat("\n==== Experiment 2 results: tail probability P(X_1 > 1.645) ====\n")
  print(format(subset(results, estimand == "tail")[, c("sweep", "d", "rho", "sampler", "n_fail", "bias", "se_bias",
                                                       "mse", "ess", "k_ess", "k_mse")], digits = 3), row.names = FALSE)
  cat("\n==== Growth of k with dimension (rho = 0) ====\n"); print(format(growth, digits = 3), row.names = FALSE)
  cat(sprintf("\nExperiment 2 run time: %.1f minutes. Output written to '%s/'.\n\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins")), OUT_DIR))
}
if (RUN_EXP2) run_experiment2()


###############################################################################
# PART 3 -- Real-Data Example: the Challenger O-Ring Data (Appendix E)
#
# Data     : 23 shuttle launches before the Challenger accident (Dalal, Fowlkes and
#            Hoadley, 1989). y = 1 if at least one O-ring showed thermal distress.
# Model    : y_i ~ Bernoulli(p_i),  logit(p_i) = beta0 + beta1 * temp_i,
#            prior beta0, beta1 ~ N(0, 100^2) independently (nearly flat)
# Target   : the posterior of (beta0, beta1), in two parameterisations of the same model:
#              raw          : temperature in degrees F         (strongly correlated)
#              standardised : temperature centred and scaled   (nearly uncorrelated)
# Samplers : AR and IMH with the same proposal g, a bivariate t fitted to the posterior;
#            RWM tuned to an acceptance rate of 0.234 (lectures, multidimensional case)
# Estimands: posterior mean of beta1, and posterior mean of p(31), the probability of
#            O-ring distress at 31 F, the temperature of the Challenger launch
# Output   : CSV files and one LaTeX table in ExpReal_output/
###############################################################################
run_real_data <- function() {

  ## ---------------------------- 3.0 Settings ---------------------------------
  ACC_TARGET  <- 0.234        # acceptance rate RWM is tuned to
  PRIOR_SD    <- 100          # prior standard deviation of beta0 and beta1
  NU          <- 4            # degrees of freedom of the t proposal
  S_G         <- 1.3          # proposal scale: S_G^2 x (Laplace covariance); c is flat near this value
  M_MARGIN    <- 1.01         # the numerically found bound M is inflated by 1%
  SEED_OFFSET <- 2e7          # keeps these seeds distinct from those of Parts 1 and 2
  OUT_DIR     <- "ExpReal_output"
  dir.create(OUT_DIR, showWarnings = FALSE)

  ## ---------------------------- 3.1 Data -------------------------------------
  # In launch order. Source: Dalal, Fowlkes and Hoadley (1989), JASA 84, 945-957.
  temp <- c(66, 70, 69, 68, 67, 72, 73, 70, 57, 63, 70, 78, 67, 53, 67, 75, 70, 81, 76, 79, 75, 76, 58)
  fail <- c( 0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  1)
  T_LAUNCH <- 31

  ## ---------------------------- 3.2 Posterior --------------------------------
  log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))   # log(1 + e^x) without overflow
  # theta = (theta1, theta2) are the coefficients for the predictor (temp - shift) / scale,
  # so that beta1 = theta2 / scale and beta0 = theta1 - theta2 * shift / scale.
  make_posterior <- function(name, shift, scale) {
    X <- cbind(1, (temp - shift) / scale)
    A <- rbind(c(1, -shift / scale), c(0, 1 / scale))          # beta = A theta
    P <- crossprod(A) / PRIOR_SD^2                             # prior precision for theta
    lf <- function(th) {                                       # log f~, vectorised over rows
      th  <- matrix(th, ncol = 2)
      eta <- th %*% t(X)
      drop(eta %*% fail) - rowSums(log1pexp(eta)) - 0.5 * rowSums((th %*% P) * th)
    }
    x2 <- X[, 2]; sum_y <- sum(fail); sum_yx <- sum(fail * x2)
    lf1 <- function(th) {                                      # the same for a single point (faster in loops)
      eta <- th[1] + th[2] * x2
      th[1] * sum_y + th[2] * sum_yx - sum(log1pexp(eta)) -
        0.5 * (P[1, 1] * th[1]^2 + 2 * P[1, 2] * th[1] * th[2] + P[2, 2] * th[2]^2)
    }
    # posterior mode by Newton's method (the log-posterior is concave)
    th <- c(0, 0)
    for (it in 1:100) {
      p  <- plogis(drop(X %*% th))
      gr <- drop(crossprod(X, fail - p)) - drop(P %*% th)
      H  <- -crossprod(X, X * (p * (1 - p))) - P
      step <- solve(H, gr); th <- th - step
      if (max(abs(step)) < 1e-12) break
    }
    S <- solve(-H)                                             # Laplace approximation to the covariance
    list(name = name, shift = shift, scale = scale, lf = lf, lf1 = lf1, mode = th, S = S, L = t(chol(S)),
         h = list(beta1 = function(th) matrix(th, ncol = 2)[, 2] / scale,
                  p31   = function(th) { th <- matrix(th, ncol = 2); plogis(th[, 1] + th[, 2] * (T_LAUNCH - shift) / scale) }))
  }

  # True posterior quantities by numerical integration on a fine grid, in coordinates
  # standardised by the Laplace approximation (z = L^{-1}(theta - mode)): the estimands,
  # their variances, the posterior standard deviations and correlation, and log Z.
  grid_points <- function(tg, half_width, n) {
    g <- seq(-half_width, half_width, length.out = n)
    Z <- as.matrix(expand.grid(g, g))
    list(TH = sweep(Z %*% t(tg$L), 2, tg$mode, "+"), Z = Z, cell = (g[2] - g[1])^2 * det(tg$L))
  }
  true_values <- function(tg) {
    G  <- grid_points(tg, 10, 601)
    lf <- tg$lf(G$TH); w <- exp(lf - max(lf))
    out <- lapply(tg$h, function(h) { v <- h(G$TH); m <- sum(w * v) / sum(w)
      c(theta = m, var_h = sum(w * v^2) / sum(w) - m^2) })
    wn <- w / sum(w); mu <- colSums(wn * G$TH)
    V  <- crossprod(G$TH * sqrt(wn)) - tcrossprod(mu)                    # true posterior covariance
    edge <- max(w[abs(G$Z[, 1]) > 9.9 | abs(G$Z[, 2]) > 9.9]) / max(w)   # mass left at the edge
    c(out, list(sd = sqrt(diag(V)), corr = cov2cor(V)[1, 2],
                logZ = log(sum(w) * G$cell) + max(lf), edge = edge))
  }

  ## ---------------------------- 3.3 Proposal for AR and IMH ------------------
  # Bivariate t with NU degrees of freedom, centred at the mode, scale matrix S_G^2 x S.
  make_proposal <- function(tg) {
    Lg <- S_G * tg$L
    lconst <- lgamma((NU + 2) / 2) - lgamma(NU / 2) - log(NU * pi) - log(det(Lg))
    lg <- function(th) {
      z <- forwardsolve(Lg, t(sweep(matrix(th, ncol = 2), 2, tg$mode)))
      lconst - (NU + 2) / 2 * log1p(colSums(z^2) / NU)
    }
    r <- function(n) sweep(matrix(rnorm(2 * n), n) %*% t(Lg) / sqrt(rchisq(n, NU) / NU), 2, tg$mode, "+")
    # M = sup f~/g: best point of a fine grid, refined with a numerical optimiser
    G <- grid_points(tg, 12, 801)
    v <- tg$lf(G$TH) - lg(G$TH); i <- which.max(v)
    opt <- optim(G$TH[i, ], function(th) -(tg$lf(th) - lg(th)), method = "BFGS", control = list(reltol = 1e-14))
    list(lg = lg, r = r, logM = max(-opt$value, v[i]) + log(M_MARGIN),
         edge = max(v[abs(G$Z[, 1]) > 11.9 | abs(G$Z[, 2]) > 11.9]) - max(v))
  }

  ## ---------------------------- 3.4 Samplers ---------------------------------
  ar_sampler <- function(tg, pr, B) {                      # exactly B proposals = B evaluations
    y  <- pr$r(B)
    lw <- tg$lf(y) - pr$lg(y) - pr$logM                    # log{ f~(y) / (M g(y)) }, must be <= 0
    keep <- log(runif(B)) <= lw
    list(x = y[keep, , drop = FALSE], acc = mean(keep), max_ratio = exp(max(lw)))
  }
  imh_sampler <- function(tg, pr, n_iter, x0) {
    y <- pr$r(n_iter); lwy <- tg$lf(y) - pr$lg(y); logu <- log(runif(n_iter))
    x <- matrix(0, n_iter, 2); cur <- x0; cur_lw <- tg$lf(x0) - pr$lg(x0); n_acc <- 0
    for (t in seq_len(n_iter)) {
      if (logu[t] <= lwy[t] - cur_lw) { cur <- y[t, ]; cur_lw <- lwy[t]; n_acc <- n_acc + 1 }
      x[t, ] <- cur
    }
    list(x = x, acc = n_acc / n_iter)
  }
  rwm_sampler <- function(tg, n_iter, x0, sigma) {
    z <- matrix(sigma * rnorm(2 * n_iter), n_iter); logu <- log(runif(n_iter))
    x <- matrix(0, n_iter, 2); cur <- x0; cur_lf <- tg$lf1(cur); n_acc <- 0
    for (t in seq_len(n_iter)) {
      y <- cur + z[t, ]; ly <- tg$lf1(y)
      if (logu[t] <= ly - cur_lf) { cur <- y; cur_lf <- ly; n_acc <- n_acc + 1 }
      x[t, ] <- cur
    }
    list(x = x, acc = n_acc / n_iter)
  }
  run_mcmc <- function(sampler, tg, pr, n_iter, x0, step)
    if (sampler == "IMH") imh_sampler(tg, pr, n_iter, x0) else rwm_sampler(tg, n_iter, x0, step)

  ## ---------------------------- 3.5 Step size and burn-in --------------------
  # RWM step: bisection on log(sigma) towards ACC_TARGET, in short runs started at the mode
  choose_step <- function(tg, n_iter = 5000, n_bisect = 16) {
    lo <- 1e-4; hi <- 1e2
    for (i in seq_len(n_bisect)) {
      mid <- sqrt(lo * hi)
      if (rwm_sampler(tg, n_iter, tg$mode, mid)$acc > ACC_TARGET) lo <- mid else hi <- mid
    }
    sqrt(lo * hi)
  }
  # Pilot chains start at mode + L z for z = (3, 3), (-3, -3), (3, -3), (-3, 3), where
  # L L' is the Laplace covariance; R-hat uses the draws of beta1.
  start_point <- function(tg, z) drop(tg$mode + tg$L %*% z)
  choose_burnin <- function(sampler, tg, pr, step) {
    starts <- list(c(3, 3), c(-3, -3), c(3, -3), c(-3, 3))
    X <- t(sapply(starts, function(z)
      tg$h$beta1(run_mcmc(sampler, tg, pr, N_PILOT, start_point(tg, z), step)$x)))
    burnin_from_pilots(X, max_b = B / 2)[["b"]]
  }

  ## ---------------------------- 3.6 Run --------------------------------------
  targets <- list(make_posterior("Raw", 0, 1),
                  make_posterior("Standardised", mean(temp), sd(temp)))
  truth <- lapply(targets, true_values)
  props <- lapply(targets, make_proposal)

  cat("---- Challenger posterior and proposal ----\n")
  for (j in seq_along(targets)) {
    tg <- targets[[j]]; tv <- truth[[j]]
    cat(sprintf("%-12s mode (%.4f, %.4f) | posterior sd (%.4f, %.4f), correlation %.3f | Laplace sd (%.4f, %.4f), correlation %.3f\n",
                tg$name, tg$mode[1], tg$mode[2], tv$sd[1], tv$sd[2], tv$corr,
                sqrt(tg$S[1, 1]), sqrt(tg$S[2, 2]), cov2cor(tg$S)[1, 2]))
    cat(sprintf("%-12s E[beta1] = %.6f, E[p(31)] = %.6f | c = %.4f | grid mass at edge %.1e\n",
                "", tv$beta1[["theta"]], tv$p31[["theta"]], exp(props[[j]]$logM - tv$logZ), tv$edge))
  }
  cat(sprintf("Largest log(f~/g) at the edge of the search grid, relative to its maximum: %.1f\n\n",
              max(sapply(props, `[[`, "edge"))))

  settings <- expand.grid(sampler = c("AR", "IMH", "RWM"), param = seq_along(targets), stringsAsFactors = FALSE)
  raw <- list(); t_start <- Sys.time()
  for (i in seq_len(nrow(settings))) {
    s <- settings$sampler[i]; j <- settings$param[i]; tg <- targets[[j]]; pr <- props[[j]]
    set.seed(SEED + SEED_OFFSET + 1e5 * i)
    t_setup <- system.time({
      step <- if (s == "RWM") choose_step(tg) else NA
      b    <- if (s == "AR") 0 else choose_burnin(s, tg, pr, step)
    })[["elapsed"]]
    reps <- do.call(rbind, lapply(seq_len(N_REP), function(r) {
      set.seed(SEED + SEED_OFFSET + 1e5 * i + r)
      secs <- system.time(out <- if (s == "AR") ar_sampler(tg, pr, B)
                                 else run_mcmc(s, tg, pr, B, start_point(tg, c(3, 3)), step))[["elapsed"]]
      x <- if (s == "AR") out$x else out$x[(b + 1):B, , drop = FALSE]
      est <- sapply(tg$h, function(h) mean(h(x)))
      ess <- sapply(tg$h, function(h) if (s == "AR") nrow(x) else ess_geyer(h(x)))
      data.frame(param = tg$name, sampler = s, rep = r, acc = out$acc, secs = secs,
                 est_beta1 = est[["beta1"]], est_p31 = est[["p31"]],
                 ess_beta1 = ess[["beta1"]], ess_p31 = ess[["p31"]],
                 max_ratio = if (s == "AR") out$max_ratio else NA)
    }))
    raw[[i]] <- cbind(reps, step = step, burnin = b, setup_secs = t_setup)
    cat(sprintf("[Real: %d/%d] %-12s %-3s  burn-in = %6d  acceptance = %.3f  (%.0f s elapsed)\n", i, nrow(settings),
                tg$name, s, b, mean(reps$acc), as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  }
  replicates <- do.call(rbind, raw)

  ## ---------------------------- 3.7 Summary ----------------------------------
  results <- do.call(rbind, lapply(raw, function(df) {
    j  <- match(df$param[1], sapply(targets, `[[`, "name"))
    do.call(rbind, lapply(c("beta1", "p31"), function(h) {
      est <- df[[paste0("est_", h)]]; ess <- df[[paste0("ess_", h)]]; ess[is.na(ess)] <- 0
      th  <- truth[[j]][[h]][["theta"]]; vh <- truth[[j]][[h]][["var_h"]]; mse <- mean((est - th)^2)
      data.frame(param = df$param[1], sampler = df$sampler[1], estimand = h, truth = th,
                 bias = mean(est) - th, se_bias = sd(est) / sqrt(length(est)), mse = mse,
                 ess = mean(ess), k_ess = B / mean(ess), k_mse = mse * B / vh,
                 acc = mean(df$acc), c_theory = if (df$sampler[1] == "RWM") NA else exp(props[[j]]$logM - truth[[j]]$logZ),
                 step = df$step[1], burnin = df$burnin[1], converged = df$burnin[1] < B / 2,
                 secs_per_run = mean(df$secs), ess_per_sec = mean(ess) / mean(df$secs),
                 setup_secs = df$setup_secs[1], max_ratio = max(df$max_ratio))
    }))
  }))
  posterior <- do.call(rbind, lapply(seq_along(targets), function(j) data.frame(
    param = targets[[j]]$name, mode1 = targets[[j]]$mode[1], mode2 = targets[[j]]$mode[2],
    sd1 = truth[[j]]$sd[1], sd2 = truth[[j]]$sd[2], corr = truth[[j]]$corr,
    laplace_sd1 = sqrt(targets[[j]]$S[1, 1]), laplace_sd2 = sqrt(targets[[j]]$S[2, 2]),
    laplace_corr = cov2cor(targets[[j]]$S)[1, 2],
    E_beta1 = truth[[j]]$beta1[["theta"]], E_p31 = truth[[j]]$p31[["theta"]],
    c = exp(props[[j]]$logM - truth[[j]]$logZ))))
  write.csv(replicates, file.path(OUT_DIR, "real_replicates.csv"), row.names = FALSE)
  write.csv(results,    file.path(OUT_DIR, "real_results.csv"),    row.names = FALSE)
  write.csv(posterior,  file.path(OUT_DIR, "real_posterior.csv"),  row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "real_sessionInfo.txt"))

  ## ---------------------------- 3.8 LaTeX table (Table 7) --------------------
  cell <- function(p, s, h, col) results[results$param == p & results$sampler == s & results$estimand == h, col]
  rows <- c()
  for (p in sapply(targets, `[[`, "name")) {
    if (length(rows) > 0) rows <- c(rows, "\\addlinespace[2pt]")
    for (s in c("AR", "IMH", "RWM")) {
      nc  <- !cell(p, s, "beta1", "converged")
      mse <- sapply(c("beta1", "p31"), function(h) paste0(fmt_mse(cell(p, s, h, "mse")),
               if (abs(cell(p, s, h, "bias") / cell(p, s, h, "se_bias")) > 2) DG else ""))
      rows <- c(rows, sprintf("%s & %s%s & %.3f & %s & %s & %s & %s & %s & %s \\\\",
        if (s == "AR") p else "", s, if (nc) DD else "", cell(p, s, "beta1", "acc"),
        if (s == "AR") "--" else format(cell(p, s, "beta1", "burnin"), big.mark = ","),
        fmt(cell(p, s, "beta1", "k_ess")), fmt(cell(p, s, "p31", "k_ess")), mse[1], mse[2],
        fmt(cell(p, s, "beta1", "ess_per_sec"))))
    }
  }
  write_tex(c("\\begin{tabular}{@{}ll r r rr rr r@{}}", "\\toprule",
    " & & & & \\multicolumn{2}{c}{$k$} & \\multicolumn{2}{c}{$\\widehat{\\MSE}$} & \\\\",
    "\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
    "Temperature & Sampler & Acceptance & Burn-in $b$ & $\\beta_1$ & $p(31)$ & $\\beta_1$ & $p(31)$ & ESS/s \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}", foot(NOTE_BOTH)), OUT_DIR, "tab_real.tex")

  ## ---------------------------- 3.9 Console report ---------------------------
  print(format(results[, c("param", "sampler", "estimand", "truth", "bias", "se_bias", "mse", "k_ess", "k_mse",
                           "acc", "c_theory", "step", "burnin", "converged", "ess_per_sec")], digits = 3), row.names = FALSE)
  ar <- results[results$sampler == "AR" & results$estimand == "beta1", ]
  cat(sprintf("\nAR acceptance against 1/c: %s\n", paste(sprintf("%.4f vs %.4f", ar$acc, 1 / ar$c_theory), collapse = "; ")))
  cat(sprintf("Largest value of f~/(M g) over all AR proposals: %.6f (must not exceed 1)\n", max(results$max_ratio, na.rm = TRUE)))
  cat(sprintf("Real-data run time: %.1f minutes. Output written to '%s/'.\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins")), OUT_DIR))
}
if (RUN_REAL) run_real_data()
