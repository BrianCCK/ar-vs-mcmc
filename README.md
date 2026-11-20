# A Comparative Study of Acceptance-Rejection Sampling and MCMC Algorithms

Acceptance–rejection (AR) sampling returns exact, independent draws, while Markov chain Monte Carlo (MCMC) returns correlated ones. This project asks which gives the more accurate Monte Carlo estimate for a **fixed number of evaluations of an unnormalised target density**, and how the answer changes with the shape, dimension and correlation of the target.

Author: Chi Kong Brian CHAN

## Approach

To compare exact and correlated samplers on one scale, every sampler is measured by

$$k = \frac{B}{\text{ESS}} \quad \text{(evaluations of } \tilde f \text{ per effective sample)},$$

which is the envelope constant $c$ for AR and the autocorrelation time $\tau$ for MCMC. For the same budget $B$, $\text{MSE} \approx \sigma_h^2 k / B$, so the sampler with the smaller $k$ is the more accurate one.

**Samplers:** acceptance–rejection (AR), the independence sampler (IMH), random-walk Metropolis (RWM) and the Gibbs sampler, plus the Metropolis-adjusted Langevin algorithm (MALA) as a gradient-based extension.

**Experiments** (budget $B = 10^5$ evaluations, 50 independent replicates per setting):

| | Targets | Varies |
|---|---|---|
| Experiment 1 | Gamma(2.5, 1), Student $t_3$, two-mode normal mixture | shape, envelope quality, distance between modes |
| Experiment 2 | $N_d(0, \Sigma_\rho)$ with $(\Sigma_\rho)_{ij} = \rho^{\lvert i-j \rvert}$ | dimension $d = 1, \dots, 50$; correlation $\rho = 0, \dots, 0.99$ |
| Real data | Bayesian logistic regression on the Challenger O-ring data | raw vs standardised predictor |

ESS is estimated with Geyer's initial positive sequence estimator (checked against `coda`), and the burn-in is chosen from four dispersed pilot chains with the Gelman–Rubin statistic.

## Main findings

- **Low dimension, good envelope: AR wins.** In one dimension AR needs 1.2 to 2.2 evaluations per effective sample, against 4.4 to 6.7 for RWM, and its draws are exact with no tuning or burn-in.
- **IMH is not uniformly cheaper than AR.** With the same proposal its cost can reach $2c - 1$; it cost more than AR in 10 of 16 comparisons in Experiment 1.
- **Dimension favours MCMC.** AR's cost grows exponentially ($c = 1.2^d$ for a proposal 20% wider than the target, about 9,100 at $d = 50$), while RWM's grows roughly linearly ($k \propto d^{0.91}$). RWM becomes cheaper beyond about 23 dimensions.
- **Correlation slows every sampler whose proposal ignores it.** At $d = 5$ and $\rho = 0.99$, AR's constant ($3.4 \times 10^5$) exceeds the whole budget, and Gibbs, the cheapest sampler here, still needs about 400 evaluations per effective sample.
- **MCMC can fail silently.** RWM on a well-separated mixture, and on the strongly correlated Challenger posterior, returned full-length chains with a normal-looking acceptance rate but biased estimates; only comparing several chains revealed the problem.
- **Real data confirms the pattern.** With a proposal fitted to the posterior, AR and IMH need fewer than two evaluations per effective sample on the Challenger data.

The best sampler depends on what is known about the target. A good envelope favours AR, exact full conditionals favour Gibbs, a gradient favours MALA, and RWM is the fallback when only $\tilde f$ is available.

## Repository contents

```
report/   AR_vs_MCMC.pdf           Report
          (LaTeX source)
code/     Experiments.R            Experiments
data/     challenger_oring.csv          Challenger O-ring data
```

The Challenger data are also written into the script, so `data/challenger_oring.csv` is provided for convenience only.

## Running the code

Requirements: R (tested with R 4.5.1). Only base R is needed; the package `coda` is optional and used only as an independent check of the ESS estimator.

```bash
cd code
Rscript MATH3428_project.R
```

The full run takes about half an hour. To run only some parts, set the flags at the top of the script:

```r
RUN_EXP1 <- TRUE   # Experiment 1: one-dimensional targets
RUN_EXP2 <- TRUE   # Experiment 2: dimension and correlation, with MALA
RUN_REAL <- TRUE   # Real-data example: Challenger O-ring data
```

Each part writes its results to its own folder in the working directory:

| Folder | Contents |
|---|---|
| `Exp1_output/` | results and replicates (CSV), self-check log, LaTeX tables, mixture figure |
| `Exp2_output/` | results, replicates and growth-law fits (CSV), self-check log, LaTeX tables, cost figures |
| `ExpReal_output/` | results, replicates and true posterior summaries (CSV), LaTeX table |

**Reproducibility.** All random seeds are fixed (derived from 3428), and each setting has its own seed, so a part gives the same results whether or not the other parts are run. Rerunning the script reproduces every table and figure in the report exactly, apart from the run times, which depend on the machine.

## Data source

Dalal, S. R., Fowlkes, E. B. and Hoadley, B. (1989). Risk analysis of the space shuttle: pre-Challenger prediction of failure. *Journal of the American Statistical Association*, 84(408), 945–957.
