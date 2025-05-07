# 📊 Robust Dependence Modeling and Copula Simulation in R

This repository provides an end-to-end framework for **multivariate dependence modeling** using **t-distributions** and **copula theory**. The code is implemented in R with LaTeX integration (`.Rnw`) and demonstrates:

- Robust estimation of multivariate t-distribution parameters using profile likelihood
- Sector-wise and decade-wise comparison of financial return distributions
- Simulation and calibration of Gaussian and t-copulas
- Visualization of empirical vs. simulated copulas
- Estimation of tail dependence coefficients

> 📘 Part of EECS 545 coursework at the University of Michigan.

---

## 📁 Project Files

| File | Description |
|------|-------------|
| `Robust-Dependence-Modeling-and-Copula-Simulation.Rnw` | Main R + LaTeX source |
| `Robust-Dependence-Modeling-and-Copula-Simulation.pdf` | Compiled report |
| `stock_returns.csv` | Daily return data of tech & construction sector stocks |

---

## 🧩 Sector Profile Likelihood Estimation

- Fitted separate multivariate t-distributions for **tech** (AAPL, AMD, INTC, IBM) and **construction** (AA, CAT, F, MMM) sectors
- Estimated ν using profile log-likelihood
- Constructed 95% confidence intervals from chi-squared thresholds

**Result:**
- **Tech sector**: more heavy-tailed (lower ν), higher volatility
- **Construction sector**: more stable (higher ν), lighter tails

📈 Profile likelihood plots showed **non-overlapping CIs**, confirming significantly different distributions.

---

## 🧭 Problem 1b: Decade-wise Analysis (1980s–2010s)

Estimated ν and 95% CIs for both sectors per decade:

| Decade     | Sector       | ν̂ MLE | CI Lower | CI Upper |
|------------|--------------|--------|----------|----------|
| 1980–1989  | Tech         | 4.50   | 4.40     | 4.50     |
| 1980–1989  | Construction | 4.50   | 4.21     | 4.50     |
| 1990–1999  | Tech         | 4.45   | 4.08     | 4.50     |
| 1990–1999  | Construction | 4.50   | 4.40     | 4.50     |
| 2000–2009  | Tech         | 3.23   | 2.98     | 3.50     |
| 2000–2009  | Construction | 4.02   | 3.70     | 4.39     |
| 2010–2019  | Tech         | 3.89   | 3.57     | 4.25     |
| 2010–2019  | Construction | 4.50   | 4.17     | 4.50     |

**Insights:**
- Tech sector ν dropped significantly in 2000s, indicating heavier tails.
- Construction sector remained more stable across decades.

---

## 🔁 Problem 2: Copula Simulation and Calibration

- Implemented `sim_Gauss_copula()` and `sim_T_copula()` for generating Gaussian and t-copula samples
- Calibrated copulas for IBM & CAT using empirical Kendall’s τ
- Compared empirical vs. simulated copulas for various ν values

```r
sim_T_copula <- function(n, Sig, nu){
  Z = MASS::mvrnorm(n, mu = rep(0, ncol(Sig)), Sigma = Sig)
  T = Z / sqrt(rchisq(n, df=nu) / nu)
  U = pt(T, df=nu)
  return(U)
}
```
---
## 🧩 Tail Dependence Estimation
Plotted empirical vs. t-copula lower tail dependence:


```r
plot(alpha_seq, lambda_empirical, type="l", col="blue", lwd=2,
     xlab=expression(alpha), ylab=expression(lambda[L](alpha)))
lines(alpha_seq, lambda_t_copula, col="red", lwd=2)
legend("topright", legend=c("Empirical", "t-Copula (ν=5)"), col=c("blue", "red"), lty=1)
```
