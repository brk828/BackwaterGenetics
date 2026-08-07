# Summary: Sangnawakij et al. 2026 — Two-Way Capture-Recapture Methods with Emphasis on the Bootstrap

## Full Citation
Sangnawakij, P., R. Lerdsuwansri, P. Pijitrattana, P. Schlattmann, A. Maruotti, and D. Böhning. 2026. Two-way capture-recapture methods with emphasis on the bootstrap. *Statistical Methods & Applications* 35:1–22. https://doi.org/10.1007/s10260-026-00840-5

---

## Overview
Focuses on dual-system (two-source, two-occasion) capture-recapture estimation using the Chapman estimator, with emphasis on bootstrapping as an alternative to asymptotic normal confidence intervals. The paper is the primary reference supporting the use of bootstrap confidence intervals in the backwater population monitoring project.

Key finding: **the imputed bootstrap (Bootstrap I) provides coverage probabilities closest to the nominal level**, especially for small sample sizes. Simple non-parametric bootstrap (Bootstrap III, without imputation) and asymptotic normal intervals perform adequately for large sample sizes but undercover when sample sizes are small or when the recapture cell ($n_{11}$) is small.

---

## Two-Source Setting
Data take the form of a 2×2 contingency table:

|  | Source 2: Observed | Source 2: Not observed |
|---|---|---|
| **Source 1: Observed** | $n_{11}$ | $n_{10}$ |
| **Source 1: Not observed** | $n_{01}$ | $n_{00}$ (unknown) |

In the mark-recapture context: Source 1 = marking period, Source 2 = capture period, so $n_{11}$ = recaptures (R), $n_{10}$ = marked but not recaptured (M − R), $n_{01}$ = captured but not previously marked (C − R), and $n_{00}$ = unmarked and uncaptured (unknown).

---

## Estimators

### Lincoln-Petersen
$$\hat{N}_{LP} = \frac{n_{1\cdot} \, n_{\cdot 1}}{n_{11}}$$

Undefined if $n_{11} = 0$; positively biased for small sample sizes.

### Chapman (1951)
$$\hat{N}_C = \frac{(n_{1\cdot}+1)(n_{\cdot 1}+1)}{n_{11}+1} - 1$$

Nearly unbiased; exactly unbiased when $n_{1\cdot} + n_{\cdot 1} \geq N$. Reduces to the Lincoln-Petersen ratio estimate when sample sizes are large. The asymptotic variance is:

$$\widehat{Var}(\hat{N}_C) = \frac{(n_{1\cdot}+1)(n_{\cdot 1}+1) \, n_{10} \, n_{01}}{(n_{11}+1)^2 (n_{11}+2)}$$

Normal-approximation (Wald) confidence intervals based on this variance undercover when $n_{11}$ is small.

### Bias-Corrected Chapman
$$\hat{N}_{BC} = \frac{\hat{N}_C}{1 - \hat{\lambda}}$$

where $\hat{\lambda} = \exp\{-(n_{1\cdot}+1)(n_{\cdot 1}+1)/\hat{N}_C\}$. Provides a small-sample correction for residual downward bias in the Chapman estimator. Effect is negligible for moderate to large sample sizes.

---

## Bootstrap Methods Compared

### Bootstrap I — Imputed Bootstrap (recommended)
1. Compute event probabilities from $\hat{N}_C$ including the imputed missing cell.
2. Resample from a Multinomial($\hat{N}_C$, $\hat{p}$) distribution.
3. Compute Chapman estimate from each resample.
4. Use percentile CI from B bootstrap replicates.

Preferred because it accounts for uncertainty in the unobserved cell and yields the best coverage probabilities.

### Bootstrap II — Double Bootstrap with Imputation
An extended version of Bootstrap I that incorporates additional uncertainty by re-estimating the missing cell in each bootstrap iteration. More computationally intensive; performs similarly to Bootstrap I.

### Bootstrap III — Simple (Non-Parametric) Bootstrap Without Imputation
1. Resample from observed data only (Multinomial($n$, $\tilde{p}$) where $n = n_{11} + n_{10} + n_{01}$).
2. Compute Chapman estimate from each resample.
3. Use percentile CI from B replicates.

Simpler; performs well for large samples but undercoverers for small $n_{11}$.

---

## Simulation Findings
- Imputed bootstrap (Bootstrap I) consistently achieves coverage probabilities closest to the nominal 95% level across sample sizes and source dependence scenarios.
- Simple bootstrap (Bootstrap III) and normal-approximation intervals perform adequately only for larger $n_{11}$.
- Even moderate positive dependence between sources can degrade coverage, but the imputed bootstrap is the most robust.

---

## Relevance to This Project
The `recapr::ciChapman()` function with `method = "boot"` uses a bootstrap approach that resamples recaptures from a Binomial distribution with parameters $n = C$ and $p = R/C$, which is a simplified parametric bootstrap. This is conceptually similar to Bootstrap III (without imputation), in that the unobserved cell is not imputed. Sangnawakij et al. (2026) indicate that the imputed bootstrap (Bootstrap I) is preferred for small sample sizes. For the backwater monitoring data — where recapture numbers $R$ sometimes approach the minimum threshold of 4 — this distinction is relevant: reported confidence intervals should be interpreted with the understanding that coverage may be slightly below nominal, particularly in years with very few recaptures.

The Chapman point estimate $\hat{N}_C$ produced by `ciChapman()` correctly applies the `−1` in the formula: for M, C, and R values, $\hat{N} = (M+1)(C+1)/(R+1) - 1$.

**Special case note:** When $R = C$ (all captured fish were previously marked), the Chapman formula simplifies algebraically to $\hat{N} = M$, which coincidentally equals the simple ratio $MC/R = M$. This is not a formula error — the $-1$ is present and the formula gives the correct result; it is a mathematical identity that holds only when $R = C$.
