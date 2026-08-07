# Summary: Pierce et al. — Chapter 11: Estimating Animal Abundance

## Full Citation
Pierce, B.L., R.R. Lopez, and N.J. Silvy. 2012. Chapter 11: Estimating Animal Abundance. Pages 284–321 *in* Silvy, N.J. (ed.), *The Wildlife Techniques Manual*, 7th ed. The Wildlife Society, Bethesda, MD.

---

## Overview
A broad-audience survey chapter covering the full range of animal abundance estimation methods, written for undergraduate wildlife students and practicing biologists. Emphasizes conceptual foundations, survey design, and method selection over mathematical derivation. Primary reference used in this project for grounding the Chapman mark-recapture estimator in standard wildlife practice.

---

## Key Methods Covered
- Definitions: population, census, population estimate, population estimator, closed/open population, detection probability
- Precision vs. accuracy (rifle-target analogy from Overton and Davis 1969)
- Survey design considerations (spatial extent, experimental units, sample units, sampling intensity, power analysis)
- **Indices** (relative abundance only; recommend detection probability accompany all index data)
- **Census / total counts** (drive counts, aerial photography, spot mapping)
- **Counts on fixed-area plots** (strip counts, point counts)
- **Mark-recapture methods** (Lincoln-Petersen, Chapman modification)

---

## Lincoln-Petersen and Chapman Estimators

The Lincoln-Petersen estimator for a single mark-recapture occasion is:

$$\hat{N}_{LP} = \frac{n_1 n_2}{m}$$

where $n_1$ = number marked and released, $n_2$ = number captured in the second sample, and $m$ = number of marked individuals in the second sample.

**Chapman (1951)** proposed a nearly unbiased modification:

$$\hat{N}_C = \frac{(n_1+1)(n_2+1)}{m+1} - 1$$

with variance:

$$s^2_{\hat{N}} = \frac{(n_1+1)(n_2+1)(n_1-m)(n_2-m)}{(m+1)^2(m+2)}$$

The chapter recommends using the Chapman estimator as the standard because it reduces the positive bias of the Lincoln-Petersen estimator when sample sizes are small.

---

## Assumptions of Mark-Recapture Methods (Chapter 11)
1. The population is closed (no births, deaths, immigration, or emigration during the sampling period).
2. All individuals have the same probability of capture in both samples.
3. Marks are not lost and are correctly identified.
4. The second sample is drawn randomly from the population.

---

## Relevance to This Project
- Provides the methodological context and original citation basis for the Chapman estimator used in all backwater population estimates.
- The formula `N̂ = (M+1)(C+1)/(R+1) - 1` used in `BWMarkRecaptureEstimates.R` and `PopulationMonitoring.qmd` is the Chapman (1951) estimator presented in this chapter.
- The notation differs by convention only: $n_1$ = M (marked), $n_2$ = C (captured), $m$ = R (recaptured).
