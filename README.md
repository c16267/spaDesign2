<img src="man/figures/logo.png" align="right" height="150" alt="spaDesign2 logo" />

# spaDesign2

> **Pilot-driven power analysis and sample-size planning for multisample spatial transcriptomics.**

![version](https://img.shields.io/badge/version-1.0.1-5E8CB6)
![R](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-5E8CB6)
![license](https://img.shields.io/badge/license-GPL%20(%E2%89%A5%203)-5E8CB6)
![platforms](https://img.shields.io/badge/platforms-Visium%20%7C%20Visium%20HD%20%7C%20Stereo--seq-5E8CB6)

`spaDesign2` answers the first question of every comparative spatial-transcriptomics
study — **how many samples per group do I need?** — by learning a cohort-level
generative model from a small **pilot** dataset and running a full
*generate → recover → test* Monte-Carlo pipeline. Power is reported per endpoint and
per group sample size (*K*), and summarized into a concrete sample-size recommendation.

Unlike single-cell power tools, `spaDesign2` is built for the structure that defines
spatial data: it models **domain composition**, **spatial geometry**, and **spatially
correlated expression** jointly, re-discovers spatial domains in every synthetic
replicate, and evaluates the *same* tests an analyst would actually run — so the
reported power reflects the real analysis pipeline, not an idealized oracle.


---

## Why spaDesign2

- **Pilot-based, single-pilot friendly.** Learns the generator from real pilot
  sections; designs with as little as one section per group are supported through
  automatic pseudo-replication.
- **A three-layer generative model** that respects spatial biology:
  composition (baseline-anchored Binomial–logit), geometry
  (Fisher–Gaussian kernel mixture, FGKMM), and expression
  (Gaussian process with nearest-neighbor approximation, NNGP).
- **Honest *generate → recover → test* power.** Spatial domains are recovered in each
  replicate with a pilot-guided BANKSY extension (**pBANKSY**), so clustering error is
  *inside* the power estimate rather than assumed away.
- **Two complementary endpoints**, both with a minimum-effect **TREAT** margin:
  a spatially-aware differential-expression test (**SaLFC**) and a baseline-anchored
  compositional test (**LOR**).
- **Calibrated, monotone reporting.** Per-group sample sizes are recommended from
  Monte-Carlo estimates smoothed with shape-constrained additive models (SCAM),
  guaranteeing power curves that are monotone in *K*.
- **Cross-platform.** Validated on **10x Visium**, **Visium HD**, and **Stereo-seq**.

---

## Installation

`spaDesign2` links to compiled code (`Rcpp` / `RcppArmadillo`) and uses a few
Bioconductor packages, so install those first.

```r
# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "SingleCellExperiment", "SpatialExperiment", "spatialLIBD"))

# spaDesign2 (development version from GitHub)
# install.packages("remotes")
remotes::install_github("c16267/spaDesign2")    # replace with the repo path
```

Or from a local source tarball:

```r
install.packages("spaDesign2_1.0.1.tar.gz", repos = NULL, type = "source")
```

**Requirements.** R ≥ 4.1.0 and a C++ toolchain (for `Rcpp`/`RcppArmadillo`).
Heavy steps are parallelized via `n_cores`.

---

## Method at a glance

A `spaDesign2` analysis has two phases: **fit** the generator once from pilot data,
then **sweep** sample sizes and effect sizes to map power.

**1. Cohort-level generative model (learned from the pilot).**

| Layer | What it captures | Model |
|-------|------------------|-------|
| Composition | Per-domain abundance, case vs. control | Baseline-anchored **Binomial–logit** GLM (intercept β₀, effect β₁) |
| Geometry | Shape and placement of spatial domains | **Fisher–Gaussian kernel mixture (FGKMM)**, pooled across samples |
| Expression | Domain means + spatial covariance + between-sample variance | **Gaussian process (NNGP)** with per-gene parameters |

**2. Generate → recover → test (per Monte-Carlo replicate).**

```
simulateSpaDesign2()  ─►  pBANSKY()  ─►  [rearrangeSyntheticToPilot()]  ─►  SaLFC() / LOR()
  synthetic cohort        recover         morph onto pilot geometry          endpoint test
  at effect size δ, K     domains (ĥd)     (SaLFC only)                       + TREAT margin τ
```

Repeating this across a grid of (`K`, effect size) yields a rejection-rate surface;
the **null row** (effect = 0) returns the empirical type-I error / FDP.

**3. Summarize → recommend.** `summary()` and `plotPowerCurve*()` fit a monotone SCAM
(monotone-increasing in *K*) to the surface and report the **minimum *K*** reaching a
target power (default 0.8) for each effect size.

---

## Quick start

The package ships a downsampled human-DLPFC pilot (10x Visium) so the full workflow
runs out of the box.

```r
library(spaDesign2)

## 1 — Pilot data  ->  spaDesign2 object ---------------------------------------
data(mini_pilot_data_list)
obj <- createSpaDesign2Object(mini_pilot_data_list)   # auto pseudo-reps if K = 1

## 2 — Feature selection (domain-informative + stable/null genes) --------------
obj <- featureSelection(obj,        max_num_gene = 20)
obj <- featureSelectionStable(obj,  max_num_gene = 50)

## 3 — Labeling / null / spike-in gene sets ------------------------------------
sets <- makeCustomGeneSets(
  obj,
  G_svg_base      = obj@params_expression$top_genes,
  G_stable        = obj@params_expression$stable_genes,
  target_domain   = "WM",
  reference_domain = "Layer6"
)

## 4 — Fit the three-layer generative model ------------------------------------
obj <- estimateCompositionParams_v2(obj, target_domain = "WM", reference_domain = "Layer6")
obj <- estimateGeometryParams_v2  (obj, n_cores = 4)
obj <- estimateExpressionParams_v2(obj, target_domain = "WM", reference_domain = "Layer6",
                                   genes_to_use = sets$G_svg, n_cores = 4)

## 5 — (optional) tune spatial hyperparameters as a stable power oracle --------
tuning <- tuneSpaDesignLambdas(
  obj, target_domain = "WM", reference_domain = "Layer6",
  G_DE        = sets$G_spike,
  scenario_H1 = list(DE_lfc = c(0.30, 0.35, 0.40),
                     target_prop_case = c(0.65, 0.70, 0.75))
)

## 6a — Power for the spatial DE endpoint (SaLFC) ------------------------------
powSaLFC <- evaluatePowerSaLFC(
  obj, K_grid = 3:8, lfc_grid = c(0.3, 0.4, 0.5),
  G_svg = sets$G_svg, G_null = sets$G_null, G_spike = sets$G_spike,
  target_domain = "WM", reference_domain = "Layer6",
  n_sim = 30, tuning_obj = tuning
)
summary(powSaLFC, target_power = 0.8)   # min K per effect + FDR diagnostic
plotPowerCurveSaLFC(powSaLFC)

## 6b — Power for the composition endpoint (LOR) ------------------------------
powLOR <- evaluatePowerLOR(
  obj, K_grid = 3:8, delta_grid = c(0.60, 0.70, 0.80),
  genes = sets$G_svg, target_domain = "WM", reference_domain = "Layer6",
  effect_type = "target_prop_case", n_sim = 30, tuning_obj = tuning
)
summary(powLOR, target_power = 0.8)     # min K per effect + empirical size
plotPowerCurveLOR(powLOR)
```

> **Tip.** To prototype a single design point, call
> `simulateSpaDesign2() -> pBANSKY() -> SaLFC()/LOR()` directly (see `?SaLFC`, `?LOR`).

---

## The two endpoints

| | **SaLFC** — spatial differential expression | **LOR** — domain composition |
|---|---|---|
| **Question** | Is a gene's Target-vs-Reference log-fold-change different between case and control? | Does the target domain's abundance differ between case and control? |
| **Per-sample statistic** | z = x̄_T − x̄_R with spatial sampling variance τ², empirical-Bayes moderated (`limma::squeezeVar`) | Lₖ = log((yₖ + ε)/(bₖ + ε)) on target/reference spot counts |
| **Group comparison** | Welch *t* or Wilcoxon, Satterthwaite df | Welch *t* or Wilcoxon |
| **Null / margin** | TREAT against `lfc_threshold` (pilot plug-in τ_g) | TREAT centered on pilot baseline δ₀ = β₁, margin τ |
| **Multiplicity** | per-gene; BH (FDR) or Bonferroni | single test per (Target, Reference) pair |
| **Reported** | power & FDP over (K, DE_lfc) | power & empirical size over (K, effect) |

Both TREAT thresholds are estimated once from the pilot
(`estimate_lfc_threshold_from_pilot`, `estimate_lor_tau_from_pilot`) and shared across
replicates, so type-I error is controlled relative to a biologically meaningful
minimum effect rather than a point null.

---

## Function reference

| Stage | Functions |
|-------|-----------|
| **Object & data** | `createSpaDesign2Object`, `spaDesign2-class`, `syntheticData()`, `testingResult()`, `addTestingResult()` |
| **Feature / gene-set selection** | `featureSelection`, `featureSelectionStable`, `makeCustomGeneSets` |
| **Generative model (3 layers)** | `estimateCompositionParams_v2`, `estimateGeometryParams_v2`, `estimateExpressionParams_v2` |
| **Simulate & recover domains** | `simulateSpaDesign2`, `pBANSKY`, `rearrangeSyntheticToPilot` |
| **Endpoints & TREAT thresholds** | `SaLFC`, `LOR`, `estimate_lfc_threshold_from_pilot`, `estimate_lor_tau_from_pilot` |
| **Power evaluation & reporting** | `evaluatePowerSaLFC`, `evaluatePowerLOR`, `summary()`, `plotPowerCurveSaLFC`, `plotPowerCurveLOR`, `print()` |
| **Hyperparameter tuning** | `tuneSpaDesignLambdas` (+ `summary`/`plot` for `spaDesign2Tuning`) |
| **Visualization** | `plot(<spaDesign2>)` — pilot / synthetic spatial domains |

---

## Bundled data

| Object | Description |
|--------|-------------|
| `mini_pilot_data_list` | Downsampled human DLPFC pilot (Maynard et al. 2021, via `spatialLIBD`) in input format |
| `mini_spaDesign2_obj`  | Object created from the pilot (unfitted) |
| `mini_obj_features`    | After `featureSelection` / `featureSelectionStable` |
| `mini_obj_fitted`      | Fully fitted (composition + geometry + expression) — ready for power evaluation |
| `mini_custom_genes`    | Labeling (`G_svg`), null (`G_null`), and spike-in (`G_spike`) gene sets for "WM" vs "Layer6" |

---

## Citation

If you use `spaDesign2`, please cite the package (a methods manuscript is in
preparation):

```
Shin J, Xie J, Chung D (2026). spaDesign2: Pilot-Based Power Analysis and
Sample-Size Determination for Multi-Sample Spatial Transcriptomics.
R package version 1.0.1.
```

```bibtex
@Manual{spaDesign2,
  title  = {spaDesign2: Pilot-Based Power Analysis and Sample-Size Determination
            for Multi-Sample Spatial Transcriptomics},
  author = {Jungmin Shin and Juan Xie and Dongjun Chung},
  year   = {2026},
  note   = {R package version 1.0.1}
}
```

---

## References

- Mukhopadhyay S., Li D., Dunson D.B. (2020). Fisher–Gaussian kernels.
  *JRSS-B.* [doi:10.1111/rssb.12390](https://doi.org/10.1111/rssb.12390)
- Saha A., Datta A. (2018). BRISC: nearest-neighbor Gaussian processes.
  *Stat.* [doi:10.1002/sta4.184](https://doi.org/10.1002/sta4.184)
- Singhal V. *et al.* (2024). BANKSY. *Nature Genetics.*
  [doi:10.1038/s41588-024-01664-3](https://doi.org/10.1038/s41588-024-01664-3)
- McCarthy D.J., Smyth G.K. (2009). Testing relative to a fold-change threshold
  (TREAT). *Bioinformatics.* [doi:10.1093/bioinformatics/btp053](https://doi.org/10.1093/bioinformatics/btp053)
- Maynard K.R. *et al.* (2021). DLPFC spatial transcriptomics.
  *Nature Neuroscience* 24:425–436. (data via `spatialLIBD`, Pardo et al. 2022)

---

## Authors & license

Jungmin Shin (aut, cre · `c16267@gmail.com`) · Juan Xie (aut) · Dongjun Chung (aut)

Released under **GPL (≥ 3)**.
