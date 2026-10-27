## Phase 1: Historical Context — The Genesis of Batch Effects

### The Pre-History: When “One Experiment” Was a Lie

Before high-throughput technologies, an omics experiment was small. You ran a gel, a column, or a single hybridization. You could process all your samples side-by-side, with the same reagents, on the same day, in the same room. The technical variation existed, but it was *local* and *visible*: you knew if a gel ran weirdly because you were standing right there.

Then came microarrays in the late 1990s and early 2000s, and later RNA-seq in the late 2000s. Suddenly, you could measure thousands of genes at once, but you could only process a handful of samples per array or per sequencing lane. A study of 200 samples could not be run on one chip, in one day, by one person. It had to be split into **batches**.

> **The original problem:** High-throughput measurement made it impossible to process all samples under identical conditions. Technical variation from reagent lots, array scanners, sequencing runs, personnel, and even the day of the week became *systematic* and *confounded* with the biological question.

### The Birth of a Named Enemy: 2001–2007

The term **“batch effect”** entered the omics vocabulary as researchers noticed a disturbing pattern: samples processed together clustered together in PCA plots, regardless of their biological phenotype. A tumor sample run on Monday looked more like a normal sample run on Monday than like a tumor sample run on Tuesday.

Early efforts relied on **normalization** (e.g., quantile normalization for microarrays). But normalization assumes the *global distribution* of expression is the same across samples — an assumption that fails when batches differ in their technical bias.

A critical conceptual shift happened around 2005–2007: batch effects were no longer seen as mere “noise to be normalized away.” They were recognized as **systematic, estimable, and removable** if modeled explicitly.

The landmark moment arrived in **2007** with the publication of **ComBat** by Johnson, Li, and Rabinovic: *“Adjusting batch effects in microarray expression data using empirical Bayes methods.”* [reference:0] ComBat introduced a crucial idea:

> **Empirical Bayes shrinkage.** Instead of estimating batch parameters for each gene independently (which would be noisy with small batches), ComBat borrows information across genes to stabilize the estimates.

This is the **“Janitor vs. Detective”** distinction in its earliest form:

| Role                             | Analogy                                       | What it does                                                                                     |
| -------------------------------- | --------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| **Normalization (Janitor)**      | Wipes every surface equally                   | Adjusts global distribution; assumes all samples should look alike at the bulk level.            |
| **Batch correction (Detective)** | Finds the *specific* footprints of each batch | Models batch-specific shifts and scales, then subtracts them while preserving biological signal. |

ComBat became the de facto standard because it worked, it was implemented in the **`sva` R package**, and it handled both **known batches** (explicit covariates) and, through related methods like **SVA (Surrogate Variable Analysis)**, **unknown latent factors**. 

### The RNA-seq Detour: Why ComBat Wasn’t Enough (2010s)

ComBat was designed for **continuous, approximately Gaussian microarray intensities**. RNA-seq data are **counts** — non-negative integers with a mean–variance relationship (Negative Binomial). Applying ComBat directly to raw counts is statistically inappropriate: it can produce negative “expression” values and distort the count distribution.

For years, the practical workaround was to apply ComBat to **log-transformed, normalized counts** (e.g., `log2(CPM + 1)` or VST-transformed data). This worked reasonably well for many bulk RNA-seq studies, but it was a square peg in a round hole.

The proper solution arrived in **2020** with **ComBat-Seq** (Zhang, Parmigiani, and Johnson), which models raw counts with a **Negative Binomial** distribution and returns **adjusted counts** that remain integers.  ComBat-Seq is now implemented in the same **`sva` R package**, alongside the original ComBat.

### The Python Migration: 2020s

For two decades, batch correction was overwhelmingly an **R** story: `sva`, `limma` (`removeBatchEffect`), `RUV`, and later `Seurat` for single-cell integration. Python users had to call R via `rpy2` or accept suboptimal workarounds.

The situation changed with **pyComBat** (2023), a Python implementation of ComBat and ComBat-Seq based on the *same mathematical framework*, distributed as part of the **`InMoose`** package. [reference:3] It provides similar correction performance at comparable or better speed, and it finally gives Python-native bioinformatics pipelines access to the gold-standard empirical Bayes approach. [reference:4]

### The Unresolved Tension: Over-Correction and the “Set Bias” Problem

From the very beginning, batch correction carried a philosophical risk:

> **If you remove too much technical variation, you may also remove biological signal.**

This is the **“over-correction”** problem. If batch is perfectly confounded with phenotype (e.g., all tumors run in batch 1, all normals in batch 2), no statistical method can separate them. ComBat will “correct” the difference — and in doing so, erase the biology.

A related practical issue is **“set bias”**: the original ComBat uses all samples to estimate batch parameters. If you add new samples later and re-run ComBat on the combined dataset, the corrected values for the *original* samples can change. [reference:5] This makes biomarker discovery non-reproducible across dataset versions. Proposed solutions include **M-ComBat** (adjusting to a predefined reference batch) and **reference-batch ComBat** workflows. [reference:6]

### Timeline Summary

| Year       | Milestone                                                  | Why it mattered                                                                                              |
| ---------- | ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| ~2001–2005 | Batch effects formally recognized in microarray studies    | Technical variation shown to confound biological interpretation                                              |
| **2007**   | **ComBat** (Johnson et al.)                                | Empirical Bayes framework; became the standard for microarray batch adjustment                               |
| 2007–2008  | **SVA** (Leek & Storey)                                    | Detection and adjustment of *unknown* latent batch effects                                                   |
| 2010s      | RNA-seq adoption; ComBat applied to log-transformed counts | Practical but statistically imperfect                                                                        |
| **2020**   | **ComBat-Seq** (Zhang et al.)                              | Negative Binomial model for raw RNA-seq counts; returns adjusted counts                                      |
| **2023**   | **pyComBat / InMoose**                                     | Python-native ComBat and ComBat-Seq; bridges the R–Python gap                                                |
| 2020s      | Single-cell and multi-omics batch correction               | New methods (MNN, Harmony, Seurat integration) for cell-level and cross-platform heterogeneity [reference:7] |

### The Core Historical Lesson

Batch effects were not a “bug” in early omics analysis. They were an **inevitable consequence of scaling up measurement**. The field’s response — from normalization to explicit batch modeling to distribution-aware count correction to cross-language implementations — reflects a deeper truth:

> **Technical variation is not random noise. It is structured. And structured variation can be modeled, estimated, and removed — but only if you respect the data’s generative distribution and the study design.**

---

**Ready for Phase 2? Any questions or wonderings so far?**