# Differential Expression in R: edgeR, DESeq2, and limma-voom

## 1. The big picture

Differential-expression analysis asks whether a gene's expression differs systematically between experimental conditions, relative to the uncertainty in its estimate, while accounting for the experimental design and testing many genes at once.

**The experimental design defines the question; the model defines how the data and uncertainty are represented; the contrast defines the comparison being tested.** The package is only one part of the analysis.

| Method | Historical origin / typical use | Core model | What it does with RNA-seq counts |
|---|---|---|---|
| **edgeR** | Developed for digital gene-expression count data | Negative-binomial generalized linear models (GLMs); commonly a quasi-likelihood (QL) workflow for replicated studies | Models counts directly; estimates normalization factors and dispersion; tests coefficients or contrasts |
| **DESeq2** | Successor to DESeq for count-based differential expression | Negative-binomial GLMs | Models counts directly; estimates size factors and dispersions; supports Wald and likelihood-ratio tests, plus optional log-fold-change shrinkage |
| **limma-voom** | Extends limma's linear-model framework to RNA-seq | Weighted linear models on log-counts-per-million (logCPM) | `voom()` transforms counts to logCPM and estimates observation-level precision weights; limma fits weighted linear models and moderates residual variances |

The microarray-era limma framework was designed around expression intensities and linear models. RNA-seq introduced discrete counts and a different mean–variance relationship, inspiring edgeR and DESeq. **limma-voom bridges these approaches**: it uses counts to estimate weights, then uses weighted linear models.

## 2. Why the models differ

### Negative-binomial models: edgeR and DESeq2

For gene *i* in sample *j*, the count model is commonly expressed as:

\[
K_{ij} \sim \operatorname{NB}(\mu_{ij}, \alpha_i), \qquad
\operatorname{Var}(K_{ij}) = \mu_{ij} + \alpha_i\mu_{ij}^{2}.
\]

- `K_ij`: observed count.
- `mu_ij`: expected count, related to library scaling and experimental covariates.
- `alpha_i`: gene-specific dispersion, representing variability beyond basic count-sampling variation.

A log link relates the expected expression to predictors such as condition and batch. The models do not simply log-transform raw counts and then treat them as normally distributed.

### Linear models: limma-voom

Voom calculates logCPM values and estimates precision weights to account for the abundance-dependent variance of RNA-seq measurements. limma then fits weighted linear models and uses empirical Bayes moderation of gene-wise residual variances.

Logging is useful because expression changes are often multiplicative and fold changes are naturally expressed on the log scale. This is **not simply because RNA-seq involves PCR**; the transformation and voom weights are used to represent the statistical behavior of the measured counts.

## 3. Empirical Bayes and borrowing information

With few biological replicates, each gene contains limited information for estimating its variability. These methods use information across genes to stabilize gene-level estimates. This is **borrowing information across genes**, not selecting only stable samples.

- **edgeR:** empirical Bayes methods help estimate dispersion; the QL workflow also accounts for uncertainty in gene-wise variability.
- **DESeq2:** estimates dispersions using gene-wise information and a fitted mean–dispersion trend; optional `lfcShrink()` methods stabilize noisy log-fold-change estimates.
- **limma:** `eBayes()` moderates gene-wise residual variances to improve inference from t- and F-statistics.

These procedures can reduce the impact of noisy parameter estimates, but they do not replace biological replication, repair confounding, or guarantee that a model is correct. Dispersion/variance moderation and log-fold-change shrinkage are related but distinct operations.

## 4. Input data and normalization

### Standard inputs

- Use a **gene-by-sample raw count matrix** for the standard edgeR and DESeq2 count models.
- Use raw counts as the starting point for the standard limma-voom workflow.
- TPM, FPKM/RPKM, and already log-transformed values are **not substitutes for raw counts** in the standard negative-binomial workflows. Voom is not a generic function for arbitrary transformed expression matrices.
- If the quantifier is Salmon or another transcript-level tool, use a supported import workflow (for example, `tximport` where appropriate) rather than relabeling TPM as counts.

### Normalization is not batch correction

| Method | Typical normalization | Basic idea |
|---|---|---|
| edgeR | TMM via `calcNormFactors()` | Estimates sample scaling factors while reducing composition bias |
| DESeq2 | Median-ratio size factors | Uses ratios to gene-wise reference values to estimate sample scaling |
| limma-voom | Often TMM followed by `voom()` | TMM estimates scaling; voom estimates logCPM and observation-level precision weights |

These methods address relative library size/composition differences. **They do not automatically remove batch effects or other confounding.** Ordinary sequencing counts are compositional, so absolute global RNA changes may require suitable external references or spike-ins and an analysis designed for that target.

## 5. Design matrices, factors, and contrasts

Example metadata columns: `condition` (`Control`, `Cancer`) and `batch` (`B1`, `B2`).

```r
design <- model.matrix(~ batch + condition, data = samples)
colnames(design)
```

With standard treatment coding and `Control` as the reference level, `conditionCancer` usually represents the cancer-versus-control effect adjusted for batch.

Useful R formula syntax:

| Formula | Meaning |
|---|---|
| `~ condition` | Condition effect relative to a reference level |
| `~ batch + condition` | Additive effects; condition effect is adjusted for batch |
| `~ patient_id + condition` | Paired/repeated measurements, with patient-specific baseline effects |
| `~ genotype * condition` | Main effects plus genotype-by-condition interaction |
| `~ genotype + condition + genotype:condition` | Equivalent expansion of the previous formula |
| `~ condition - 1` | No intercept; estimates a coefficient for each condition level under standard coding |

Check alignment before modeling:

```r
stopifnot(!anyDuplicated(samples$sample_id))
stopifnot(all(samples$sample_id %in% colnames(counts)))
counts <- counts[, samples$sample_id, drop = FALSE]
rownames(samples) <- samples$sample_id
stopifnot(identical(colnames(counts), rownames(samples)))
```

Also inspect whether the model matrix is estimable:

```r
design <- model.matrix(~ batch + condition, data = samples)
stopifnot(qr(design)$rank == ncol(design))
```

A full-rank check is useful but does not replace design inspection or guarantee sufficient replication.

### Confounding versus interaction

These are different concepts:

- **Confounding / non-identifiability:** two effects cannot be separated from the observed design. Example: every control is in batch B1 and every cancer sample is in B2. Adding `~ batch + condition` cannot recover the two separate effects.
- **Interaction / effect modification:** the effect of one factor differs across levels of another factor.
- **Statistical independence:** a separate probability concept; it is not the definition of interaction.

For a WT/KO by control/treatment experiment:

```r
~ genotype + condition
```

assumes the treatment log-fold change is the same in both genotypes. In contrast:

```r
~ genotype * condition
```

allows treatment response to differ by genotype. With WT and Control as reference levels, the interaction coefficient represents the difference between the treatment effect in KO and the treatment effect in WT, on the model's log-expression scale.

Example:

| Genotype | Control | Treatment | Treatment fold change |
|---|---:|---:|---:|
| WT | 100 | 200 | 2x |
| KO | 100 | 400 | 4x |

The log2 treatment effects are 1 for WT and 2 for KO; the interaction is `2 - 1 = 1`. By contrast, WT 100→200 and KO 300→600 both have a 2x change, so the interaction is zero even though KO has a higher baseline.

Do not infer a group difference merely because an effect is significant in one group and not significant in another. Test the interaction directly.

### Paired designs

For before/after samples from the same participants, a common design is:

```r
~ patient_id + condition
```

This accounts for stable patient-specific baseline differences while estimating a common within-patient condition effect. A patient-by-condition interaction is not automatically necessary; if there is only one before and one after observation per patient, it can overparameterize or saturate the model. Interaction modeling requires a scientific rationale and adequate information to estimate effect heterogeneity.

## 6. Typical analysis workflows

The examples assume `counts` is a numeric gene-by-sample raw-count matrix; `samples` contains the matching sample metadata; `condition` has levels `Control` and `Cancer`; and `batch` exists. If there is no batch variable, use `~ condition` and adapt coefficient names accordingly.

### edgeR: quasi-likelihood workflow

```r
library(edgeR)

design <- model.matrix(~ batch + condition, data = samples)
y <- DGEList(counts = counts)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "conditionCancer")
res_edgeR <- topTags(test, n = Inf)$table
```

Typical sequence: filter → normalize → estimate dispersion → fit QL model → test coefficient/contrast.

### DESeq2: integrated negative-binomial workflow

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ batch + condition
)
dds <- DESeq(dds)
res_deseq2 <- results(
    dds,
    contrast = c("condition", "Cancer", "Control")
)
res_deseq2 <- as.data.frame(res_deseq2)
```

`DESeq()` carries out the core size-factor estimation, dispersion estimation, model fitting, and testing workflow. `results()` defines the comparison explicitly.

### limma-voom: weighted linear models

```r
library(edgeR)
library(limma)

design <- model.matrix(~ batch + condition, data = samples)
y <- DGEList(counts = counts)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)
res_limma <- topTable(
    fit,
    coef = "conditionCancer",
    number = Inf,
    sort.by = "P"
)
```

`voom()` estimates logCPM and precision weights; `lmFit()` fits weighted models; `eBayes()` moderates residual variances; `topTable()` extracts results.

**When comparing packages:** align the samples, gene definitions, design, contrast direction, and preferably filtering strategy. Results may still differ because the statistical procedures differ. A longer significant-gene list does not establish that a method is better.

## 7. Interpreting results

Common result columns:

| Meaning | edgeR | DESeq2 | limma |
|---|---|---|---|
| Estimated log2 fold change | `logFC` | `log2FoldChange` | `logFC` |
| Raw p-value | `PValue` | `pvalue` | `P.Value` |
| Adjusted p-value | `FDR` | `padj` | `adj.P.Val` |

For Cancer relative to Control:

- `log2FC = +1`: estimated expression is 2x higher in Cancer.
- `log2FC = -1`: estimated expression is half as high in Cancer.
- `log2FC = 0`: no estimated difference.

Adjusted p-values account for multiple testing across genes. Effect size, statistical significance, and predictive performance are distinct properties. For example, DESeq2's optional `lfcShrink()` stabilizes effect-size estimates; it should not be confused with dispersion estimation or automatically substituted into every testing step.

## 8. Basic QC and validation checklist

1. **Sample identity and metadata:** verify IDs, sample order, group labels, duplicates, missing values, and gene identifiers.
2. **Library and expression QC:** inspect library sizes, detected genes, sample correlations, PCA/MDS, and potential outliers.
3. **Design:** confirm sufficient biological replication, model estimability, and reasonable representation of groups across technical batches.
4. **Filtering and normalization:** document the low-expression filtering method and normalization assumptions.
5. **Results:** review effect sizes, adjusted p-values, variance/dispersion diagnostics, and possible outlier influence.
6. **Sensitivity and biological validation:** investigate consequential disagreements between methods; validate key findings independently when possible. Packages applied to the same data are not independent replications.

For small-RNA sequencing, make sure each row represents a consistently defined feature (for example, a mature miRNA, a specific tRNA-derived fragment, or a sequence-defined feature). Annotation, multi-mapping, and read-assignment policies must be consistent across samples; a statistical model cannot repair inconsistent feature definitions.

## 9. Bash automation and the safety rule

When wrapping the analysis in Bash, always double-quote variable expansions, for example `"${sample_id}"` and `"${fastq}"`. This prevents common word-splitting and pathname-expansion bugs.

```bash
#!/usr/bin/env bash
set -euo pipefail

sample_file="sample_ids.txt"

if [[ ! -r "${sample_file}" ]]; then
    printf 'ERROR: cannot read %s\n' "${sample_file}" >&2
    exit 1
fi

while IFS= read -r sample_id || [[ -n "${sample_id:-}" ]]; do
    [[ -n "${sample_id}" ]] || continue
    fastq="data/raw/${sample_id}_R1.fastq.gz"

    if [[ ! -s "${fastq}" ]]; then
        printf 'ERROR: missing or empty file: %s\n' "${fastq}" >&2
        exit 1
    fi

    printf 'OK: %s\n' "${fastq}"
done < "${sample_file}"
```

- `[[ -s "${fastq}" ]]` is true when the file exists and has nonzero size.
- `while IFS= read -r` preserves whitespace and backslashes in input lines.
- `|| [[ -n "${sample_id:-}" ]]` also processes the last line when the file lacks a trailing newline.
- **Glob and regex are different:** `*.fastq.gz` is a shell glob for filename matching; `^sample_.*\.fastq\.gz$` is a regular expression for matching text. Neither is the same as an R model formula such as `~ batch + condition`.

A reproducible project should record input checksums, software/package versions, reference genome and annotation versions, parameters, logs, and output locations.

## 10. Practical rules to remember

1. **Design first:** a statistical package cannot fix a confounded experiment.
2. **Raw counts for count models:** do not substitute TPM/FPKM or arbitrary transformed values into standard negative-binomial workflows.
3. **Normalization is not batch correction:** these solve different problems.
4. **Main effects are not interactions:** interaction tests whether an effect changes across levels of another factor.
5. **Contrasts define direction:** explicitly specify Cancer versus Control and inspect coefficient names.
6. **Empirical Bayes stabilizes estimates; it does not replace replication.**
7. **Do not choose the method because it discovers the most genes.** Evaluate assumptions, diagnostics, robustness, and biological validity.
8. **Double-quote Bash variables:** `"${var}"` is the default safe habit.

## Further reading

- edgeR user guide and workflow: <https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf>
- DESeq2 vignette: <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>
- limma user guide: <https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>
- Law et al., *voom: Precision weights unlock linear model analysis tools for RNA-seq read counts*, Genome Biology (2014): <https://doi.org/10.1186/gb-2014-15-2-r29>
- Love et al., *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*, Genome Biology (2014): <https://doi.org/10.1186/s13059-014-0550-8>
- Robinson et al., *edgeR: a Bioconductor package for differential expression analysis of digital gene expression data*, Bioinformatics (2010): <https://doi.org/10.1093/bioinformatics/btp616>
