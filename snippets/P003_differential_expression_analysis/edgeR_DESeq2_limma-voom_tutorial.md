
## Grammar Parts
### 0. the basic setup

```R
# Rows = genes; columns = samples
counts <- matrix(
    c(120, 135, 130, 320, 350, 340,
       12,  27,  20,  15,  29,  22,
     1200, 1450, 1300, 1300, 1410, 1350),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
        c("GeneA", "GeneB", "GeneC"),
        c("Control1", "Control2", "Control3",
          "Cancer1", "Cancer2", "Cancer3")
    )
)

# Rows = samples; columns = sample characteristics
coldata <- data.frame(
    condition = factor(
        c("Control", "Control", "Control",
          "Cancer", "Cancer", "Cancer"),
        levels = c("Control", "Cancer")
    ),
    row.names = colnames(counts)
)

coldata$batch <- factor(
    c("B1", "B1", "B2", "B1", "B2", "B2")
)

# counts: what was measured
# coldata: under which condition each measurement was obtained

stopifnot(
    identical(colnames(counts), rownames(coldata))
)


```

### 1. differential analysis
- define the design matrix

A crucial limitation is that the model cannot separate condition from batch if they are perfectly confounded. For example, if every control was sequenced in batch B1 and every cancer sample in B2, the data cannot independently identify those two effects.

```R
design <- model.matrix(
    ~ batch + condition,
    data = coldata
)

# inspect the design matrix
colnames(design)
design
```

- 1.1 edgeR: Function-by-function syntax

```R
library(edgeR)

# 1. Create the count-data object
y <- DGEList(counts = counts)

# 2. Remove genes with insufficient expression
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]

# 3. Estimate normalization factors
y <- calcNormFactors(y, method = "TMM")

# 4. Estimate gene-wise dispersions
y <- estimateDisp(y, design, robust = TRUE)

# 5. Fit the quasi-likelihood model
fit <- glmQLFit(y, design, robust = TRUE)

# 6. Test the cancer effect relative to control
test <- glmQLFTest(
    fit,
    coef = "conditionCancer"
)

# 7. Extract ranked results
res_edgeR <- topTags(test, n = Inf)$table
```

Understand the verbs rather than memorizing the code.
DGEList() stores counts and associated library information.
filterByExpr() identifies genes with sufficient expression for the experimental design. Filtering low-information genes reduces the burden of testing uninformative hypotheses.
calcNormFactors() estimates scaling factors, with TMM as the default normalization method. It does not transform your count matrix into the final expression values.
estimateDisp() estimates the variability beyond the basic count-sampling process.
glmQLFit() fits the quasi-likelihood model.
glmQLFTest() tests the coefficient of interest using a quasi-likelihood F-test.
topTags() extracts and ranks the results, including adjusted significance values.
For real experiments, the design matrix must correspond exactly to the sample order in the count object.


- 1.2 DESeq2: A different API for a related model
```R
library(DESeq2)

# 1. Build the experiment object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ batch + condition
)

# 2. Example of a simple low-count filter
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 3. Run normalization, dispersion estimation,
#    model fitting, and statistical testing
dds <- DESeq(dds)

# 4. Compare cancer against control
res_deseq2 <- results(
    dds,
    contrast = c("condition", "Cancer", "Control")
)

# 5. Inspect the most significant genes
res_deseq2 <- res_deseq2[
    order(res_deseq2$padj),
]
```

In DESeq2, contrast = c("condition", "Cancer", "Control") explicitly requests cancer relative to control. The sign of the reported log2 fold change follows that direction.
DESeq2 also separates significance testing from optional fold-change shrinkage. For example, lfcShrink() can provide stabilized effect-size estimates for visualization and ranking. Those estimates should not be confused with the unshrunken coefficients used in every testing workflow.

- 1.3 limma-voom: Counts first, weighted linear models afterward

```R
library(edgeR)
library(limma)

# 1. Store counts
y <- DGEList(counts = counts)

# 2. Filter low-expression genes
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]

# 3. Normalize library composition
y <- calcNormFactors(y)

# 4. Transform counts and estimate precision weights
v <- voom(y, design, plot = TRUE)

# 5. Fit weighted linear models
fit <- lmFit(v, design)

# 6. Moderate gene-wise variances
fit <- eBayes(fit, robust = TRUE)

# 7. Extract the condition effect
res_limma <- topTable(
    fit,
    coef = "conditionCancer",
    number = Inf,
    sort.by = "P"
)
```

Notice the sequence:
voom() → lmFit() → eBayes() → topTable()
The roles are distinct:
- voom() estimates precision weights for logCPM values.
- lmFit() fits the gene-wise linear models.
- eBayes() moderates the residual variance estimates across genes.
- topTable() extracts the ranked differential-expression results.
The fact that this workflow uses linear models does not mean that it ignores the count-generating process. Voom is the component that connects the RNA-seq counts to the weighted linear-model framework.


- 1.4 Comparison of results: Equivalent interpretations

The three methods report different column names, but many quantities have equivalent interpretations.

| Meaning | edgeR | DESeq2 | limma |
|:---:|:---:|:---:|:---:|
| Estimated expression effect | logFC | log2FoldChange | logFC |
| Raw significance | PValue | pvalue | P.Value |
| Adjusted significance | FDR | padj | adj.P.Val |

filter significant genes
```R
# An example reporting filter, not a complete analysis
sig_edgeR <- subset(
    res_edgeR,
    FDR < 0.05 & abs(logFC) > 1
)
```


## Gotcha Parts
### 1. The master comparison table
| Feature | edgeR | DESeq2 | limma-voom |
|:---:|:---:|:---:|:---:|
| Main model | Negative-binomial GLM |	Negative-binomial GLM |	Weighted linear model |
| Primary input | Raw counts | Raw counts | Raw counts | Raw counts, transformed internally |
| Typical normalization | TMM | Median-ratio size factors |	Often TMM, followed by voom |
| Variability estimation |NB dispersion; QL dispersion in QL workflow |	Gene-wise dispersion estimates moderated toward a fitted trend | Mean–variance weights plus moderated residual variances |
| Common testing approach |	Quasi-likelihood | F-test	Wald test or likelihood-ratio test | Moderated t-test or F-test |
| Fold-change shrinkage	| Available through related methods; not identical to DESeq2's workflow	| Explicit lfcShrink() options |	Not part of the ordinary eBayes() procedure in the same way |
| Experimental designs |	Flexible GLMs | Flexible GLMs | Flexible linear models |
| Main conceptual strength |	Count modeling and quasi-likelihood inference |	Integrated count-modeling workflow | Flexible design modeling and weighted linear-model machinery |

### 2. The NB modeling does not mean identical statistical testing
Do not assume that two packages must generate identical p-values because they use negative-binomial models.

The edgeR-specific nuance
edgeR's quasi-likelihood approach models **additional gene-specific variability** beyond the negative-binomial mean–variance relationship. It uses this information in the quasi-likelihood test, providing an additional layer of uncertainty assessment. 

The DESeq2-specific nuance
DESeq2 typically uses a Wald test for a specified coefficient or contrast. A likelihood-ratio test can instead compare a full model with a reduced model when the hypothesis concerns multiple terms.

### 3. Empirical Bayes does not mean the same thing in all three packages
Three kinds of statistical borrowing

edgeR: dispersion estimation
Genes contribute to the estimation of dispersion distributions and trends. In the quasi-likelihood workflow, additional gene-specific variability is considered.


DESeq2: dispersion and fold-change shrinkage
Dispersion estimates are moderated toward a fitted mean–dispersion relationship. Optional fold-change shrinkage stabilizes noisy effect estimates, particularly when the evidence is weak.


limma: residual-variance moderation
Gene-wise residual variances are moderated using information across genes, improving the stability of moderated t- and F-statistics.

Notice that none of these procedures selects a collection of stable samples and restricts the comparison to those samples. Instead, the model uses all suitable biological replicates while borrowing information across genes to estimate statistical parameters.
This is also why empirical Bayes should not be described simply as an anti-overfitting mechanism. It stabilizes estimates and improves inference under the model assumptions, but it cannot create biological replication or correct a confounded experiment.

### 4. Raw counts, TPM, CPM, and logCPM are not interchangeable inputs

| Data representation | edgeR NB workflow | DESeq2 NB workflow | limma-voom |
|:---:|:---:|:---:|:---:|
| Raw gene-level integer counts	| Yes | Yes | Yes, with a count-aware preprocessing workflow |
|TPM|Not as input to the standard NB count model|	Not as input to the standard NB count model |	Not a substitute for raw counts in the standard voom workflow|
|FPKM/RPKM | 	Not as input to the standard NB count model	| Not as input to the standard NB count model	| Not a substitute for raw counts in the standard voom workflow|
|LogCPM |	Not a replacement for counts in the standard NB workflow	| Not a replacement for counts in the standard NB workflow	|Yes, voom generates logCPM and precision weights internally|

The negative-binomial likelihood is defined for count observations. TPM and FPKM/RPKM are abundance-normalized quantities, not raw count observations. Feeding them into an ordinary negative-binomial count model violates the intended input model.
Likewise, voom is not merely a function that accepts any logarithmic expression matrix. Its precision-weight estimation is designed around the relationship between sequencing counts, library sizes, and logCPM.
For example, if your upstream quantifier is Salmon, a supported workflow can import transcript-level estimates and produce appropriate gene-level count estimates using tools such as tximport. You should not simply rename a TPM column to counts and pass it to DESeq2.

### 5. More significant genes does not necessarily mean a better analysis
The difference could arise from effect estimates, variance or dispersion estimation, test-statistic calibration, filtering, outlier handling, or the data-generating conditions.
Benchmark studies have found that method performance and agreement vary with sample size, dispersion, the proportion of truly differentially expressed genes, outliers, and other experimental characteristics. There is no universal ranking that applies to every experiment. 

The correct approach is not to choose whichever package produces the longest gene list.
Instead, assess whether the workflow is appropriate for the design, whether assumptions are plausible, whether effect sizes are stable, whether diagnostics raise concerns, and whether biological conclusions are robust.
For a high-stakes study, reasonable sensitivity analyses may include comparing the three methods on the same properly filtered count matrix and design, then investigating meaningful disagreements instead of taking a majority vote.
Do not treat the three gene lists as independent replications: they analyze the same data and share many assumptions.

- A practical decision guide

Choose edgeR when count-modeling flexibility and quasi-likelihood inference are central.
It is a well-established option for bulk RNA-seq with replicated groups, multi-factor designs, and contrasts. Start with its quasi-likelihood GLM workflow for a conventional replicated study.

Choose DESeq2 when you want an integrated count-modeling workflow and convenient effect-size shrinkage.
Its structured object model, normalization and dispersion estimation, testing, independent filtering, outlier diagnostics, and optional LFC shrinkage make it a practical general-purpose choice.

Choose limma-voom when the flexibility of the limma linear-model framework is especially useful.
It supports complex experimental designs and a broad statistical ecosystem through weighted linear models. It is particularly attractive for multi-factor experiments and larger studies.