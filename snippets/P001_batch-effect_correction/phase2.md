## Phase 2: Deep Grammatical Breakdown — The Syntax of Batch Correction

In Phase 1 we established the **why**. Now we learn the **how**: the exact function signatures, the arguments that matter, and the mental model behind each parameter. Think of this phase as learning the **grammar** of the language — not just vocabulary, but which words go where, and why.

> **The organizing principle:** Every batch correction tool answers three questions:
> 1. **What am I correcting?** (raw counts vs. log-normalized values)
> 2. **What must I preserve?** (the biological condition of interest)
> 3. **What is the batch variable?** (known labels vs. unknown latent factors)

The syntax of each function is a direct expression of its answer to those three questions.

---

### 1. `sva::ComBat()` — The Original Empirical Bayes Workhorse

**Package:** `sva` (Bioconductor)
**Input:** **Normalized, continuous** expression data (log-CPM, VST, TPM, microarray intensities). **Never raw counts.**

```r
ComBat(
  dat,                    # genes x samples numeric matrix (normalized, continuous)
  batch,                  # factor or vector of batch labels (only ONE batch variable)
  mod = NULL,             # model matrix for biological covariates to PRESERVE
  par.prior = TRUE,       # TRUE = parametric empirical Bayes; FALSE = non-parametric
  prior.plots = FALSE,    # TRUE = diagnostic plots of the prior
  mean.only = FALSE,      # TRUE = correct location only, not scale
  ref.batch = NULL,       # reference batch for reproducible correction
  BPPARAM = bpparam("SerialParam")  # BiocParallel backend
)
```

#### The three arguments you must understand

| Argument | What it does | Why it matters |
|----------|-------------|----------------|
| `dat` | The expression matrix, **genes in rows, samples in columns** | ComBat assumes the data are already **cleaned and normalized**. Raw counts violate the Gaussian assumption. |
| `batch` | The batch labels. **Only one batch variable is allowed per call.** | If you have two crossed batch factors (e.g., sequencing lane and operator), you must combine them into a single factor or use `limma::removeBatchEffect(batch, batch2)`. |
| `mod` | A **model matrix** encoding the biological condition(s) you want to **protect from correction**. | If you omit `mod`, ComBat may remove biological signal that happens to correlate with batch. **Always supply `mod` when you have a known condition.** |

#### The `mod` model matrix — the “shield” syntax

```r
# metadata is a data.frame with columns: batch, condition
mod <- model.matrix(~ condition, data = metadata)

# Now ComBat knows: "Remove batch, but do NOT touch the condition effect."
corrected <- ComBat(
  dat       = log_normalized_expr,
  batch     = metadata$batch,
  mod       = mod,
  par.prior = TRUE
)
```

> **Golden rule:** `mod` is your **biological shield**. It tells ComBat which variation is *wanted*. Without it, you are asking ComBat to remove “everything that looks like batch,” including biology that is confounded with batch.

#### `par.prior` — parametric vs. non-parametric

| `par.prior` | When to use | Trade-off |
|-------------|-------------|-----------|
| `TRUE` (default) | Most datasets, especially when you have **many genes** and **moderate batch sizes** | Fast; assumes a parametric prior (normal/inverse-gamma). |
| `FALSE` | When the parametric assumption is clearly violated; small number of genes | Slower; estimates the prior non-parametrically. |

#### `mean.only` — location-only correction

- `mean.only = FALSE` (default): corrects both **location (mean)** and **scale (variance)** differences between batches.
- `mean.only = TRUE`: corrects only the mean shift; leaves the variance untouched.

> Use `mean.only = TRUE` when you suspect the scale differences are **biological** rather than technical, or when the number of samples per batch is very small (scale estimation becomes unstable).

#### `ref.batch` — the “set bias” solution

```r
corrected <- ComBat(
  dat       = expr,
  batch     = batch,
  mod       = mod,
  ref.batch = "Batch1"   # all other batches are adjusted TO this reference
)
```

> When `ref.batch` is `NULL`, ComBat uses a **grand mean** across all batches. Adding new samples later and re-running ComBat can shift the corrected values of the original samples. Setting `ref.batch` anchors the correction to a fixed reference, improving **reproducibility across dataset versions**.

---

### 2. `sva::ComBat_seq()` — The Count-Aware Variant

**Package:** `sva` (Bioconductor)
**Input:** **Raw RNA-seq counts** (non-negative integers). The Negative Binomial model is used instead of the Gaussian model.

```r
ComBat_seq(
  counts,                 # genes x samples matrix of RAW integer counts
  batch,                  # batch labels
  group = NULL,           # biological condition to PRESERVE (vector or factor)
  covar_mod = NULL,       # additional covariates to preserve (model matrix)
  full_mod = TRUE,        # TRUE = include group in the model; FALSE = intercept-only
  shrink = FALSE,         # shrinkage on the batch effect estimates
  shrink.disp = FALSE,    # shrinkage on the dispersion estimates
  gene.subset.n = NULL    # number of genes for empirical Bayes (when shrink = TRUE)
)
```

#### The critical difference from `ComBat`

| Aspect | `ComBat` | `ComBat_seq` |
|--------|----------|--------------|
| Input | Log-normalized continuous values | **Raw integer counts** |
| Distribution | Gaussian | **Negative Binomial** |
| Output | Continuous corrected matrix | **Integer corrected counts** |
| Biological protection | `mod` (model matrix) | `group` (vector/factor) **or** `covar_mod` (matrix) |
| Typical use | Microarray, log-CPM, VST | **Bulk RNA-seq raw counts** |

#### The `group` argument — the count-world equivalent of `mod`

```r
# group is a simple vector, NOT a model matrix
corrected_counts <- ComBat_seq(
  counts = as.matrix(raw_counts),
  batch  = metadata$batch,
  group  = metadata$condition,   # preserve the condition effect
  full_mod = TRUE
)
```

> `group` is the **biological condition** you want ComBat-Seq to preserve. It is a vector/factor, not a model matrix. Internally, ComBat-Seq converts it into a design matrix for the Negative Binomial GLM. If you have multiple covariates, use `covar_mod` (a model matrix) instead, and set `group = NULL`.

#### `full_mod = TRUE` vs. `FALSE`

- `full_mod = TRUE` (default): the full model includes **both batch and group**. This is the correct choice when you want to preserve the biological condition.
- `full_mod = FALSE`: the null model is intercept-only. **Use only when you have no biological condition to preserve** (rare, and dangerous if misused).

> **Warning:** If `batch` and `group` are **perfectly confounded** (all cases in one batch, all controls in another), ComBat-Seq cannot distinguish technical from biological variation. It will silently remove the treatment effect. Always check your experimental design for confounding **before** correction.

---

### 3. `limma::removeBatchEffect()` — The Linear-Model Remover (Visualization Only)

**Package:** `limma` (Bioconductor)
**Input:** Log-expression values (continuous). **Not raw counts.**

```r
removeBatchEffect(
  x,                      # numeric matrix of log-expression values (genes x samples)
  batch = NULL,           # factor or vector of batch labels
  batch2 = NULL,          # optional SECOND independent series of batches
  covariates = NULL,      # numeric matrix of continuous covariates to adjust for
  design = matrix(1, ncol(x), 1)  # design matrix of effects to PRESERVE
)
```

#### The anatomy of the call

| Argument | Role | Analogy |
|----------|------|---------|
| `x` | The data to correct | The canvas |
| `batch` | Primary technical grouping | The footprint to erase |
| `batch2` | Secondary technical grouping (e.g., operator, hybridization time) | A second, independent footprint |
| `covariates` | **Continuous** technical effects (e.g., RIN score, sequencing depth) | A gradient to subtract |
| `design` | Biological effects to **keep** | The mask that protects biology from erasure |

> **`design` is the “do-not-erase” mask.** It is the limma equivalent of ComBat’s `mod`. If you omit it, `removeBatchEffect` will remove **all** variation associated with batch, including biological signal that is correlated with batch.

#### The canonical call for visualization

```r
library(limma)

# design must encode the biological condition you want to preserve
design <- model.matrix(~ condition, data = metadata)

# Correct log-normalized expression for visualization
corrected_viz <- removeBatchEffect(
  x       = log_normalized_expr,
  batch   = metadata$batch,
  design  = design
)

# Now use corrected_viz for PCA, heatmaps, MDS — NOT for DE testing
```

> **limma’s own documentation is explicit:** *“This function is not intended to be used prior to linear modelling. For linear modelling, it is better to include the batch factors in the linear model.”* The corrected matrix is for **seeing**, not for **testing**.

#### `batch` vs. `batch2` — when to use two

```r
# Two independent technical factors: sequencing lane AND operator
corrected <- removeBatchEffect(
  x      = expr,
  batch  = metadata$lane,
  batch2 = metadata$operator,
  design = design
)
```

> The two batch effects are assumed to be **additive**. If they interact (e.g., one operator produces worse data on one lane), a single combined factor is safer than two additive terms.

---

### 4. `RUVSeq` — The Control-Gene Anchor

**Package:** `RUVSeq` (Bioconductor)
**Input:** RNA-seq counts (or log counts, depending on method). **Requires control features or replicate samples** to anchor the unwanted variation estimate.

RUVSeq has **three methods**, each with a distinct syntax and a distinct source of “anchor” information.

#### 4a. `RUVg()` — Negative control genes (e.g., ERCC spike-ins, housekeeping genes)

```r
RUVg(
  x,          # genes x samples matrix (counts or log-counts)
  cIdx,       # index or names of NEGATIVE CONTROL genes (assumed constant across biology)
  k,          # number of unwanted variation factors to estimate
  drop = 0,   # number of features to drop from the top of the residual matrix
  center = TRUE,
  round = TRUE,
  epsilon = 1,
  tolerance = 1e-8,
  isLog = FALSE
)
```

> **`cIdx` is the anchor.** It tells RUVg: “These genes should not vary with biology. Any variation in them is unwanted.” If you have ERCC spike-ins, use them. If not, use a set of housekeeping genes that are reliably stable in your tissue.

#### 4b. `RUVs()` — Replicate / negative control samples

```r
RUVs(
  x,          # genes x samples matrix
  cIdx,       # index of CONTROL SAMPLES (replicates within groups)
  k,          # number of unwanted variation factors
  scIdx,      # matrix mapping replicate groups (use makeGroups())
  ...
)
```

> **`cIdx` is the anchor.** Here, the anchor is **samples**, not genes: technical replicates for which the biological condition is constant. The `makeGroups()` helper constructs the `scIdx` mapping:

```r
scIdx <- makeGroups(metadata$replicate_group)
ruv_s <- RUVs(counts, cIdx = control_sample_indices, k = 2, scIdx = scIdx)
```

#### 4c. `RUVr()` — Residuals from a first-pass GLM

```r
RUVr(
  x,          # genes x samples matrix
  cIdx,       # index of genes to use (can be all genes, or a stable subset)
  k,          # number of unwanted variation factors
  residuals,  # matrix of residuals from a first-pass GLM (e.g., edgeR)
  ...
)
```

> **`residuals` is the anchor.** RUVr does **not** require negative control genes or replicate samples. Instead, it uses the residuals from a first-pass GLM that includes the biological covariates of interest. The unwanted variation is estimated from what the biological model **fails to explain**.

#### What RUVSeq returns

All three methods return a list containing:
- `$normalizedCounts`: the corrected matrix
- `$W`: the estimated unwanted variation factors
- `$alpha`: the estimated coefficients

---

### 5. `sva::sva()` — Surrogate Variable Analysis (Unknown Latent Factors)

**Package:** `sva` (Bioconductor)
**Input:** Normalized expression data. **Discovers** latent factors when you do not know the batch labels.

```r
sva(
  dat,        # genes x samples normalized expression matrix
  mod,        # FULL model matrix (includes the biological condition)
  mod0,       # NULL model matrix (includes only the intercept, or nuisance covariates)
  n.sv = NULL, # number of surrogate variables (if NULL, estimated automatically)
  method = "irw",  # "irw" (iterative reweighted) or "two-step"
  vfilter = NULL,  # number of top-variance genes to use
  ...
)
```

#### The two model matrices — the heart of SVA

```r
# Full model: condition is the biological signal we want to model
mod  <- model.matrix(~ condition, data = metadata)

# Null model: intercept only (or nuisance covariates only)
mod0 <- model.matrix(~ 1, data = metadata)

# Discover surrogate variables
svobj <- sva(normalized_expr, mod = mod, mod0 = mod0)
```

> **`mod` and `mod0` define the contrast.** SVA finds the directions in the data that are **not explained** by the biological model (`mod`) but are **explained** by the null model (`mod0`) — i.e., the latent structure that is neither biology nor intercept. Those directions become the surrogate variables.

#### Choosing `n.sv` automatically

```r
n.sv <- num.sv(normalized_expr, mod, method = "leek")
svobj <- sva(normalized_expr, mod = mod, mod0 = mod0, n.sv = n.sv)
```

> `num.sv()` estimates the number of latent factors. If you set `n.sv = NULL`, `sva()` will call this internally. For small datasets, automatic estimation can be unstable; inspect the scree plot or the `num.sv` output.

#### How to use the surrogate variables

The surrogate variables are **not** a corrected matrix. They are **covariates** that you add to your downstream model:

```r
# Add SVs to the design matrix for DE testing
design_with_sv <- cbind(mod, svobj$sv)
fit <- lmFit(normalized_expr, design_with_sv)
```

> **SVA is a detection-and-modeling tool, not a removal tool.** It tells you *what* the latent factors are; you decide whether to model them (preferred for DE) or regress them out (for visualization).

---

### 6. Design-Matrix Inclusion — The “No Correction” Path (Preferred for DE)

**Packages:** `DESeq2`, `edgeR`, `limma`
**Input:** Raw counts. **No matrix is modified.** Batch is a term in the statistical model.

#### DESeq2 syntax

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ batch + condition   # batch FIRST, condition LAST
)
```

> **Why batch first?** DESeq2 documentation recommends putting the **variable of interest last** in the design formula. The order affects the interpretation of the coefficients and the default results table. `~ batch + condition` means: “Test condition, controlling for batch.”

#### edgeR syntax

```r
library(edgeR)

# Create a design matrix with batch and treatment
design <- model.matrix(~ Batch + Treatment, data = metadata)

# Fit the GLM
y <- DGEList(counts = counts, group = metadata$condition)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "Treatmenttreated")  # test treatment, adjust for batch
```

> In edgeR, `model.matrix(~ Batch + Treatment)` creates columns for the intercept, batch levels, and treatment levels. The `coef` argument selects which coefficient to test. Batch coefficients are estimated but not tested — they are **nuisance parameters**.

#### limma-voom syntax

```r
library(limma)

design <- model.matrix(~ batch + condition, data = metadata)
v <- voom(counts, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef = "conditiontreated")
```

> The same principle: batch is in the design, so the linear model **estimates and absorbs** the batch effect while testing the condition coefficient.

---

### 7. `pyComBat` (InMoose) — The Python Bridge

**Package:** `inmoose` (Python)
**Functions:** `pycombat_norm()` (ComBat equivalent) and `pycombat_seq()` (ComBat-Seq equivalent).

#### `pycombat_norm()` — normalized continuous data

```python
from inmoose.pycombat import pycombat_norm

corrected = pycombat_norm(
    data,        # numpy array or pandas DataFrame (genes x samples)
    batch,       # batch labels
    covar_mod=None,  # covariates to preserve (model matrix or list)
    ...
)
```

#### `pycombat_seq()` — raw RNA-seq counts

```python
from inmoose.pycombat import pycombat_seq

corrected_counts = pycombat_seq(
    counts,             # raw count matrix (genes x samples)
    batch,              # batch labels
    covar_mod=None,     # covariates to preserve (categorical only)
    shrink=False,       # shrinkage on parameter estimates
    shrink_disp=False,  # shrinkage on dispersion
    gene_subset_n=None, # genes for empirical Bayes (if shrink=True)
    ref_batch=None,     # reference batch
    na_cov_action='raise'  # 'raise', 'remove', or 'fill'
)
```

> **Critical Python-specific note:** `covar_mod` in `pycombat_seq()` must be **categorical**, not continuous. If you pass continuous covariates, the function will raise an error or misbehave. Encode them as factors before calling.

#### The Python mental model

`pyComBat` is a **direct port** of the R functions. The arguments map almost one-to-one:

| R | Python |
|---|--------|
| `ComBat(dat, batch, mod)` | `pycombat_norm(data, batch, covar_mod)` |
| `ComBat_seq(counts, batch, group)` | `pycombat_seq(counts, batch, covar_mod)` |
| `ref.batch` | `ref_batch` |
| `par.prior` | Not exposed (parametric only in pyComBat norm) |

> **Caveat:** `pycombat_norm()` does not expose the `par.prior` switch. It uses the parametric empirical Bayes path only. If you need non-parametric adjustment, you must use R’s `ComBat(par.prior = FALSE)`.

---

### 8. The Full Syntax Map — At a Glance

| Function | Input | Preserve biology via | Batch via | Output | DE-safe? |
|----------|-------|---------------------|-----------|--------|----------|
| `ComBat()` | Log-normalized continuous | `mod` (model matrix) | `batch` | Corrected continuous matrix | **No** (use for viz/ML) |
| `ComBat_seq()` | Raw integer counts | `group` (vector) or `covar_mod` | `batch` | Corrected **integer** counts | **No** (use for viz/clustering) |
| `removeBatchEffect()` | Log-expression | `design` (model matrix) | `batch` (+ `batch2`) | Corrected log-expression | **No** (documented: viz only) |
| `RUVg()` | Counts / log-counts | Control genes (`cIdx`) | Estimated factors (`k`) | `$normalizedCounts` | **No** (use factors in model) |
| `RUVs()` | Counts | Replicate samples (`cIdx`, `scIdx`) | Estimated factors (`k`) | `$normalizedCounts` | **No** |
| `RUVr()` | Counts | GLM residuals | Estimated factors (`k`) | `$normalizedCounts` | **No** |
| `sva()` | Normalized expression | `mod` / `mod0` | **Unknown** (discovered) | Surrogate variables (`$sv`) | **Yes** (add SVs to design) |
| DESeq2 `design = ~ batch + condition` | Raw counts | Condition term in formula | Batch term in formula | **No matrix modified** | **Yes** |
| edgeR `model.matrix(~ Batch + Treatment)` | Raw counts | Treatment coefficient | Batch coefficient | **No matrix modified** | **Yes** |
| `pycombat_norm()` | Normalized continuous | `covar_mod` | `batch` | Corrected matrix | **No** |
| `pycombat_seq()` | Raw counts | `covar_mod` (categorical) | `batch` | Corrected counts | **No** |

---

### 9. The Grammar in One Sentence

> **`ComBat` and `ComBat_seq` answer “make me a clean matrix.” `removeBatchEffect` answers “make me a clean picture.” `RUVSeq` answers “make me a clean matrix when I have an anchor.” `sva` answers “tell me what the hidden batch is.” Design-matrix inclusion answers “test my hypothesis without touching the data.”**

Each function’s syntax is a contract: you tell it what to protect (`mod`, `group`, `design`, `cIdx`, `covariates`) and what to remove (`batch`, `k`, latent factors). The contract is only as good as the biological shield you provide.

---

**Ready for Phase 3, the “Gotcha” section — where we compare the tools that look similar but behave differently? Any questions or wonderings so far?**
```