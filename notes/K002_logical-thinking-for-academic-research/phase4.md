## Phase 3 — Point-by-Point Divergence: The "Gotcha" Section

Logical thinking is a **family resemblance** concept. Several neighbors look almost identical but fail differently. The table below is the one to keep.

### Table 1 — The central trio: logical thinking vs. statistical significance vs. rhetorical persuasion

| Dimension | Logical thinking (research) | Statistical significance | Rhetorical persuasion |
|---|---|---|---|
| **Core question** | "Is the inference warranted given the design?" | "Is the observed effect unlikely under \(H_0\)?" | "Is the claim convincing to an audience?" |
| **What it operates on** | The whole chain: boundary → claim → design → selection → interpretation | One number: \(p\), or a corrected \(p\) | Language, framing, narrative, authority |
| **Failure mode** | Unfalsifiable or confounded claim | Significant but trivial, or significant by multiplicity | Persuasive but unsupported |
| **Can be satisfied while the others fail?** | Yes — a logically sound study can be underpowered or unpersuasive | Yes — \(p < 0.05\) in a perfectly confounded design | Yes — a beautifully written paper with a broken design |
| **What it protects against** | Overreach and self-deception | Random fluctuation under a null | Being ignored |
| **What it does NOT protect against** | Bad measurement, wrong model | Confounding, selection, interpretation | Truth |

> **The most important row is the last one.** Statistical significance is a *narrow* property of a *single* comparison under a *specified* null. It says nothing about whether the comparison was fair, whether the null was the right one, or whether the claim matches the design. A study can be statistically significant and logically invalid at the same time. This is not a paradox; it is the normal state of a confounded experiment.

### Table 2 — Hypothesis-driven vs. data-driven vs. hybrid

| Dimension | Hypothesis-driven | Data-driven | Hybrid (discovery → validation) |
|---|---|---|---|
| **Direction** | Theory → data | Data → theory | Data → hypothesis → new data |
| **Strength** | Sharp discrimination; known error rate; interpretable | Finds the unexpected; no theory blinders | Both, if the split is honest |
| **Vulnerability** | Confirmation bias; misses unanticipated structure; "theory-laden" measurement | Overfitting; multiplicity; post-hoc storytelling | Leakage across the split; "validation" on a non-independent set |
| **Provenance** | Claim fixed before data | Claim emerges from data | Claim fixed before *validation* data |
| **Correct label** | Confirmatory | Exploratory (hypothesis-generating) | Confirmatory, *for the validation phase only* |
| **Correct language** | "We tested whether…" | "We observed…; this suggests…" | "We discovered… and confirmed in an independent cohort" |
| **Fatal misuse** | Forcing data into a preferred model | Reporting exploration as confirmation | Calling a re-split of the same data "independent" |

### Table 3 — Confirmatory vs. exploratory vs. descriptive vs. predictive

These four are constantly confused. They are **different claim types** with different validity criteria.

| Claim type | Question | Validity criterion | Omics example |
|---|---|---|---|
| **Descriptive** | What is in the data? | Accuracy of measurement and summary | "The cohort contains 412 differentially methylated regions." |
| **Exploratory** | What patterns might be real? | Plausibility + explicit labeling as hypothesis-generating | "Clustering suggests an immune subgroup." |
| **Confirmatory** | Is the pre-specified claim supported? | Pre-specification + error control + discrimination | "The pre-specified 12-gene signature separates cases from controls, FDR < 0.05." |
| **Predictive** | Does it generalize to new samples? | Out-of-sample performance (calibration, AUC, Brier) | "The classifier achieves 0.81 AUC on a held-out cohort." |

> **A finding can be descriptive and true, exploratory and true, and still not confirmatory.** The claim type determines what evidence is needed. Mixing them is the single most productive source of inflated papers.

### Table 4 — Neighbor concepts that look like logical thinking but aren't

| Neighbor | What it actually does | Where it diverges from logical thinking |
|---|---|---|
| **Critical thinking (general)** | Evaluates arguments for internal consistency | Can be perfectly consistent and empirically unfalsifiable |
| **Skepticism** | Doubts claims | Doubt without a discriminating test is just contrarianism |
| **Rigor** | Procedural care and precision | Rigor in service of a trivial or confounded question is wasted |
| **Reproducibility** | Same data + same analysis → same result | A reproducible confounded analysis is still confounded |
| **Replicability** | New data + same design → same conclusion | Stronger than reproducibility, still not causation |
| **Statistical power** | Probability of detecting an effect if present | Power without a valid contrast detects the confounder efficiently |
| **Pre-registration** | Fixes the analysis plan in advance | Can be gamed; a plan written after pilot data is theater |
| **Peer review** | Social gatekeeping | Catches some errors; cannot audit data it never sees |

### The single most common beginner confusion

> **Believing that "data-driven" means "hypothesis-free," and therefore "unbiased."**

This is the deepest and most damaging confusion, and it deserves its own paragraph.

The reasoning feels like this: *"I didn't impose a hypothesis, so I didn't impose a bias. The data spoke for themselves."*

The mechanics of why this is wrong:

1. **A hypothesis was still present — it was just implicit.** You chose the normalization, the covariates, the clustering method, the distance metric, the number of clusters, the significance threshold, the gene-set database, the pathway to highlight. Every one of those choices encodes an assumption about what matters. The hypothesis did not vanish; it migrated into the pipeline.
2. **The absence of a pre-specified hypothesis removes the *constraint*, not the *bias*.** It increases your freedom to select a flattering result. Freedom of analysis is freedom to overfit.
3. **The data cannot "speak."** They are consistent with many stories. The selection of which story to tell is a human act, and it is exactly where the logical work must happen.

The correct discipline for data-driven work is not "have no hypothesis." It is:

- **Declare the mode**: "This is exploratory."
- **Freeze the pipeline** as far as possible before looking at the outcome.
- **Report the search space**: how many tests, how many cluster resolutions, how many covariates tried.
- **Validate on data that did not participate in the discovery.**

> **The janitor is not less biased than the detective. The janitor is biased in a different direction — toward whatever is visually salient. Logical thinking is the chain-of-custody protocol that applies to both.**

### Second-most-common confusion (worth naming)

> **"Absence of evidence is evidence of absence."**

A non-significant result is not evidence that there is no effect. It is evidence that *this design, at this sample size, with this measurement, could not distinguish the effect from noise* — which may mean the effect is small, or the assay is noisy, or the study is underpowered, or the assay failed. Without a **positive control**, a null result cannot even distinguish "no biology" from "broken experiment."

### What this distinction is NOT

Careful — overcorrecting is its own failure mode. The distinctions above do **not** mean:

- **Not** "hypothesis-driven is always better." Hypothesis-driven research with a wrong hypothesis produces confident nulls and misses real biology. Data-driven research discovers things hypothesis-driven work cannot anticipate.
- **Not** "correlational research is worthless." Well-designed, replicated, confounder-adjusted associational evidence is often the best available evidence, and sometimes the only ethical option.
- **Not** "exploratory work is illegitimate." Exploration is how every field generates its next generation of hypotheses. The error is mislabeling, not exploring.
- **Not** "pre-registration guarantees validity." A pre-registered confounded design is still confounded.
- **Not** "statistical significance is meaningless." It is meaningful *given* a valid design, a correct null, and a declared error budget. It is meaningless when transplanted into a design that cannot support it.
- **Not** "skepticism means rejecting claims." Skepticism means asking what would count for and against, then checking whether the study was built to answer that.
- **Not** "rigor is bureaucracy." Rigor is the set of artifacts that let a stranger audit your reasoning. If it feels like paperwork, that is often because the reasoning was never made explicit.

> **The overcorrection to avoid:** swinging from "any significant result is a discovery" to "no result is trustworthy unless it comes from a randomized prospective trial." The first is credulity; the second is paralysis. Logical thinking is the middle path: **match the strength of the claim to the strength of the design.**

### Decision rule (preview of Phase 4)

Before you reach for a conclusion, ask which of these you actually have:

| If your evidence is… | Then your strongest defensible claim is… |
|---|---|
| A pattern in one dataset | "We observed an association in this dataset." |
| A replicated pattern across independent datasets | "The association replicates." |
| A replicated pattern + confounder adjustment + dose–response | "A robust association, consistent with a causal effect." |
| Exogenous variation (MR, natural experiment) + defended assumptions | "Evidence of a causal effect under stated assumptions." |
| Randomized perturbation + replication | "Causal in this system." |

Any claim stronger than its row is an overreach — and overreach is a **logical** error, not a statistical one.

---

Ready for Phase 4? Any questions or wonderings so far?