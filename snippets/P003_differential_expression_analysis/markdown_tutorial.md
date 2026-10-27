# Tutorial Summary: Markdown in VSCode for NGS/Omics

> From plain-text notes to a fully reproducible, Git-tracked documentation system that lives alongside your Snakefile and config — with safe, portable, and citable Markdown.
>
> **Date:** 2026-10-08  
> **Author:** Racey Ye  
> **Tutor:** DeepSeek

---

## 1. Learning Log by Phase

### Phase 1: Historical Context

**Key ideas:**
- Markdown was created by John Gruber and Aaron Swartz in **2004** to make web writing readable as plain text while still converting to HTML.
- Before Markdown: HTML was verbose, WYSIWYG editors produced bloated non-portable output.
- Markdown was intentionally vague, leading to many flavors (CommonMark, GitHub Flavored Markdown, MultiMarkdown).
- VSCode (2015) became the editor of choice because of built-in preview, Git integration, extensions, and terminal access.

**My questions:**
- "Is Markdown a simpler method for displaying code snippets with annotations, which was once only achievable via annotations or screenshots?"
  - **Answer:** Yes. The annotation pattern (numbered comments in code plus an explanation list below) reproduces the pedagogical effect of screenshots, while staying copy-pasteable, diffable, and searchable.

**My wonderings:**
- How does Markdown handle complex documents that used to require LaTeX?
  - **Answer:** Markdown with MathJax (via extensions) handles light math. For publication-quality output, Markdown is the source and Pandoc converts to LaTeX/PDF.

---

### Phase 2: Deep Grammatical Breakdown

**Key syntax learned:**

| Feature | Syntax | Example |
|---------|--------|---------|
| Heading | `# Title` (space required) | `## QC` |
| Bold | `**text**` or `__text__` | `**fastp**` |
| Italic | `*text*` or `_text_` | `*fastp*` |
| Bullet list | `- item` (space required) | `- FastQC` |
| Ordered list | `1. item` (space required) | `1. fastp` |
| Task list | `- [ ] item` / `- [x] done` | `- [x] Trim reads` |
| Inline code | `` `code` `` | `` `samtools sort` `` |
| Fenced code | ` ```lang ... ``` ` | ` ```bash ... ``` ` |
| Link | `[text](url)` | `[STAR](https://...)` |
| Image | `![alt](path)` | `![QC](qc.png)` |
| Table | `\| col \| col \|` | see Phase 3 |
| Blockquote | `> text` | `> Warning` |
| Horizontal rule | `---` | section separator |
| Escaping | `\*` | `\*not italic\*` |
| Footnote | `text[^1]` | GFM only |
| Mermaid | ` ```mermaid ... ``` ` | pipeline diagrams |

**Key pattern — annotated code:**

````markdown
```bash
fastp \
  -i sample_R1.fq.gz \   # (1)
  -I sample_R2.fq.gz \   # (2)
  -o clean_R1.fq.gz \    # (3)
  -O clean_R2.fq.gz \    # (4)
  --thread 8             # (5)
```

1. Input R1
2. Input R2
3. Output R1 after trimming
4. Output R2 after trimming
5. Use 8 threads
````

**My questions:**
- "Should I use backticks instead of double quotes in bash for integrity?"
  - **Answer:** No. This is a category error. Backticks in bash are **command substitution** (execution), not quoting. Double quotes `"$var"` are the safety mechanism. Markdown backticks and bash double quotes are unrelated. See Section 3.

**My wonderings:**
- How do I keep code and prose separate?
  - **Answer:** Fenced code blocks for commands, tables for parameters, lists for explanation. Never mix.

---

### Phase 3: Point-by-Point Divergence

**Key comparisons:**

| Feature A | Feature B | Rule |
|-----------|-----------|------|
| `*text*` (italic) | `**text**` (bold) | One star = italic, two = bold |
| `# Heading` | `#Heading` | Space is mandatory |
| `- item` | `-item` | Space is mandatory |
| `\|` in table | escaped pipe | Must escape pipes |
| `*.fastq` in prose | `` `*.fastq` `` | Wrap in backticks |
| `1.` auto-numbered | Manual numbering | Markdown renumbers automatically |
| Inline link | Reference link | Reference needs `[ref]: url` defined |

**Common mistakes identified:**

1. **Missing space after marker** — `#Heading`, `-item`, `>quote` all become literal text.
2. **Unclosed fenced code blocks** — the rest of the document is swallowed.
3. **Bare code-like text** — `*.fastq`, `sample_id`, `__init__` get mangled by emphasis.
4. **Unescaped pipes** — break table columns.
5. **Assuming VSCode preview equals GitHub rendering** — different feature sets.

**My questions:**
- Why does `__init__` break in Markdown?
  - **Answer:** `__...__` is bold syntax. Wrap in backticks: `` `__init__` ``.

**My wonderings:**
- How does Markdown distinguish Glob `*` from Regex `*`?
  - **Answer:** It doesn't. Both are literal text unless wrapped in backticks. The distinction only matters when you interpret them (filesystem vs. text).

---

### Phase 4: Extended Usage & Real-World Application

**Key examples learned:**
- NGS repo layout (`README.md`, `CHANGELOG.md`, `docs/`, `notes/`, `config/`, `envs/`, `workflow/`)
- `README.md` template with Mermaid pipeline diagram
- Sample sheet as Markdown table (`docs/samples.md`)
- Annotated command documentation (`docs/pipeline.md`)
- Lab notebook entries (`notes/YYYY-MM-DD.md`)
- Methods section for a paper (`docs/methods.md` plus Pandoc export)
- Snakemake rule docstrings linking to Markdown anchors
- VSCode snippets, tasks, extensions for NGS Markdown
- QC report inline documentation with images

**My questions:**
- "Could everything in bioinformatics be recorded as Markdown, including code snippets, methods, version control?"
  - **Answer:** Markdown covers documentation, methods, notes, and rationale. It does **not** replace YAML for configuration, TSV for sample sheets, Snakefile for workflow logic, or Git for version control. Each layer answers a different question.

- "If I set Snakemake as the pipeline controller, is describing the Snakefile itself enough?"
  - **Answer:** No. The Snakefile describes **what runs and in what order**, not **why** parameters were chosen, **what** samples are, **which** tool versions are needed, or **how** to install and run. You need README, config, environment.yml, and docs.

**My wonderings:**
- How do YAML, JSON, Markdown, and Snakefile relate?
  - **Answer:** YAML configures; Markdown explains; Snakefile orchestrates; Git tracks everything as text.

- Can Markdown replace LaTeX for publication?
  - **Answer:** Only via Pandoc conversion. Markdown is the source; LaTeX/PDF is the output.

- Should I use backticks in bash instead of double quotes?
  - **Answer:** No. Category error. See Section 3.

---

## 2. Mastery Check: My Answers and Corrections

**Final Score: 61 / 100**

| Q | My Answer | Correct Answer | Explanation / Correction |
|---|-----------|----------------|--------------------------|
| Q1 | B | B | Correct |
| Q2 | B | B | `Ctrl+K V` = side-by-side preview |
| Q3 | B | B | Space after `#` is mandatory |
| Q4 | B | B | `\|` or backticks for literal pipes |
| Q5 | C | C | Markdown explains; YAML configures; Git tracks |
| Q6 | ` ``` ` | ` ``` ` | Correct |
| Q7 | `**` or `__` | `**` or `__` | Accepted |
| Q8 | whitespace | space | Accepted |
| Q9 | even | even | Fence parity |
| Q10 | quotes `''` | backticks `` ` `` | Critical error. See Section 3. |
| Q11 | partial | full | Regex `*` is a quantifier, not "any char" |
| Q12 | added dots | wrap in backticks | Separator uses `-` and `:`, not `.` |
| Q13 | correct | correct | Executable vs explanatory truth |
| Q14 | prose + 2 blocks | 1 block + annotations | Missed annotation pattern |
| Q15 | partial | full lab entry | Missing Objective, Results table, Task list, Wonderings |
| Q16 | tree + brief | full mapping table | Missing explicit file-to-role mapping |

---

## 3. Mistakes & Misunderstandings

### Mistake 1: Quotes vs. Backticks in Markdown Prose

- **What I thought:** Quotes (`''` or `""`) protect underscores from Markdown emphasis.
- **What is actually true:** Only **backticks** suppress Markdown interpretation. Quotes are literal characters with no special power.
- **Correct mental model:** Backticks = "this is code, don't touch it." Quotes = "this is a quotation."

### Mistake 2: Regex `*` Semantics

- **What I thought:** In regex, `*` means "any single character."
- **What is actually true:** In regex:
  - `.` = any single character
  - `*` = zero or more of the **preceding** element (a quantifier)
  - `.*` = any character, zero or more times
- **Correct mental model:** Glob `*` is greedy filename matching. Regex `*` is a quantifier on the previous token. Markdown `*` is emphasis. Three different `*`s.

### Mistake 3: Table Separator Syntax

- **What I thought:** Dots (`.`) could serve as separator characters.
- **What is actually true:** The separator row uses **hyphens** (`-`) and optional **colons** (`:`) only.
- **Correct mental model:** `|:---|:---:|---:|` = left, center, right alignment.

### Mistake 4: Annotated Code Pattern

- **What I thought:** Prose + separate code blocks = annotation.
- **What is actually true:** Annotation = **numbered comments inside the code block** + **numbered list below**.
- **Correct mental model:** The code block stays runnable; annotations live outside it.

### Mistake 5: Lab Notebook Completeness

- **What I thought:** A title, hypothesis, and command are enough.
- **What is actually true:** A complete entry needs Title, Objective, Hypothesis, Commands, Results table, Interpretation, Task list, Wonderings.
- **Correct mental model:** A lab notebook entry answers "what, why, what happened, what next."

### Mistake 6: Bash Backticks vs. Double Quotes

- **What I thought:** Backticks in bash serve the same "protect this" role as backticks in Markdown.
- **What is actually true:** In bash, backticks are **command substitution** — they **execute** the enclosed command and substitute its output. They are the opposite of "protecting." Double quotes `"$var"` are the safety mechanism.
- **Correct mental model:** Markdown backticks = display. Bash backticks = execute (deprecated, use `$(...)`). Bash double quotes = safety.

---

## 4. Open Wonderings & Next Steps

**Wonderings:**
- How does `Ctrl+K V` differ from `Ctrl+Shift+V`?
  - **Answer:** `Ctrl+K V` = side-by-side preview. `Ctrl+Shift+V` = preview replaces current tab.
- How do Mermaid diagrams render across VSCode, GitHub, GitLab?
  - **Answer:** All three support Mermaid, but styling differs. Preview on target platform before publishing.
- When to use `samples.tsv` vs. embedding sample info in `config.yaml`?
  - **Answer:** `samples.tsv` for tabular, one-row-per-sample data. `config.yaml` for global parameters and metadata that applies across samples.
- Should I use backticks in bash instead of double quotes?
  - **Answer:** No. Different languages, different purposes. See Section 3.

**Topics to revisit:**
- Pandoc export workflows (Markdown to PDF, DOCX, LaTeX).
- Quarto for reproducible reports (Markdown + code + output in one document).
- Mermaid diagrams for pipeline visualization.
- Pre-commit hooks for Markdown linting.
- `markdownlint` rule set customization.

**Suggested next exercises:**
1. Write a `README.md` for an existing NGS project using the Phase 4 template.
2. Convert a `methods.md` to PDF via Pandoc.
3. Add a pre-commit hook that checks fence parity in all `.md` files.
4. Create a VSCode snippet for the annotated code block pattern.
5. Write a full lab notebook entry with all seven sections.
6. Rewrite 5 sentences with bare `sample_id`, `*.fastq`, `__init__` using proper backticks.

---

## 5. Golden Rules & Cheat Sheet

### Golden Rules

> **Markdown Rule 1:** A marker without a space is just text. `#`, `-`, `>`, `1.` all need a following space.
>
> **Markdown Rule 2:** Count your fences. Every ` ``` ` must have a partner. Odd count = broken document.
>
> **Markdown Rule 3:** Wrap code-like text in backticks. Globs, regex, variable names, snake_case identifiers, dunder methods, filenames.
>
> **Bash Rule 1:** Double-quote your variables: `"${var}"`. Unquoted variables word-split and glob-expand.
>
> **Bash Rule 2:** Backticks in bash are **execution**, not quoting. Use `$(...)` for command substitution.
>
> **Separation Rule:** Markdown explains. YAML configures. Snakefile orchestrates. TSV tabulates. Git tracks all of them.
>
> **Documentation Rule:** Write once in Markdown. Export to PDF/DOCX/HTML via Pandoc. Never edit exports directly.
>
> **Reproducibility Rule:** Every command in a README must be copy-pasteable and runnable.

### Cheat Sheet — Markdown Syntax

| Element | Syntax |
|---------|--------|
| Heading | `# H1` ... `###### H6` |
| Bold | `**text**` |
| Italic | `*text*` |
| Inline code | `` `code` `` |
| Fenced code | ` ```bash ... ``` ` |
| Link | `[text](url)` |
| Image | `![alt](path)` |
| Bullet | `- item` |
| Numbered | `1. item` |
| Task | `- [ ] todo` / `- [x] done` |
| Table | `\| a \| b \|` |
| Blockquote | `> text` |
| Rule | `---` |
| Escape | `\*` |
| Footnote | `text[^1]` |
| Mermaid | ` ```mermaid ... ``` ` |

### Cheat Sheet — Bash Quoting

| Syntax | Behavior | Use |
|--------|----------|-----|
| `$var` | Splits, globs (dangerous) | Almost never |
| `"$var"` | Safe expansion | **Default** |
| `'$var'` | Literal, no expansion | Fixed strings |
| `` `cmd` `` | Executes (deprecated) | Never |
| `"$(cmd)"` | Executes + quotes | Preferred |

### Cheat Sheet — NGS Project Layout

```text
project/
├── README.md              ← Front door
├── CHANGELOG.md           ← Version history
├── Snakefile              ← Workflow logic
├── config/
│   ├── config.yaml        ← Parameters
│   └── samples.tsv        ← Sample sheet
├── envs/
│   └── environment.yml    ← Tool versions
├── docs/
│   ├── pipeline.md        ← Walkthrough
│   ├── methods.md         ← Paper methods
│   └── troubleshooting.md ← Known issues
├── notes/
│   └── YYYY-MM-DD.md      ← Lab notebook
└── results/
    └── figures/           ← Referenced from docs
```

### Cheat Sheet — VSCode Markdown

| Action | Shortcut |
|--------|----------|
| Side-by-side preview | `Ctrl+K V` |
| Preview in current tab | `Ctrl+Shift+V` |
| Bold | `Ctrl+B` |
| Italic | `Ctrl+I` |
| Insert link | `Ctrl+K` |
| Toggle task | `Alt+C` |
| Command palette | `Ctrl+Shift+P` |
| Find across `.md` | `Ctrl+Shift+F` |

---

## 6. Final Mental Model

Markdown is the **lab notebook** of bioinformatics: a plain-text, Git-friendly, human-readable format that documents *what* you did and *why*, using simple symbols (`#`, `-`, backtick, `|`) that render to HTML, PDF, or DOCX on demand. It is the **explanatory truth** sitting alongside the **executable truth** of your Snakefile and the **configuration truth** of your YAML. VSCode is its natural habitat: live preview, Git integration, extensions for Mermaid and Pandoc, and snippets that turn repetitive documentation into one-keystroke patterns. The three rules that keep you safe are: (1) **space after every marker**, (2) **backticks around every code-like token** (globs, regex, snake_case, dunders, filenames), (3) **even number of code fences**. And never confuse Markdown backticks (display) with bash backticks (execution) or bash double quotes (safety). Each language has its own rules; the bioinformatics engineer's job is to keep them in sync and never mix their syntaxes.

---

**End of tutorial summary.**