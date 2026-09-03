# Sustained MYC activation reveals head-constrained transcriptome distortion

This is the reviewer-facing analysis repository for the manuscript **“Sustained MYC activation reveals head-constrained transcriptome distortion”** (iScience manuscript **ISCIENCE-D-26-10474**). It contains compact processed inputs where available, analysis and figure-generation code, source-data tables or accession-level provenance records, and manuscript-numbered display files required to audit the reported results.

The repository is organized by final manuscript display rather than as a general-purpose software package. [`FIGURE_MAP.tsv`](FIGURE_MAP.tsv) is the authoritative map from manuscript panels to code, inputs, source data, and outputs.

## Scientific overview

The study uses a controlled Tet3G/MYC–ER perturbation system to examine how MYC-associated transcription changes when activation is sustained. Three linked questions organize the analysis:

1. **Does sustained MYC output remain proportional across the ranked transcriptome?**  
   Distributional and rank-space analyses distinguish early near-parallel behavior from a late state in which proportionality is selectively attenuated at the expression head. This late state is termed **transcriptome distortion**.

2. **Can the late state be explained by the engineered background?**  
   Dox-only, Tam-only, priming-response, promoter-motif, and fixed promoter-context analyses distinguish the Dox-established baseline—termed **doxycycline priming**—from distortion.

3. **What analytical behavior accompanies the biological transition?**  
   Independent readouts converge on the late state, while differential-expression frameworks diverge in a reproducible, expression-head-localized manner.

The repository supports the reported state definition, localization, and sensitivity analyses. It does not claim to identify a single molecular bottleneck responsible for distortion.

## Reviewer quick start

```bash
git clone https://github.com/pop3www/GMS_Ranking.git
cd GMS_Ranking
git checkout v1.0.0-iscience-review
python3 scripts/validate_repo.py --root . --strict
```

For an audit without rerunning the full workflow, begin with:

```text
README.md           scientific overview, repository map, and all run instructions
final_outputs/      final manuscript-numbered figures
FIGURE_MAP.tsv      panel → code → input → source data → output
README.md           execution and provenance notes
README.md           GEO and processed-data access
*/source_data/      bundled source tables underlying the displays
source_data/GEO_TRACE.tsv  GEO traceability when a processed table is not duplicated
```

All reviewer navigation is consolidated in this root README.

## Experimental conditions and principal contrasts

| Label | Experimental state |
|---|---|
| `Ctrl` | Untreated matched control |
| `Tam` | Tamoxifen-triggered MYC–ER activation without Dox priming |
| `D` | Dox/Tet3G background with MYC–ER OFF |
| `DT` | Tamoxifen-triggered MYC–ER activation with Dox held constant |
| `D_L_CPT`, `DT_L_CPT` | Dox-matched MYC–ER OFF/ON pair under low CPT |
| `D_H_CPT`, `DT_H_CPT` | Dox-matched MYC–ER OFF/ON pair under high CPT |

Principal MYC activation contrasts hold the Dox background constant:

```text
DT vs D
DT_L_CPT vs D_L_CPT
DT_H_CPT vs D_H_CPT
```

Dox-only and Tam-only controls are represented by `D vs Ctrl` and `Tam vs Ctrl`. Analyses compare early (~4 h) and sustained (~24 h) responses.

## Final figure map

| Final display | Primary module(s) | Final output |
|---|---|---|
| Figure 1 | `figure1/` | `final_outputs/Figure1.pdf` |
| Figure 2 | `figure2/` | `final_outputs/Figure2.pdf` |
| Figure 3A–C and Figure S5 | `figure3/differential_expression/` | component of `final_outputs/Figure3.pdf`; `FigureS5.pdf` |
| Figure 3D–E | `figure3/scalar_metrics/` | component of `final_outputs/Figure3.pdf` |
| Figure 3 composite | `figure3/final_assembly/` | `final_outputs/Figure3.pdf` |
| Figure 4 | `figure4/` | `final_outputs/Figure4.pdf` |
| Figure 5 and Figure S2 | `figure5/` | `Figure5.pdf`; `FigureS2.pdf` |
| Figures 6–7 and Figure S3 | `figure6_7/` | `Figure6.pdf`; `Figure7.pdf`; `FigureS3.pdf` |
| Figure S1 | `supplement/figureS1/` | `final_outputs/FigureS1.pdf`; standardized inputs traced in `source_data/GEO_TRACE.tsv` |
| Figure S4 | `supplement/figureS4/` | `final_outputs/FigureS4.pdf` |

Use [`FIGURE_MAP.tsv`](FIGURE_MAP.tsv) for panel-level provenance and exact entry points.

## Analysis inputs and provenance

The downstream analyses start from compact RSEM-derived gene-level matrices in `data/processed/`. RSEM expected counts form the unified gene-by-sample expression matrix. For DESeq2, edgeR, and limma/voom, non-negative expected counts were rounded before method-specific normalization and testing. RankProd used logCPM derived from the same rounded matrix. PenDA and RankCompV3 results used in the integrated comparison are included as frozen, precomputed outputs with manifests documenting provenance and contrast orientation.

Method-specific signed outputs are rank-normalized within method before cross-method Spearman concordance and inter-method spread calculations. The **baseline-expression slope, β** is the robust regression coefficient relating gene-wise log2 fold change to matched-baseline expression. Negative values indicate attenuation or inversion of baseline-scaled proportionality—not “negative amplification.” Gene-span sensitivity analyses use GENCODE transcription-unit spans calculated as gene end minus gene start plus one.

The primary sample-level audit uses ridge-logistic models with feature selection performed within training folds. Leave-one-sample-out evaluation measures sample-held-out separation; leave-one-context-out evaluation tests transfer to a withheld CPT context.

## Data access

The standardized datasets underlying this study are deposited in NCBI GEO:

- **Bulk RNA-seq:** [GSE318584](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318584)
- **ChIP-seq:** [GSE318271](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318271)
- **SuperSeries:** [GSE318586](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318586)

Only the bulk RNA-seq and ChIP-seq datasets listed above underlie this manuscript. Additional records linked within the SuperSeries were not analyzed here. If the GEO records remain private during peer review, anonymous reviewer tokens are supplied confidentially through Editorial Manager and are never committed to GitHub.

Large raw sequencing files and public genome references are not versioned in this repository. For Figure S1, the processed promoter-profile tables are not duplicated in the concise release; the ChIP-seq source is traced to GSE318271 and the expression-rank source to GSE318584 in `supplement/figureS1/source_data/GEO_TRACE.tsv`.

## Figure S1 data boundary

The final Supplementary Figure 1 display is included in `final_outputs/FigureS1.pdf`. The exact processed input-normalized promoter-profile tables are not bundled in this concise release and are therefore **not claimed as bundled source data**. Instead, `supplement/figureS1/source_data/GEO_TRACE.tsv` records the complete accession-level provenance: MYC/RNAPII ChIP-seq is traced to **GSE318271**, and RNA-seq baseline-rank assignment is traced to **GSE318584** within SuperSeries **GSE318586**. If the GEO records remain private during review, anonymous reviewer tokens are supplied confidentially through Editorial Manager and are never stored in GitHub. Reproduction of Figure S1 begins from these standardized GEO datasets.

## Repository organization

| Path | Purpose |
|---|---|
| `core/` | Shared expression-matrix, metadata, baseline-rank, validation, and plotting utilities |
| `figure1/` | Distributional and Δ-density analyses |
| `figure2/` | Reagent controls, Dox priming, FIMO motif scans, and Figure 2 assembly |
| `supplement/figureS1/` | Fixed input-normalized MYC/RNAPII promoter-context analysis; GEO provenance is recorded when processed profiles are not bundled |
| `figure3/differential_expression/` | Unified DE comparison, signed scores, method spread, and Figure S5 concordance |
| `figure3/scalar_metrics/` | Up fraction and baseline-expression slope used in final Figure 3D–E |
| `figure3/final_assembly/` | Final merged Figure 3 assembly |
| `figure4/` | Rank-shift geometry and Dox-background concordance used in final Figure 4 |
| `figure5/` | Sample-level classifier audits and Figure S2 controls used in final Figure 5 |
| `figure6_7/` | Module, enrichment, leading-edge, and mitochondrial analyses used in Figures 6–7 and Figure S3 |
| `supplement/figureS2/` | Navigation and execution wrapper for the classifier controls generated with Figure 5 |
| `supplement/figureS3/` | Navigation and execution wrapper for the mitochondrial analysis generated with Figures 6–7 |
| `supplement/figureS2/` | Navigation and execution wrapper for the classifier controls generated with Figure 5 |
| `supplement/figureS3/` | Navigation and execution wrapper for the mitochondrial analysis generated with Figures 6–7 |
| `supplement/figureS4/` | Gene-span sensitivity analysis used in Figure S4 |
| `supplement/figureS5/` | Pairwise method-concordance analysis generated from the Figure 3 differential-expression workflow |
| `supplement/figureS5/` | Pairwise method-concordance analysis generated from the Figure 3 differential-expression workflow |
| `config/` | Frozen sample sheets, contrasts, and parameters |
| `data/processed/` | Compact processed matrices used by downstream analyses |
| `final_outputs/` | Manuscript-numbered review PDFs |
| `scripts/` | Validation and top-level execution |

## Reproduction

Run all commands from the repository root. First validate the release:

```bash
python3 scripts/validate_repo.py --root . --strict
```

The validator checks repository structure, code syntax, private paths or credentials, manuscript-numbered outputs, shared-code links, and release checksums.

Rebuild the downstream workflow from bundled processed inputs after installing the documented R and Python dependencies. Figure S1 additionally requires retrieval of its standardized ChIP-seq and RNA-seq inputs from GEO as documented above:

```bash
bash scripts/run_all.sh list
# To execute all canonical workflows after installing dependencies:
bash scripts/run_all.sh execute
```

Canonical module entry points are consolidated here:

| Final display or task | Command from repository root |
|---|---|
| Prepare/validate shared expression inputs | `bash core/scripts/install_expression_inputs.sh` |
| Figure 1 | `bash figure1/run.sh` |
| Figure 2 | `bash figure2/run.sh` |
| Figure S1 provenance check | `bash supplement/figureS1/run.sh` |
| Figure 3 | `bash figure3/run.sh` |
| Figure S5 | `bash supplement/figureS5/run.sh` |
| Figure 4 | `bash figure4/run.sh` |
| Figure 5 and Figure S2 | `bash figure5/run.sh` |
| Figures 6–7 and Figure S3 | `bash figure6_7/run.sh` |
| Figure S4 | `bash supplement/figureS4/run.sh` |

Some module scripts require the environment and public references described in `env/`, `config/`, and this README. The top-level `scripts/run_all.sh` is the authoritative execution order.

## What is intentionally excluded

Raw FASTQ/BAM/BigWig files and large public reference resources are not stored in this repository. The upstream processing protocol is described in the manuscript Methods. Downstream analyses use bundled compact matrices; where a processed table is not bundled, an accession-level provenance manifest links the analysis to GEO.

## Release integrity

A complete review release contains `Figure1.pdf` through `Figure7.pdf` and `FigureS1.pdf` through `FigureS5.pdf` in `final_outputs/`, together with `SHA256SUMS.txt`. The exact release tag and commit hash should match the manuscript, response letter, and Editorial Manager access field.

## AI-assisted code development

ChatGPT (OpenAI) was used in a supportive capacity for code troubleshooting and documentation. All suggestions were reviewed, executed, and validated by the authors against the underlying data and expected outputs. Final analytical decisions and biological interpretations were made by the authors.

## Citation and contact

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Questions concerning the analysis repository may be directed through the corresponding-author contact given in the manuscript.
