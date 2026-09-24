# xQTL AD Loci Explorer

Code and released tables behind the **AD Loci Explorer**, a browsable view of the
ADSP FunGen-xQTL release across 188 candidate Alzheimer's disease loci.

Live app: https://jenny-empawi.shinyapps.io/xQTL-AD-loci-Explore/

## Layout

```
scripts/
  build_AD_locus_table.R   locus-level build: evidence integration and tier assignment
  gene_prio_utils.R        helper functions, including the T1-T6 tier rules
  add_evidence_source.R    register a new evidence source (see 'Adding data')
  validate_outputs.R       check a release against the expectations of this build
run_pipeline.qsub          end-to-end job script for the SCC scheduler
config/                    metadata tables the build reads
data/                      released tables, one directory per build
app/                       Shiny application, runnable from a clone
archive/                   superseded and ad-hoc scripts; not part of the pipeline
DEPENDENCIES.md            R version and package versions
```

## Requirements

R 4.5.3 with the packages listed in DEPENDENCIES.md:

```r
install.packages(c("data.table","stringr","openxlsx","tidyverse","ggplot2",
                   "scales","circlize","shiny","bslib","DT","plotly",
                   "shinycssloaders"))
remotes::install_github("StatFunGen/pecotmr")   # not on CRAN
```

## Running it

```bash
export AD_LOCI_ROOT=/path/to/AD_loci_xQTL      # where the input trees live
export AD_LOCI_STAGING=/path/to/staging        # three precomputed tables, see Inputs
export AD_LOCI_OUT=$AD_LOCI_ROOT/out_$(date +%Y%m%d)    # optional

Rscript scripts/build_AD_locus_table.R                 # -> $AD_LOCI_OUT
Rscript scripts/validate_outputs.R  $AD_LOCI_OUT       # expected tables, 188 loci, 508 genes
Rscript app/build_shiny_data.R $AD_LOCI_OUT app/data.csv
```

The build checks every required input before starting and stops with a list of
whatever is missing, rather than failing part-way through a long run.

## Inputs

`config/metadata_analysis.csv` is the registry of every input the build reads.
One row per exported analysis table:

| column | meaning |
|---|---|
| context | dataset label, joined to config/contexts_metadata.csv |
| Data Type | eQTL, sQTL, pQTL, ... |
| Cohort | ROSMAP, MSBB, KNIGHT, ... |
| Modality | assay or analysis modality |
| Method | how the table was produced; selects the loader in the build |
| Path | location of the table, relative to `$AD_LOCI_ROOT` |
| summary_file, summary_file_ad | written by the build at run time; leave blank |
| variant_level_method | whether the method contributes variant-level evidence |

Paths under `$AD_LOCI_ROOT` refer to the FunGen-xQTL release, Synapse
`syn68872650`. The layout on Synapse does not necessarily mirror these relative
paths; email jaempawi@bu.edu if you need help locating a particular table. Some
rows carry `<user>` and `<collaborator>` placeholders where
tables were exported from per-user analysis directories; replace these with the
directory names in your own copy.

Three precomputed tables are too large to distribute with the code and are read
from `$AD_LOCI_STAGING`:

| file | size | source |
|---|---|---|
| gwas_variants_cor0.5.csv.gz | 295 MB | Synapse `syn75082260` |
| res_APOE_interaction_summ.csv.gz | 45 MB | on request |
| res_msex_interaction_summ.csv.gz | 1.6 MB | on request |

The two interaction summaries are not redistributed publicly. To request them,
open an issue on this repository or email jaempawi@bu.edu.

The correlation table drives the LD-based credible-set extension, so the
released variant sets cannot be reproduced without it.

## Adding data

To add a new evidence source you register it in `config/metadata_analysis.csv`.
The build auto-detects any Method there that points at per-context toploci files
and has no dedicated loader, so registration is the only step.

For a ColocBoost-style colocalization summary, `scripts/add_evidence_source.R` does the
conversion and registration for you:

```bash
Rscript scripts/add_evidence_source.R --name transmap \
        --bed trans_xQTL_only_colocalization_summary_table.bed
Rscript scripts/build_AD_locus_table.R
```

It is safe to re-run: conversion is skipped when outputs already exist, and rows
for the same `--name` are replaced rather than duplicated.

### Input format

The `--bed` file is tab-delimited with a header, one row per variant, and event
lists delimited by `;`:

| # | column | notes |
|---|---|---|
| 1 | #chr | |
| 2 | start | |
| 3 | end | |
| 4 | a1 | |
| 5 | a2 | |
| 6 | variant_ID | chr:pos:ref:alt |
| 7 | region_ID | |
| 8 | event_ID | semicolon-delimited list; context is the part before `_cis_` or `_trans_` |
| 9 | cos_ID | must contain `cosN`; rows with `cos0` are skipped |
| 10 | vcp | becomes PIP in the emitted file |
| 11 | cos_npc | becomes purity |
| 14 | zscore | semicolon-delimited, positionally matched to event_ID |

It emits one toploci file per context with the columns the build expects:
`#chr start end a1 a2 variant_ID region_ID event_ID cs_coverage_0.95 purity PIP z`.

A source in any other format needs its own loader in `scripts/build_AD_locus_table.R`;
follow one of the existing `Method ==` blocks.

## Confidence tiers

`scripts/gene_prio_utils.R` assigns `top_confidence` per row during the build. T1-T5
describe genes with localized AD-xQTL support, ordered by strength of evidence.
**T6** covers genes supported only by gene-level TWAS/MR or cTWAS evidence, with
no localized xQTL signal; because the tier chain is evaluated over rows of the
xQTL overlap table, such a gene never reaches the final branch, so the build
assigns T6 directly from the gene-level tables produced earlier in the same run.

A gene is reported at its strongest tier. The app takes its tiers from the
release, so the table and the Explorer cannot diverge.

## Running the app

```r
shiny::runApp("app")
```

`app/data.csv` ships with the repository, so the Explorer runs without the
pipeline. Note that `app/build_shiny_data.R` carries a number of gene-level
columns forward from the existing `app/data.csv`, so that file is an input as
well as an output.

## Data

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.
