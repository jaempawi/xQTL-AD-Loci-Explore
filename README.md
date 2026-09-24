# xQTL AD Loci Explorer

Code and released tables behind the **AD Loci Explorer**, a browsable view of the
ADSP FunGen-xQTL release across 188 candidate Alzheimer's disease loci.

Live app: https://jenny-empawi.shinyapps.io/xQTL-AD-loci-Explore/

## Layout

```
build_AD_locus_table.R   locus-level build: evidence integration and tier assignment
gene_prio_utils.R        helper functions, including the T1-T6 tier rules
validate_outputs.R       checks a release against the expectations of this build
inputs_manifest.tsv      every input the build reads, and where it comes from
config/                  metadata tables the build reads
data/188_loci/           released locus- and variant-level tables
app/                     Shiny application, runnable from a clone
archive/                 superseded and ad-hoc scripts; not part of the pipeline
DEPENDENCIES.md          R version and package versions
```

## Running it

Two stages. Set the roots first:

```bash
export AD_LOCI_ROOT=/path/to/AD_loci_xQTL        # where the input trees live
export AD_LOCI_STAGING=/path/to/staging          # three precomputed tables, see below
export AD_LOCI_OUT=$AD_LOCI_ROOT/out_$(date +%Y%m%d)   # optional
```

```bash
Rscript build_AD_locus_table.R                   # -> $AD_LOCI_OUT
Rscript validate_outputs.R  $AD_LOCI_OUT         # sanity-check the release
Rscript app/build_shiny_data.R $AD_LOCI_OUT app/data.csv
```

The build checks every required input before it starts and stops with the list
of anything missing, rather than failing part-way through a long run.

## Inputs

`inputs_manifest.tsv` lists each input, the path it is read from, and its source.
Paths under `AD_LOCI_ROOT` come from the FunGen-xQTL release on Synapse
(`syn68872650`). Two entries in `config/metadata_analysis.csv` carry `<user>` and
`<collaborator>` placeholders where tables were exported from per-user analysis
directories; substitute your own layout.

Three tables are too large to distribute with the code and are read from
`AD_LOCI_STAGING`:

| file | source |
|---|---|
| `gwas_variants_cor0.5.csv.gz` | Synapse `syn75082260` |
| `res_APOE_interaction_summ.csv.gz` | available on request |
| `res_msex_interaction_summ.csv.gz` | available on request |

The correlation table drives the LD-based credible-set extension, so the released
variant sets cannot be reproduced without it.

## Confidence tiers

`gene_prio_utils.R` assigns `top_confidence` per row during the build. T1-T5
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
pipeline.

## Data

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.
