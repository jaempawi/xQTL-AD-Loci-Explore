# xQTL AD Loci Explorer

Code and released tables behind the **AD Loci Explorer**, a browsable view of the
ADSP FunGen-xQTL release across 188 candidate Alzheimer's disease loci.

Live app: https://jenny-empawi.shinyapps.io/xQTL-AD-loci-Explore/

## Layout

```
app/          Shiny application (runnable as-is)
pipeline/     Scripts that build the locus-level tables
  staging/    Metadata and helper functions the main script reads
  jobs/       SGE submission scripts for the cluster
  dev/        Ad-hoc inspection scripts; not pipeline steps
data/188_loci/  Released locus- and variant-level tables for this build
docs/         Deployment notes
```

## Requirements

R with `data.table`, `stringr`, `openxlsx`, `shiny`, `ggplot2`, and `pecotmr`;
Python 3 with `pandas`. The Excel summary step uses `openxlsx`.

## Running the pipeline

Both path roots are environment variables, so the pipeline is not tied to one
filesystem:

```bash
export AD_LOCI_ROOT=/path/to/AD_loci_xQTL     # where the input trees live
export AD_LOCI_OUT=$AD_LOCI_ROOT/out_$(date +%Y%m%d)   # optional; defaults to out_<today>
```

Stage order:

1. `pipeline/fetch_from_hpc.sh` — stage input tables locally.
2. `pipeline/add_evidence_source.R` — register an evidence source in
   `staging/metadata_analysis.csv`. Registering a source there is the only step
   needed to add one; the main script picks it up automatically.
3. `pipeline/complete_ADlocus_level_summary_fixed.R` — the main build. Reads the
   metadata, assembles fine-mapping, colocalization, TWAS/MR and cTWAS evidence
   per variant and gene, assigns confidence tiers T1-T6, and writes the
   per-variant tables plus the unified Excel summary into `$AD_LOCI_OUT`.
4. `pipeline/tier_assign_202609.py` and `pipeline/classify_confidence.py` —
   gene-level tier and confidence assignment over the tables from step 3.
5. `pipeline/coloc_validate.py` — colocalization sanity checks.
6. `app/build_shiny_data.R` — collapses the release into `app/data.csv`, the
   single table the Explorer reads.

### Confidence tiers

T1-T5 describe genes with localized AD-xQTL support, ordered by strength of
evidence. **T6** covers genes supported only by gene-level TWAS/MR or cTWAS
evidence, with no localized xQTL signal. Because the tier chain is evaluated
over rows of the xQTL overlap table, a gene with no localized support never
reaches the final branch; the main script therefore assigns T6 directly from
the gene-level XWAS/MR and cTWAS tables produced earlier in the same run.

## Running the app

```r
shiny::runApp("app")
```

`app/data.csv` is included, so the Explorer runs from a clone without the
pipeline.

## Data

`data/188_loci/` holds the locus summary and the unified variant-level table
for this build. Upstream inputs — GWAS fine-mapping exports, xQTL
colocalization results and the LD reference — are not redistributed here; they
are released through the AD Knowledge Portal on Synapse (`syn68872650`).

Derived from ADSP/NIAGADS study data. Downstream use of the underlying
individual-level and controlled-access resources is governed by their own data
use terms.
