# xQTL-AD-Loci-Explore

Interactive Shiny app for browsing xQTL evidence at Alzheimer's disease GWAS loci.

**Live app:** https://jenny-empawi.shinyapps.io/xQTL-AD-loci-Explore/

---

## Overview

This app visualises results from the FunGen-xQTL consortium ADSP Functional Genomics pipeline, linking AD GWAS variants to molecular QTL evidence across brain cell types.

### Data contents (current build)
- **4,817 rows** — one per variant × gene pair with xQTL evidence
- **Sources:** Brain eQTL · Brain pQTL · Brain gpQTL · Excitatory / Inhibitory / Oligodendrocyte / OPC / Astrocyte / Microglia / Bulk Immune eQTL · MiGA eQTL (GFM, GTS, SVZ, THA)
- **Trans xQTL:** 2,263 loci with at least one trans-acting gene

### Features
- Filter by chromosome, GWAS study, significance level, confidence level, cell type, and functional evidence (TWAS / MR / cTWAS)
- Toggle individual cell-type columns on/off
- Side detail panel per variant — shows evidence contexts, cell type breakdown, and trans xQTL target genes
- Download filtered results as CSV

---

## Repository structure

| File/Folder | Description |
|-------------|-------------|
| `app.R` | Shiny app source (R) |
| `data.csv` | Processed display data (generated from flat file pipeline) |
| `DEPLOY.md` | Step-by-step deployment instructions |
| `data/` | Jupyter notebook + Excel summary table |

---

## Updating the data

Re-run `AD_loci_xQTL_table_updated.ipynb` on the cluster with the latest flat file, then re-deploy to shinyapps.io:

```r
rsconnect::deployApp(appDir = ".", appName = "xQTL-AD-loci-Explore", account = "jenny-empawi")
```

See `DEPLOY.md` for full instructions.

---

## Citation / Contact

FunGen-xQTL Consortium · ADSP Functional Genomics
