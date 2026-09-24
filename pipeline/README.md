# Pipeline

See the repository README for stage order and environment variables.

`staging/` holds the metadata tables and `gene_prio_utils.R`, which the main
script sources. `metadata_analysis.csv` is the registry of input tables: each
row points at an exported analysis table and the Method that produced it.

`dev/` contains ad-hoc scripts used to inspect intermediate tables. They are not
part of the pipeline and are kept only for reference.
