# Pipeline

See the repository README for stage order and environment variables.

`staging/` holds the metadata tables and `gene_prio_utils.R`, which the main
script sources. `metadata_analysis.csv` is the registry of input tables: each
row points at an exported analysis table and the Method that produced it.

## Paths in this released copy

Absolute paths have been replaced with `$AD_LOCI_ROOT`, and two placeholders
appear in `staging/metadata_analysis.csv`:

- `<user>` and `<collaborator>` in the `Path` column, where tables were exported
  from per-user analysis directories. Substitute the directory names for your
  own layout before running.
- The `summary_file` and `summary_file_ad` columns are blank. They are
  populated at run time with paths under `$AD_LOCI_OUT` and are not input
  configuration.

The `Path` column is relative to `$AD_LOCI_ROOT`. The upstream tables it points
at are released through Synapse, not in this repository.

`dev/` contains ad-hoc scripts used to inspect intermediate tables. They are not
part of the pipeline and are kept only for reference.
