# Reproducibility records

For each analysis, record:

- the biological or technical question;
- input files, sample set, provenance and checksums;
- output files and formats;
- software names, exact versions and installation sources;
- complete commands and parameter settings;
- reference assembly, annotation, TE-library release and coordinate convention;
- Conda environment or container image and digest, where applicable;
- execution date, logs, QC criteria and exclusions.

Use configurable input and output paths rather than hard-coded user-specific paths. Validate required inputs, report software versions and provide a usage example. Keep the exact script used for manuscript results distinguishable from exploratory or superseded scripts.

Individual, tissue, sequencing library and sequencing run are distinct identifiers. Technical libraries from one individual are not independent biological replicates. For coordinate-based results, record the assembly/haplotype and whether coordinates are 0-based half-open or 1-based closed.

The templates under `metadata/templates/` are blank examples. They must be populated from the actual data records; their presence does not establish that an analysis was performed.

Report source-level audits, static syntax checks, synthetic tests, data-dependent computational reruns and experimental validation separately. A passing syntax check does not establish numerical correctness, and a checksum does not establish biological validity.
