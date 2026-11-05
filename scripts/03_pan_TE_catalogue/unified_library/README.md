# Unified curated repeat library

## Purpose

This workflow merges one or more genome-specific curated repeat libraries into a nonredundant unified consensus library. The same script is used regardless of the number or identity of input libraries; actual inputs are supplied through a runtime TSV manifest and are not committed to the repository.

## Method

1. Read every FASTA listed in the manifest.
2. Convert sequences to uppercase and calculate MD5 hashes.
3. Deduplicate exact sequence matches only. No CD-HIT or near-identity clustering is performed.
4. Reconcile classifications for exact duplicate sequences:
   - identical fine-level classifications: retain the original classification;
   - different fine-level labels but the same broad TE class: assign the corresponding broad-class fallback (`DNA/Unknown`, `LTR/Unknown`, `LINE/Unknown`, `SINE/Unknown` or `RC/Helitron`);
   - one known TE classification versus `Unknown`, with exactly one known fine-level classification: retain the known classification (`KNOWN_VS_UNKNOWN_KEEP_KNOWN`);
   - agreement on nonTE or Unknown: retain that status;
   - other conflicting broad classes: assign `Unknown` and record the conflict.
5. Sort unique sequences by MD5 before assigning unified IDs, making the output ordering deterministic for the same input sequence set.

## Input manifest

Copy `input_libraries.tsv.example` to a local manifest and replace the generic paths. The manifest must be tab-separated and contain exactly these required columns:

```text
source_id    library_fasta
```

`source_id` is a generic source label used only for provenance in the output mapping table. Add one row for each input library. Relative FASTA paths are resolved relative to the manifest location.

## Run

```bash
python3 scripts/03_pan_TE_catalogue/unified_library/01_build_unified_library.py \
    --manifest /path/to/input_libraries.tsv \
    --output-dir /path/to/unified_library_output \
    --output-prefix Unified_v5_2 \
    --id-prefix UNI_v52
```

Python 3 standard library only; no third-party Python package is required.

## Outputs

| Output | Description |
|---|---|
| `Unified_v5_2.curated_repeat.fa` | Nonredundant unified consensus library |
| `Unified_v5_2.sequence_map.tsv` | Unified ID, sequence MD5, final class, broad class, length, source-record count, decision and source-member provenance |
| `Unified_v5_2.summary.tsv` | Input counts, unique sequence count, duplicate groups, classification decisions and final broad-class counts |
| `Unified_v5_2.class_conflicts.tsv` | Exact sequence groups with broad-class disagreements, including known-versus-Unknown cases; the sequence map records the final resolution |

The FASTA and TSV outputs are generated from the supplied input libraries. Input FASTAs and production-specific manifests are not included in this repository.