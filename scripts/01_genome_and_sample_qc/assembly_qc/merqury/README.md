# Merqury assembly-evaluation scripts

These path-portable scripts implement the two commands used for k-mer-based assembly evaluation.

## Software

- Meryl 1.4.1
- Merqury 1.4.1

Both executables must be available in `PATH` on the compute node. The scripts contain no cluster-specific installation paths, user names, or private project directories.

## Workflow

1. Generate a k=21 count database from PacBio HiFi reads with `01_meryl_k21_count.sh`.
2. Run Merqury with the Meryl database and the final assembly FASTA using `02_merqury_assembly_evaluation.sh`.

The Meryl command used 36 threads. The first script requests 36 CPUs per task, matching the recorded Meryl command. The second script retains a 10-CPU Slurm request, consistent with the submitted Merqury job configuration; the Merqury command itself is unchanged.

## Required environment variables

| Script | Required variables | Optional variables |
|---|---|---|
| `01_meryl_k21_count.sh` | `HIFI_READS` | `MERYL_DB`, `THREADS` |
| `02_merqury_assembly_evaluation.sh` | `MERYL_DB`, `ASSEMBLY_FASTA` | `OUTPUT_PREFIX` |

The assembly FASTA used for the reported Merqury run was checked and confirmed to be the final assembly FASTA. The whole-assembly QV output is retained at [`result.qv`](../../../../01_genome_and_sample_qc/assembly_qc/merqury/results/result.qv); its QV value is 55.6649 (55.66 in the manuscript). Chromosome-level QV values and k-mer completeness are already reported in Supplementary Table S3 and are not duplicated as separate output tables here.
