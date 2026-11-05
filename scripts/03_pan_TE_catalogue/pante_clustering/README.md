# Pan-TE clustering scripts

## Requirements

- BLAST+ / BLASTn v2.15.0+ (`makeblastdb` and `blastn` available in `PATH`)
- Python 3 (standard library only)
- SLURM for the all-versus-all BLASTn array

## Working directory

Set `PANTE_WORKDIR` to an absolute path for the Pan-TE working directory. The scripts use `00_input/` for the five class-specific FASTA files, `01_blast/` for BLAST databases and all-versus-all results, and `02_cluster/` for initial and refined cluster tables and QC files.

## Execution order

Run these commands from the repository root, replacing input paths with those in your environment:

```bash
export PANTE_WORKDIR=/path/to/PanTE_work

python scripts/03_pan_TE_catalogue/pante_clustering/01_split_unified_by_class.py \
    --input /path/to/Unified_v5_2.curated_repeat.fa \
    --output-dir "$PANTE_WORKDIR/00_input"

cd scripts/03_pan_TE_catalogue/pante_clustering
sbatch --export=ALL 02_run_pante_all_vs_all.sh
cd ../../..

python scripts/03_pan_TE_catalogue/pante_clustering/03_build_pante_clusters.py \
    --workdir "$PANTE_WORKDIR"

# Optional: graph QC of initial components (>=10 members)
python scripts/03_pan_TE_catalogue/pante_clustering/04_audit_large_clusters.py \
    --workdir "$PANTE_WORKDIR"

python scripts/03_pan_TE_catalogue/pante_clustering/05_refine_pante_representative_clusters.py \
    --workdir "$PANTE_WORKDIR"
```

Submit the SLURM array from this script directory so its `.out` and `.err` files are written beside the script. The array has five tasks, one per broad TE class.

## Script descriptions

| Script | Function |
|---|---|
| `01_split_unified_by_class.py` | Splits the unified curated consensus library into DNA, LTR, LINE, SINE and RC FASTA files. |
| `02_run_pante_all_vs_all.sh` | Runs within-class all-versus-all BLASTn searches. |
| `03_build_pante_clusters.py` | Integrates HSPs, applies high-confidence edge criteria and constructs initial connected components. |
| `04_audit_large_clusters.py` | Optional graph QC of initial components containing at least 10 members; does not alter cluster membership. |
| `05_refine_pante_representative_clusters.py` | Refines initial components using representative-centered greedy clustering. |

The original source-script names and public-release changes are recorded in [../SOURCE_SCRIPT_INVENTORY.md](../SOURCE_SCRIPT_INVENTORY.md).
