# Chromosome-wise phylogenetic validation

The analysis is performed independently for chr01–chr24 using seven taxon proteomes: Mamb, Ppek, Calb, Cery, Cide, BT-M and BT-P.

Workflow: OrthoFinder single-copy orthologue identification (DIAMOND search); MUSCLE protein alignment; Gblocks conserved-region filtering with the original settings `-b4=5 -b5=h -t=p`; SeqKit sorting and FASTA formatting; concatenation of aligned sequences; manual taxon-ID standardization and matrix inspection; ProtTest 3.4 model assessment; RAxML 8.2.8 maximum-likelihood inference.

**Supermatrix construction:** the original SeqKit step sorts each `*.muscle-gb` alignment with `seqkit sort`, converts it to single-line FASTA with `seqkit seq -w 0`, and combines the aligned files with `paste -d ' '` followed by removal of spaces (`sed -i 's/ //g'`) to generate chromosome-level `all.fa`. Taxon sequence IDs are then manually standardized/calibrated, and the resulting matrix is checked for consistency and suitability before downstream phylogenetic analysis. This manual calibration and matrix QC are intentional steps; no separate concatenation software is required.

ProtTest settings retained from the original command: `-all-distributions -F -AIC -BIC -tc 0.5 -threads 36`. RAxML settings retained: `-f a -x 12345 -p 12345 -# 1000 -m PROTGAMMAJTTF`, with the chromosome-specific Cide outgroup label.

All 24 chromosomes, including chr01, follow the same general analysis procedure. The original scripts had a loop boundary that could start at chr02 for ProtTest/RAxML; for the released workflow, the chromosome range is consistently chr01–chr24 as specified for this analysis.

Set the paths and executable variables in the scripts before submission. The `all.fa` matrix must pass manual taxon-ID calibration and inspection before it is renamed/converted to the format required by ProtTest and RAxML.