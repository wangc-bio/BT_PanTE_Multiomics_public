# Synteny analysis

Four gene sets are compared pairwise: BSB, TC, BT-M and BT-P, producing six pairwise comparisons. GFF3 annotations are converted to BED using mRNA features and ID attributes; duplicate BED entries are removed; corresponding CDS sequences are extracted. JCVI/MCscan is used for orthology/synteny analysis. Anchor filtering retains blocks with `--minspan=30`; the historical analysis used C-score 0.80 in the MCscan stage. The original plot uses `jcvi.graphics.karyotype` with the supplied `seqids` and `layout` configuration files, 12 × 10 inch canvas, PDF/SVG output and 600 dpi.

The six pairs are BSB–BT-P, BSB–BT-M, BSB–TC, BT-M–BT-P, BT-M–TC and BT-P–TC. Provide the study-specific `seqids` and `layout` files alongside this workflow when publishing the corresponding figure. JCVI version should be recorded from the actual analysis environment; it is not inferred here.

Scripts require the input GFF3/CDS files and a JCVI conda environment. Set the documented environment variables before submission.
