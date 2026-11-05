# Homologous gene-pair identification

Reciprocal best-hit protein pairs between the BT-M and BT-P haplotypes are identified using MMseqs2 `easy-rbh`. The historical settings are retained: sensitivity 7.5, E-value 1e-5, minimum sequence identity 0.5, coverage 0.6 and `--cov-mode 0`. The tabular output includes query/target IDs, identity, alignment and sequence lengths, query/target coverage, E-value and bit score.

The reported count of 21,017 pairs is a study result and should be re-confirmed against the final RBH output before it is included in the release summary. MMseqs2 version should be recorded from the original environment.
