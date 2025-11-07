# ribo_rna_pipeline
# Ribo/RNA Joint Pipeline (RiboCode + DESeq2 TE + Pause detection)

This repo performs end-to-end analysis for matched Ribo-seq + RNA-seq:
- Download GEO raw FASTQs (Singh et al., Elkon et al. examples)
- Adapter trimming, rRNA removal
- STAR alignment
- P-site assignment (riboWaltz)
- ORF calling (RiboCode)
- Pause-site detection (codon z-score method)
- Gene-level counts and DESeq2 interaction model for differential Translation Efficiency (TE)
- QC plots and output TE lists

## Quickstart
1. Edit `config.yaml` and `samples.tsv`.
2. Create conda env:
   ```bash
   conda env create -f environment.yml
   conda activate ribo_env
