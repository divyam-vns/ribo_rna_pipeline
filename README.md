# ribo_rna_pipeline
# 1. Ribo/RNA Joint Pipeline (RiboCode + DESeq2 TE + Pause detection)

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

## Run pipeline stepwise: 
   ```bash
bash scripts/00_download_geo.sh config.yaml
bash scripts/01_trim_and_rRNA_removal.sh config.yaml
bash scripts/02_align_star.sh config.yaml
Rscript scripts/03_psite_riboWaltz.R config.yaml
bash scripts/04_orf_detection.sh config.yaml
python3 scripts/05_pause_detection.py config.yaml
bash scripts/06_counts_featureCounts.sh config.yaml
Rscript scripts/07_deseq2_TE.R config.yaml

## Open analysis/RiboSeq_DESeq2_TE_Analysis.Rmd for plots and summary.

**Notes**

Configure config.yaml carefully (paths and indexes).
Recommended >= 3 biological replicates.
RiboCode and riboWaltz outputs used by downstream steps.
