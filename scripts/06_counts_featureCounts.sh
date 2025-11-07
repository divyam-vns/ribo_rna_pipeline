#!/usr/bin/env bash
CONFIG=$1
python - <<PY
import yaml, pandas as pd, os, subprocess
cfg=yaml.safe_load(open("$CONFIG"))
proj=cfg['project_dir']; samples=pd.read_csv(cfg['sample_table'], sep='\\t')
rna_bams=[]; ribo_bams=[]
for i,row in samples.iterrows():
  base_rna = os.path.basename(row['rna_fastq']).replace('.fastq.gz','')
  base_ribo= os.path.basename(row['ribo_fastq']).replace('.fastq.gz','')
  rna_bam = os.path.join(proj,'step02_star',base_rna,'Aligned.sortedByCoord.out.bam')
  ribo_bam= os.path.join(proj,'step02_star',base_ribo,'Aligned.sortedByCoord.out.bam')
  rna_bams.append(rna_bam); ribo_bams.append(ribo_bam)
outdir=os.path.join(proj,'step06_counts'); os.makedirs(outdir, exist_ok=True)
# RNA gene-level
cmd1 = ["featureCounts","-T", str(cfg['threads']), "-a", cfg['gtf'], "-o", os.path.join(outdir,"rna_counts.txt")] + rna_bams
print("Running:", " ".join(cmd1)); subprocess.run(cmd1)
# Ribo counts on CDS
cmd2 = ["featureCounts","-T", str(cfg['threads']), "-t", "CDS", "-g", "gene_id", "-a", cfg['gtf'], "-o", os.path.join(outdir,"ribo_counts.txt")] + ribo_bams
print("Running:", " ".join(cmd2)); subprocess.run(cmd2)
PY
