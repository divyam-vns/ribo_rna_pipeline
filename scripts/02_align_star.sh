#!/usr/bin/env bash
CONFIG=$1
python - <<PY
import yaml, pandas as pd, os, subprocess
cfg=yaml.safe_load(open("$CONFIG"))
samples = pd.read_csv(cfg['sample_table'], sep='\\t')
outbase = os.path.join(cfg['project_dir'], "step02_star")
os.makedirs(outbase, exist_ok=True)
for idx,row in samples.iterrows():
  for col in ['rna_fastq','ribo_fastq']:
    base = os.path.basename(row[col]).replace('.fastq.gz','')
    fq = os.path.join(cfg['project_dir'], "step01_trim", base + ".norrna.fastq.gz")
    outdir = os.path.join(outbase, base)
    os.makedirs(outdir, exist_ok=True)
    cmd = [
      "STAR","--runThreadN",str(cfg['threads']), "--genomeDir",cfg['star_genome_dir'],
      "--readFilesIn", fq, "--readFilesCommand","zcat",
      "--outFileNamePrefix", outdir + "/", "--outSAMtype","BAM","SortedByCoordinate",
      "--outFilterMultimapNmax","1","--alignEndsType","EndToEnd"
    ]
    print("Running STAR:", " ".join(cmd))
    subprocess.run(cmd)
print("STAR alignment done")
PY
