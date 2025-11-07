#!/usr/bin/env bash
# Run RiboCode (example). Adjust RiboCode call to your installation.
CONFIG=$1
python - <<PY
import yaml, pandas as pd, os, subprocess
cfg=yaml.safe_load(open("$CONFIG"))
proj=cfg['project_dir']
psite = os.path.join(proj,"step03_riboQC","psite_table.tsv")
outdir = os.path.join(proj,"step04_ribocode")
os.makedirs(outdir, exist_ok=True)
# RiboCode command (change according to RiboCode CLI on your system)
cmd = ["RiboCode","-i",cfg['gtf'],"-r",psite,"-o",outdir]
print("Running:", " ".join(cmd))
subprocess.run(cmd)
PY


## If RiboCode is not available, use RiboTaper or RiboTISH. Outputs: translated ORFs list.
