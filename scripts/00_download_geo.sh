#!/usr/bin/env bash
# Usage: bash 00_download_geo.sh config.yaml
CONFIG=$1
if [ -z "$CONFIG" ]; then echo "Usage: $0 config.yaml"; exit 1; fi
python - <<PY
import yaml, os, subprocess
cfg=yaml.safe_load(open("$CONFIG"))
out=cfg['fastq_dir']
os.makedirs(out, exist_ok=True)
# Example GEO/SRR list - replace with actual SRR IDs or extend to parse GEO
srrs = ["SRRXXXXXX","SRRYYYYYY"]
for s in srrs:
    cmd = ["prefetch", s]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd)
    cmd2 = ["fasterq-dump", s, "-O", out, "--split-files", "--gzip"]
    subprocess.run(cmd2)
PY
echo "Download done (edit script to set actual SRR IDs)."


##Edit srrs with the runs you want (retrieve from GEO accession pages).
