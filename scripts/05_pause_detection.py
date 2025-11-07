#!/usr/bin/env python3
import yaml, pandas as pd, os
config = yaml.safe_load(open('config.yaml'))
proj = config['project_dir']
psite = pd.read_csv(os.path.join(proj,"step03_riboQC","psite_table.tsv"), sep='\t')
# Expected psite_table has columns: sample, transcript, position (nt)
# Map position to codon within transcript (assuming CDS coordinates available).
# Simple codon occupancy z-score by transcript:
psite['codon_pos'] = (psite['position'] // 3).astype(int)
grp = psite.groupby(['transcript','codon_pos']).size().reset_index(name='counts')
grp['mean'] = grp.groupby('transcript')['counts'].transform('mean')
grp['std'] = grp.groupby('transcript')['counts'].transform('std').replace(0,1)
grp['z'] = (grp['counts'] - grp['mean']) / grp['std']
pauses = grp[grp['z'] > 3]
os.makedirs(os.path.join(proj,"step05_pause"), exist_ok=True)
pauses.to_csv(os.path.join(proj,"step05_pause","pause_sites.tsv"), sep='\t', index=False)
print("Pause detection done. Found", len(pauses), "sites")

## This is a simple reproducible caller — for publication use RUST or more robust methods.
