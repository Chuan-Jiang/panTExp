#!/bin/bash
cd /lustre/research/dawli/chujiang/TE_exp/00_Codes
git add -A
git commit -m "Auto-backup: $(date '+%Y-%m-%d %H:%M:%S')" || echo "No changes to commit"
git push origin HEAD