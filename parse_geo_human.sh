#!/bin/bash -l

#$ -P bf528         # REQUIRED: Replace with your SCC project
#$ -l h_rt=12:00:00             # Max runtime (12 hours)
#$ -pe omp 8                    # Request X cores
#$ -l mem_per_core=2G           # memory per core
#$ -N parse_geo_batch            # Job name
#$ -j y                         # Merge stdout/stderr
#$ -o log.txt                   # Output log 
#$ -cwd                         # Run job from current working directory
#$ -m bea           # Send email at (b)egin and (e)nd
#$ -M yuansh24@bu.edu  # Replace with your actual BU email

# --- Setup ---
module load R/4.4.0             # Load your preferred R version

# Optional: Print basic info
echo "Job ID: $JOB_ID"
echo "Running on host: $(hostname)"
echo "Using $NSLOTS cores and $((NSLOTS * 2))G memory"
echo "Current working directory: $(pwd)"
echo "Current date and time: $(date)"
echo "R version: $(R --version)"

# --- Measure resources used and run script---
/usr/bin/time -v Rscript parse_geo_batch.R gse_human.txt output/human

