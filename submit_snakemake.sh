#!/usr/bin/bash
#SBATCH --time 35:00:00
#SBATCH --partition exacloud
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log

snakemake -j 100 --rerun-incomplete --use-conda --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem {cluster.mem}" -s Snakefile --latency-wait 60
