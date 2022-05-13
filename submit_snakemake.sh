#!/usr/bin/bash
#SBATCH --time 35:00:00
#SBATCH --partition exacloud
#SBATCH --qos normal
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log

# Use with cluster.json.qos. You will need to change the name to cluster.json.

snakemake -j 750 --rerun-incomplete --use-conda --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N} -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem-per-cpu {cluster.mem} --qos {cluster.qos} --exclude='exanode-4-15,exanode-4-30,exanode-6-27,exanode-6-47,exanode-2-44,exanode-8-[22-24],exanode-3-[0-6],exanode-1-[32-44],exanode-0-[26-29]'" -s Snakefile --latency-wait 1000 --keep-going
