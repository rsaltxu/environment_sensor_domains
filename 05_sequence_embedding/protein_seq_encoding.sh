#!/bin/bash

#SBATCH -J HK_seq            # Job name
#SBATCH -o myjob_%j.o       # Name of stdout output file(%j expands to jobId)
#SBATCH -e myjob_%j.e       # Name of stderr output file(%j expands to jobId)
#SBATCH --gres=gpu:a100:1   # Request 1 GPU of 2 available on an average A100 node
#SBATCH -c 32               # Cores per task requested
#SBATCH -t 5-00:00:00         # Run time (hh:mm:ss) - 10 min
#SBATCH --mem-per-cpu=2G    # Memory per core demandes (96 GB = 3GB * 32 cores)

conda activate esmfold
python $HOME/scripts/encoders/protein_seq_encoding_test.py $STORE/UpSet_HKs.fasta $STORE/UpSetHK_seq_encoding.csv --format csv
echo "done"                 # Write this message on the output file when finished