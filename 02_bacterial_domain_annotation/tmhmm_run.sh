#!/bin/bash
#SBATCH -o myjob_%j.o       # Name of stdout output file(%j expands to jobId)
#SBATCH -e myjob_%j.e       # Name of stderr output file(%j expands to jobId)
#SBATCH -c 5 #number of cores
#SBATCH --mem-per-cpu=3952M
#SBATCH -t 1-00:00:00 # execution time
#SBATCH --mail-type=begin #Envía un correo cuando el trabajo inicia
#SBATCH --mail-type=end #Envía un correo cuando el trabajo finaliza

tmhmm /data/tmhmm//HK_db_BLASTp_ids.fasta