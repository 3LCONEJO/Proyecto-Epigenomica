#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N matrixscan_pir_perm_Mac1
#$ -o logs/matrixscan_pir_perm_Mac1.o
#$ -e logs/matrixscan_pir_perm_Mac1.e
#$ -m e
#$ -M tu_correo@dominio.com
. /etc/profile.d/modules.sh
module load rsat/8sep2021

rsat matrix-scan   -i "data/PIR_Mac1.fa"   -bgfile "data/bg_Mac1.txt"   -m "data/Mac1_H3K27me3_perm.txt"   -matrix_format transfac   -bg_pseudo 0.01   -lth score "5"   -return sites,pval,rank   -o "results/quick_pir_perm_Mac1.txt"
