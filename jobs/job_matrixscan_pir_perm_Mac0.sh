#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N matrixscan_pir_perm_Mac0
#$ -o logs/matrixscan_pir_perm_Mac0.o
#$ -e logs/matrixscan_pir_perm_Mac0.e
#$ -m e
#$ -M tu_correo@dominio.com
. /etc/profile.d/modules.sh
module load rsat/8sep2021

rsat matrix-scan   -quick \ 
  -i "data/PIR_Mac0.fa"   -bgfile "data/bg_Mac0.txt"   -m "data/Mac0_perm.txt"   -matrix_format transfac   -bg_pseudo 0.01   -lth score "5"   -return sites,pval,rank   -o "results/quick_pir_perm_Mac0.txt"
