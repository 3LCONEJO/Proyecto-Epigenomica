#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N matrixscan_prom_orig_Mac0
#$ -o logs/matrixscan_prom_orig_Mac0.o
#$ -e logs/matrixscan_prom_orig_Mac0.e
#$ -m e
#$ -M tu_correo@dominio.com
. /etc/profile.d/modules.sh
module load rsat/8sep2021

rsat matrix-scan   -quick \ 
  -i "data/promotores_Mac0.fa"   -bgfile "data/bg_Mac0.txt"   -m "data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt"   -matrix_format transfac   -bg_pseudo 0.01   -lth score "5"   -return sites,pval,rank   -o "results/quick_prom_orig_Mac0.txt"
