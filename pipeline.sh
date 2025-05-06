#Uso: ./pipeline.sh <data_table> <cell_type> <hg19_fa_url|ruta> <bigWig_file> <histone_mark> <jaspar_matrix> <score_threshold>
# --------------------------------------------------------
set -euo pipefail
if [ $# -lt 7 ]; then
  echo "Uso: $0 <data_table> <cell_type> <hg19_fa_url|ruta> <bigWig_file> <histone_mark> <jaspar_matrix> <score_threshold>" >&2
  exit 1
fi

# Parámetros
DATA_TABLE=$1
CELL_TYPE=$2
HG19_FASTA=$3
BIG_WIG=$4
HISTONE=$5
JASPAR_MATRIX=$6
SCORE_THRES=$7

# Crear directorios si no existen
mkdir -p data results logs jobs

# Función de progreso con hélice de ADN animada
show_progress() {
  local pid=$1
  local task_name=$2
  local delay=0.2

  tput civis
  echo -ne "\033[s"

  local frames=()

  # Define cada frame individualmente para evitar el problema de sintaxis
  frames[0]="
      ..........0..........
      .........o-0.........
      ........o---0........
      .......o-----0.......
      .......0-----0.......
      .......0-----o.......
      ........0---o........
      .........0-o.........
      ..........0.........."
      
  frames[1]="
      .........o-0.........
      ........o---0........
      .......o-----0.......
      .......0-----0.......
      .......0-----o.......
      ........0---o........
      .........0-o.........
      ..........0..........
      .........o-0........."
      
  frames[2]="
      ........o---0........
      .......o-----0.......
      .......0-----0.......
      .......0-----o.......
      ........0---o........
      .........0-o.........
      ..........0..........
      .........o-0.........
      ........o---0........"
      
  frames[3]="
      .......o-----0.......
      .......0-----0.......
      .......0-----o.......
      ........0---o........
      .........0-o.........
      ..........0..........
      .........o-0.........
      ........o---0........
      .......o-----0......."
      
  frames[4]="
      .......0-----0.......
      .......0-----o.......
      ........0---o........
      .........0-o.........
      ..........0..........
      .........o-0.........
      ........o---0........
      .......o-----0.......
      .......0-----0......."
      
  frames[5]="
      .......0-----o.......
      ........0---o........
      .........0-o.........
      ..........0..........
      .........o-0.........
      ........o---0........
      .......o-----0.......
      .......0-----0.......
      .......0-----o......."
      
  frames[6]="
      ........0---o........
      .........0-o.........
      ..........0..........
      .........o-0.........
      ........o---0........
      .......o-----0.......
      .......0-----0.......
      .......0-----o.......
      ........0---o........"
      
  frames[7]="
      .........0-o.........
      ..........0..........
      .........o-0.........
      ........o---0........
      .......o-----0.......
      .......0-----0.......
      .......0-----o.......
      ........0---o........
      .........0-o........."

  local n_frames=${#frames[@]}

  while kill -0 "$pid" 2>/dev/null; do
    for ((i=0; i<n_frames; i++)); do
      echo -ne "\033[u\033[H\033[J"
      printf "%s\n" "${frames[i]}"
      printf "\n[%s] En progreso...\n" "$task_name"
      sleep "$delay"
    done
  done

  tput cnorm
  echo -ne "\033[u\033[K[$task_name] ✓ Completado\n"
}

# Función para generar scripts de job de matrix-scan
generate_matrixscan_job() {
  local seq_fa=$1      # Archivo FASTA de entrada
  local bgfile=$2      # Archivo background
  local matrix=$3      # Matriz JASPAR o permutada
  local output=$4      # Archivo de salida
  local desc=$5        # Descripción amigable
  local job_name=$6    # Nombre del job

  cat <<EOF > jobs/job_${job_name}.sge
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N ${job_name}
#$ -o logs/${job_name}.o
#$ -e logs/${job_name}.e
#$ -m e
#$ -M tu_correo@dominio.com
. /etc/profile.d/modules.sh
module load rsat/8sep2021

rsat matrix-scan \
  -quick \
  -i "${seq_fa}" \
  -bgfile "${bgfile}" \
  -m "${matrix}" \
  -matrix_format transfac \
  -bg_pseudo 0.01 \
  -lth score "${SCORE_THRES}" \
  -return sites,pval,rank \
  -o "results/${output}"
EOF

  chmod +x jobs/job_${job_name}.sge
}

# 1) Extracción de BED de promotores y PIRs
echo "==> Extrayendo datos de $DATA_TABLE para $CELL_TYPE"
Rscript Parse_Data.R "$DATA_TABLE" "$CELL_TYPE" & pid=$!; show_progress $pid "Extracción de datos"

# 2) Descarga o uso local del genoma hg19.fa
if [[ "$HG19_FASTA" =~ ^https?:// ]]; then
  echo "==> Descargando genoma desde $HG19_FASTA"
  wget --timestamping "$HG19_FASTA" & pid=$!; show_progress $pid "Descarga del genoma"
  GENOME=$(basename "$HG19_FASTA" .gz)
  if [[ "$HG19_FASTA" == *.gz ]]; then
    echo "==> Descomprimiendo ${GENOME}.gz"
    gunzip -f "${GENOME}.gz" & pid=$!; show_progress $pid "Descompresión del genoma"
  fi
else
  GENOME="$HG19_FASTA"
  echo "==> Usando genoma local: $GENOME"
fi

PROM_BED="data/promotores_${CELL_TYPE}.bed"
PIR_BED="data/PIR_${CELL_TYPE}.bed"
sed -i 's/^chrMT\t/chrM\t/' "$PROM_BED" "$PIR_BED"

# 3) Extracción de secuencias FASTA
echo "==> Extrayendo secuencias FASTA"
bedtools getfasta -fi "$GENOME" -bed "$PROM_BED" -fo "data/promotores_${CELL_TYPE}.fa" & pid=$!; show_progress $pid "Extracción promotores FASTA"
bedtools getfasta -fi "$GENOME" -bed "$PIR_BED" -fo "data/PIR_${CELL_TYPE}.fa" & pid=$!; show_progress $pid "Extracción PIRs FASTA"

# Normalizar a mayúsculas
echo "==> Normalizando secuencias"
tr 'acgt' 'ACGT' < data/promotores_${CELL_TYPE}.fa > tmp.fa && mv tmp.fa data/promotores_${CELL_TYPE}.fa
tr 'acgt' 'ACGT' < data/PIR_${CELL_TYPE}.fa > tmp.fa && mv tmp.fa data/PIR_${CELL_TYPE}.fa

cat data/promotores_${CELL_TYPE}.fa data/PIR_${CELL_TYPE}.fa > data/combined_${CELL_TYPE}.fa

# 4) Señal con bigWigAverageOverBed
echo "==> Calculando señal (${HISTONE})"
bigWigAverageOverBed "$BIG_WIG" "$PROM_BED" data/prom_${CELL_TYPE}.tab & pid=$!; show_progress $pid "Señal promotores"
bigWigAverageOverBed "$BIG_WIG" "$PIR_BED" data/pir_${CELL_TYPE}.tab & pid=$!; show_progress $pid "Señal PIRs"

# 5) Gráficas con Rscript
echo "==> Generando gráficas"
Rscript Prueba_diferenciacion.R "$CELL_TYPE" "$HISTONE" > results/wilcox_results.txt & pid=$!; show_progress $pid "Generación de gráficas"

# 6) Permutación de matriz JASPAR
echo "==> Permutando matriz JASPAR"
permute-matrix -i "$JASPAR_MATRIX" -o data/${CELL_TYPE}_perm.txt -in_format transfac -out_format transfac & pid=$!; show_progress $pid "Permutación matriz"

# 7) Modelo background común
echo "==> Generando background común"
oligo-analysis -l 1 -i data/combined_${CELL_TYPE}.fa -format fasta -o data/bg_${CELL_TYPE}.txt -type dna & pid=$!; show_progress $pid "Background common"

BG=data/bg_${CELL_TYPE}.txt

# 8) Generación y envío de jobs matrix-scan
echo "==> Creando y enviando jobs de matrix-scan"
declare -a outputs

for stage in prom_orig prom_perm pir_orig pir_perm; do
  case "$stage" in
    prom_orig)
      input="data/promotores_${CELL_TYPE}.fa"; matrix="$JASPAR_MATRIX"; desc="promotores orig";;
    prom_perm)
      input="data/promotores_${CELL_TYPE}.fa"; matrix="data/${CELL_TYPE}_perm.txt"; desc="promotores perm";;
    pir_orig)
      input="data/PIR_${CELL_TYPE}.fa"; matrix="$JASPAR_MATRIX"; desc="PIR orig";;
    pir_perm)
      input="data/PIR_${CELL_TYPE}.fa"; matrix="data/${CELL_TYPE}_perm.txt"; desc="PIR perm";;
  esac
  output="quick_${stage}_${CELL_TYPE}.txt"
  outputs+=("results/${output}")
  job_name="matrixscan_${stage}_${CELL_TYPE}"
  generate_matrixscan_job "$input" "$BG" "$matrix" "$output" "$desc" "$job_name"
  qsub jobs/job_${job_name}.sge
done

# 9) Esperar a que los archivos de salida existan y no estén vacíos
echo "Esperando a que los 4 archivos de matrix-scan estén listos..."
while :; do
  all_done=true
  for f in "${outputs[@]}"; do
    if [ ! -s "$f" ]; then
      all_done=false
      break
    fi
  done
  if $all_done; then
    echo "Todos los archivos de matrix-scan están disponibles."
    break
  fi
  sleep 30
done

# 10) Análisis de enriquecimiento vía R
echo "==> Análisis de enriquecimiento de TFs"
Rscript enrichment_analysis.R \
  data/promotores_${CELL_TYPE}.fa \
  results/quick_prom_orig_${CELL_TYPE}.txt \
  results/quick_prom_perm_${CELL_TYPE}.txt \
  data/PIR_${CELL_TYPE}.fa \
  results/quick_pir_orig_${CELL_TYPE}.txt \
  results/quick_pir_perm_${CELL_TYPE}.txt \
  results/enrichment_${CELL_TYPE}

echo "  → Resultados guardados en: \
    results/enrichment_${CELL_TYPE}_promoters_stats.csv, \
    results/enrichment_${CELL_TYPE}_PIRs_stats.csv, \
    results/enrichment_${CELL_TYPE}_both_enriched.csv"

# Finalizar
echo -e "\n╔════════════════════════════════════════╗"
echo      "║       PIPELINE COMPLETADO              ║"
echo      "╚════════════════════════════════════════╝"
