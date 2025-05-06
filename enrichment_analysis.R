#!/usr/bin/env Rscript
# enrichment_analysis.R

# ---- Instalar paquetes si es necesario ----
required <- c('Biostrings','dplyr')
installed <- rownames(installed.packages())
for(p in required){ if(!p %in% installed) BiocManager::install(p) }
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

# --- Argumentos de línea de comandos ---
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 7) {
  stop("Uso: Rscript enrichment_analysis.R <prom_fasta> <prom_orig_out> <prom_perm_out> ",
       "<pir_fasta> <pir_orig_out> <pir_perm_out> <output_prefix>\n")
}
prom_fa    <- args[1]
prom_orig  <- args[2]
prom_perm  <- args[3]
pir_fa     <- args[4]
pir_orig   <- args[5]
pir_perm   <- args[6]
out_pref   <- args[7]

# --- Función para contar hits por motif ---
count_hits <- function(scan_file) {
  df <- read.table(scan_file, header=TRUE, sep="\t", comment.char=";")
  df %>% select(sequence_id, matrix_id) %>% distinct() %>%
    group_by(matrix_id) %>% summarise(hits = n()) %>% rename(motif = matrix_id)
}

# --- Contar número de secuencias en cada FASTA ---
n_prom <- length(readDNAStringSet(prom_fa))
n_pir  <- length(readDNAStringSet(pir_fa))

# --- Cargar conteos ---
prom_o <- count_hits(prom_orig)
prom_p <- count_hits(prom_perm)
pir_o  <- count_hits(pir_orig)
pir_p  <- count_hits(pir_perm)

# --- Unir con todos los motifs posibles ---
all_motifs <- unique(c(prom_o$motif, prom_p$motif, pir_o$motif, pir_p$motif))
df_prom <- tibble(motif = all_motifs) %>%
  left_join(prom_o, by="motif") %>% replace_na(list(hits=0)) %>% rename(h1=hits) %>%
  left_join(prom_p, by="motif") %>% replace_na(list(hits=0)) %>% rename(h0=hits) %>%
  rowwise() %>% mutate(
    pvalue = fisher.test(
      matrix(c(h1, h0, n_prom - h1, n_prom - h0), nrow=2),
      alternative = "greater"
    )$p.value
  ) %>% ungroup() %>% mutate(padj = p.adjust(pvalue, method="BH")) %>% arrange(padj)

df_pir <- tibble(motif = all_motifs) %>%
  left_join(pir_o, by="motif") %>% replace_na(list(hits=0)) %>% rename(h1=hits) %>%
  left_join(pir_p, by="motif") %>% replace_na(list(hits=0)) %>% rename(h0=hits) %>%
  rowwise() %>% mutate(
    pvalue = fisher.test(
      matrix(c(h1, h0, n_pir - h1, n_pir - h0), nrow=2),
      alternative = "greater"
    )$p.value
  ) %>% ungroup() %>% mutate(padj = p.adjust(pvalue, method="BH")) %>% arrange(padj)

# --- TFs enriquecidos (padj < 0.05) en ambos conjuntos ---
enr_prom <- df_prom %>% filter(padj < 0.05) %>% pull(motif)
enr_pir  <- df_pir  %>% filter(padj < 0.05) %>% pull(motif)
common   <- intersect(enr_prom, enr_pir)

# --- Guardar resultados ---
write.csv(df_prom, paste0(out_pref, "_promoters_stats.csv"), row.names=FALSE)
write.csv(df_pir,  paste0(out_pref, "_PIRs_stats.csv"),       row.names=FALSE)
write.csv(data.frame(motif=common), paste0(out_pref, "_both_enriched.csv"), row.names=FALSE)

# --- Mensajes finales ---
cat("\n--- TFs enriquecidos en PROMOTORES (padj<0.05):\n")
print(enr_prom)
cat("\n--- TFs enriquecidos en PIRs (padj<0.05):\n")
print(enr_pir)
cat("\n--- TFs comunes a ambos: \n")
print(common)
