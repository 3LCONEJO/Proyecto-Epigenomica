#!/usr/bin/env Rscript
# --------------------------------------------------------
# Uso: Rscript Prueba_diferenciacion.R <cell_type> <histone_mark>
# --------------------------------------------------------

# ---- Instalar paquetes si es necesario ----
required <- c('ggplot2')
installed <- rownames(installed.packages())
for(p in required){ if(!p %in% installed) install.packages(p, repos = 'https://cloud.r-project.org') }

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript Prueba_diferenciacion.R <cell_type> <histone_mark>")
}
cell_type <- args[1]
histone   <- args[2]

# Directorio de salida
outdir <- "data"

tasks <- c("Leyendo archivos de señal", 
           "Preparando dataframe", 
           "Generando gráfica", 
           "Realizando test estadístico")
pb <- txtProgressBar(min = 0, max = length(tasks), style = 3)
step <- 0

# 1) Leer archivos de señal
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(tasks), tasks[step]))
setTxtProgressBar(pb, step)
prom_file <- file.path(outdir, paste0("prom_", cell_type, ".tab"))
pir_file  <- file.path(outdir, paste0("pir_",       cell_type, ".tab"))

if (!file.exists(prom_file) || !file.exists(pir_file)) {
  close(pb)
  stop("No se encontraron los archivos:\n  ", prom_file, "\n  ", pir_file)
}

promotores  <- read.table(prom_file, header = FALSE)
pir         <- read.table(pir_file,  header = FALSE)
prom_signal <- promotores$V5
pir_signal  <- pir$V5

# 2) Preparar dataframe
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(tasks), tasks[step]))
setTxtProgressBar(pb, step)
signal <- c(prom_signal, pir_signal)
region <- factor(
  c(rep("Promotores", length(prom_signal)), rep("PIR", length(pir_signal))),
  levels = c("Promotores", "PIR")
)
df <- data.frame(signal = signal, region = region)

# 3) Generar gráfica
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(tasks), tasks[step]))
setTxtProgressBar(pb, step)

plot_file <- file.path(outdir,paste0("signal_", cell_type, "_", histone, ".png"))
png(filename = plot_file, width = 800, height = 600)
suppressPackageStartupMessages(library(ggplot2))

ggplot(df, aes(x = region, y = signal, fill = region)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = paste0("Señal ", histone), y = "Signal") +
  scale_y_log10()
dev.off()

# 4) Realizar test de Wilcoxon
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(tasks), tasks[step]))
setTxtProgressBar(pb, step)

test_res <- wilcox.test(prom_signal, pir_signal)
pval <- test_res$p.value
sig <- ifelse(pval < 0.05, "estadísticamente significativa", "no estadísticamente significativa")

close(pb)

cat("=== Resultados del test de Wilcoxon ===\n")
print(test_res)
cat(sprintf("\nLa diferencia entre Promotores y PIR es %s (p = %.3g).\n", sig, pval))
cat(sprintf("\nGráfica guardada en:\n  %s\n", plot_file))
