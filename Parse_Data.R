#!/usr/bin/env Rscript
# --------------------------------------------------------
# Extraer datos PCHiC con barra de progreso
# Uso: Rscript Parse_Data.R <data_table> <cell_type>
# --------------------------------------------------------

# 1) Extraer argumentos del comando
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript Parse_Data.R <data_table> <cell_type>")
}
data_table <- args[1]
cell_type  <- args[2]

# Definir pasos y barra de progreso
steps <- c("Leyendo tabla", "Filtrando datos", "Obteniendo Promotores", "Obteniendo PIRs")
pb    <- txtProgressBar(min = 0, max = length(steps), style = 3)
step  <- 0

# 2) Leer la tabla
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(steps), steps[step]))
setTxtProgressBar(pb, step)
raw_data <- read.table(data_table, header = TRUE)

# Validar que la célula existe
if (!(cell_type %in% colnames(raw_data))) {
  close(pb)
  stop(paste("La célula", cell_type, "no se encuentra en la tabla."))
}

# 3) Filtrar filas con score > 5
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(steps), steps[step]))
setTxtProgressBar(pb, step)
filtered_data <- raw_data[ raw_data[[cell_type]] > 5, ]

# 4) Obtener Promotores
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(steps), steps[step]))
setTxtProgressBar(pb, step)
promotores <- filtered_data[, c("baitChr", "baitStart", "baitEnd", "baitID")]
colnames(promotores) <- c("chr", "start", "end", "ID")
promotores <- promotores[!duplicated(promotores), ]
promotores$chr <- paste0("chr", promotores$chr)
write.table(promotores,
            file = paste0("data/promotores_", cell_type, ".bed"),
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# 5) Obtener PIRs
step <- step + 1
cat(sprintf("[%d/%d] %s...\n", step, length(steps), steps[step]))
setTxtProgressBar(pb, step)
PIR <- filtered_data[, c("oeChr", "oeStart", "oeEnd", "oeID")]
colnames(PIR) <- c("chr", "start", "end", "ID")
PIR <- PIR[!duplicated(PIR), ]
PIR$chr <- paste0("chr", PIR$chr)
write.table(PIR,
            file = paste0("data/PIR_", cell_type, ".bed"),
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# Cerrar barra de progreso
close(pb)
cat(sprintf("Proceso completado. Archivos en data/: promotores_%s.bed, PIR_%s.bed\n", cell_type, cell_type))
