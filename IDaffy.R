


install.packages("genekitr")
install.packages("clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

# Cargar el paquete biomaRt
library(genekitr)
library(clusterProfiler)



# Establecer el directorio de trabajo
#setwd("C:/Users/Johana Banana/Documents/Scrips")
setwd("C:/Users/Johana Banana/Documents/Brier/Modificados/Jacknife")

# Leer el archivo CSV con identificadores de miARN
mirna <- read.csv("GenessignificativosArticuloAtahuapa.csv")

# Lista de IDs de miARN (asegúrate de que la columna tiene el nombre correcto)
miarn_ids <- mirna$SYM

gene_mapping <- transId(id = miarn_ids, transTo = "ens", keepNA = TRUE)

# Obtener los identificadores de ENSG en el mismo orden que los originales
miarn_ensg <- gene_mapping$ensembl

# Crear un nuevo data frame con los identificadores de miARN originales y sus correspondientes identificadores de ENSG
resultado <- data.frame(IDENTIFIER = miarn_ids, ENSG = miarn_ensg)

# Guardar el resultado en un archivo CSV
write.csv(gene_mapping, file = "gene_mapping.csv", row.names = FALSE)
