# Cargar paquetes necesarios
install.packages("sets")
library(sets)

# Establecer el directorio de trabajo
setwd("C:/Users/Johana Banana/Documents/Repeticiones/intercecciones")

# Lectura de listas
PGLS_TMM_fantom <- read.csv("LM_Cuartiles_Edadescom_SUETAL.csv")
PGLS_Quartile_fantom <- read.csv("LM_TMM_Edadescom_SUETAL.csv")

# Extraer las listas de genes
genes_PGLS_TMM_fantom <- PGLS_TMM_fantom$X
genes_PGLS_Quartile_fantom <- PGLS_Quartile_fantom$X

# Unir todas las listas de genes en un único vector
all_genes <- c(genes_PGLS_TMM_fantom, genes_PGLS_Quartile_fantom)

# Contar la frecuencia de cada gen
gene_frequencies <- table(all_genes)
# Encontrar genes que se comparten entre las listas
shared_genes <- intersect(genes_PGLS_TMM_fantom, genes_PGLS_Quartile_fantom)

# Guardar los genes compartidos en un archivo CSV
write.csv(shared_genes, "genescompartidos20.csv", row.names = FALSE)
