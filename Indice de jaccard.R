# Instalar y cargar el paquete sets si no está ya instalado
if (!require("sets")) {
  install.packages("sets")
}
library(sets)

# Establecer el directorio de trabajo
setwd("C:/Users/Johana Banana/Documents/Repeticiones/JACCARD")

# Leer las listas de genes
LM_TMM_fantom <- read.csv("LM_TMM_Edadescom_FANTOM.csv")
LM_Quartile_fantom <- read.csv("LM_Cuartiles_Edadescom_FANTOM.csv")
PGLS_TMM_fantom <- read.csv("PGLS_TMM_Edadescom_FANTOM.csv")
PGLS_Quartile_fantom <- read.csv("PGLS_Cuartiles_Edadescom_FANTOM.csv")
LM_TMM_Suetal <- read.csv("LM_TMM_Edadescom_SUETAL.csv")
LM_Quartile_suetal <- read.csv("LM_Cuartiles_Edadescom_SUETAL.csv")
PGLS_TMM_suetal <- read.csv("PGLS_TMM_Edadescom_SUETAL.csv")
PGLS_Quartile_suetal <- read.csv("PGLS_Cuartiles_Edadescom_SUETAL.csv")
Lista_PGLS_fantom<- read.csv("ListacompartidaPGLS.csv")

# Extraer los genes de cada lista
genes_LM_TMM_fantom <- LM_TMM_fantom$X
genes_LM_Quartile_fantom <- LM_Quartile_fantom$GeneID
genes_PGLS_TMM_fantom <- PGLS_TMM_fantom$X
genes_PGLS_Quartile_fantom <- PGLS_Quartile_fantom$X
genes_LM_TMM_Suetal <- LM_TMM_Suetal$X
genes_LM_Quartile_suetal <- LM_Quartile_suetal$X
genes_PGLS_TMM_suetal <- PGLS_TMM_suetal$X
genes_PGLS_Quartile_suetal <- PGLS_Quartile_suetal$X
genes_lista_PGLS_fantom <-Lista_PGLS_fantom$x
# Crear una lista con todas las listas de genes
gene_lists <- list(
  LM_TMM_fantom = genes_LM_TMM_fantom,
  LM_Quartile_fantom = genes_LM_Quartile_fantom,
  PGLS_TMM_fantom = genes_PGLS_TMM_fantom,
  PGLS_Quartile_fantom = genes_PGLS_Quartile_fantom,
  LM_TMM_Suetal = genes_LM_TMM_Suetal,
  LM_Quartile_suetal = genes_LM_Quartile_suetal,
  PGLS_TMM_suetal = genes_PGLS_TMM_suetal,
  PGLS_Quartile_suetal = genes_PGLS_Quartile_suetal,
  Lista_PGLS_fantom = genes_lista_PGLS_fantom
)

# Definir una función para calcular el índice de Jaccard entre dos listas de genes
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Calcular el índice de Jaccard para cada par de listas de genes
n <- length(gene_lists)
jaccard_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(gene_lists)

for (i in 1:n) {
  for (j in 1:n) {
    jaccard_matrix[i, j] <- jaccard_index(gene_lists[[i]], gene_lists[[j]])
  }
}

# Imprimir la matriz de índices de Jaccard
print("Matriz de índices de Jaccard:")
print(jaccard_matrix)

# Guardar la matriz de índices de Jaccard en un archivo CSV
write.csv(jaccard_matrix, "jaccard_matrix.csv")




# En tu caso, sustituye jaccard_matrix por tu matriz real de índices de Jaccard
# Imprimir la matriz de índices de Jaccard
print("Matriz de índices de Jaccard:")
print(jaccard_matrix)

# Guardar la matriz de índices de Jaccard en un archivo CSV
write.csv(jaccard_matrix, "jaccard_matrix.csv")

# Crear el mapa de calor interactivo
heatmaply(jaccard_matrix,
          main = "Mapa de Calor de Índices de Jaccard entre Listas de Genes",
          fontsize = 12,
          cellnote = round(jaccard_matrix, 2),
          notecol = "black",
          margins = c(80, 100))
