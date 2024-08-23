# Establecer el directorio de trabajo
setwd("C:/Users/Johana Banana/Documents/fantom/Modifiaciones/Listas")
# Leer el archivo CSV con los datos
datos <- read.csv("UNIQUECOLs.rename.phatom data human.tissue.hCAGE.hg19.tpm.refgene.osc.csv", check.names = FALSE)

# Verificar las primeras filas para asegurarse de que la columna 'Gen' está presente
head(datos)

# Definir la función para promediar filas con el mismo nombre
promediar_filas <- function(datos) {
  datos_promedio <- aggregate(. ~ ENSG, data = datos, FUN = mean, na.rm = TRUE)
  return(datos_promedio)
}

# Aplicar la función para obtener los promedios de las filas con el mismo nombre
datos_promediados <- promediar_filas(datos)

# Guardar los datos promediados en un nuevo archivo CSV
write.csv(datos_promediados, "promedsensembl.csv", row.names = FALSE)
