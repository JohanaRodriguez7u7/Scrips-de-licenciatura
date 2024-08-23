setwd("C:/Users/Johana Banana/Documents/Headmap y Go/REPETICION/Edadescompletas/Negativo")

rm(list=ls())
library("crayon") 

nameout <- paste("All_Groups")  # put name of file TodosM4top20
filterchar <- toString(paste("Negativo_"))  # Put filter and BP characteristics for name
outfilename <- toString(paste("GO_MLSP", nameout, "_", filterchar, sep=""))

# List all files that match the pattern and contain "LM_TMM" in their names
archivos <- list.files(pattern = "GO_")
archivos <- archivos[grep("Negativo_", archivos)]

# Re-order the species into the following order which reflects their increasing divergence level from dmel: dsec, dsim, dere, dyak, dana, dpse, dper, dwil, dmoj, dvir, dgri.
# archivos <- archivos[c(7,8,2,11,1,6,5,10,4,9,3)] # Drosophila 
# archivos <- archivos[c(3,1,4,2)] # Primate species
noms <- sapply(strsplit(archivos, "_"), "[", 2) # BP 4, # MF and CC 5. La parte del nombre cambiante/significativa

tb1 <- read.table(archivos[1], sep=",", header=TRUE)
ps <- matrix(NA, dim(tb1)[1], length(archivos), dimnames=list(tb1$GO.Term, noms))
zs <- matrix(NA, dim(tb1)[1], length(archivos), dimnames=list(tb1$GO.Term, noms))
obs <- matrix(NA, dim(tb1)[1], length(archivos), dimnames=list(tb1$GO.Term, noms))
exp_fam <- matrix(NA, dim(tb1)[1], length(archivos), dimnames=list(tb1$GO.Term, noms))

for (i in c(1:length(archivos))) {
  tb1 <- read.table(archivos[i], sep=",", header=TRUE)
  # tb1 <- tb1[order(tb1$GOID),]
  ps[,i] <- tb1$FDR
  zs[,i] <- tb1$Z.score
  obs[,i] <- tb1$Gene.Families.in.the.GO
  exp_fam[,i] <- tb1$Expected.families.per.GO # Añadir el filtro de familias esperadas
}

GOn <- apply(zs, 2, function(x) length(which(x > 0))) # valores de zscore positivos (por lo tanto solo categorías sobre representadas), tanto significativos como no

GOn2 <- apply(ps, 1, function(x) length(which(x < 0.01)))
ind <- which(GOn2 == 0)

newnames <- paste(colnames(ps), GOn)
colnames(ps) <- newnames
ps[which(ps > 0.01)] <- NA
ps[which(zs < 0)] <- NA
ps[which(obs < 5)] <- NA # filter that more than 1 Gene Family is observed
ps[which(exp_fam < 1)] <- NA # filtro para familias esperadas por conjunto GO
# Filtrar GO que se expresan en más de una lista
#ps <- ps[apply(ps, 1, function(x) sum(!is.na(x)) > 1), ]

ps <- ps[-ind, ]
ps <- ps[order(-ps[,2], -ps[,1]), ] # para ordenar (ascendente) primero con el p-value de la 1a columna, luego por el de la 2a y finalmente el de la 3a
ps <- ps[apply(ps, 1, function(x) any(!is.na(x))),] # TRUE para las filas que no son todas NA

#***
## Save "ps" and manually change the order as desired for the graphic
## First invert the order to see it as the final output graphic
ps <- ps[seq(dim(ps)[1], 1),]
tablename <- toString(paste(outfilename, "_Table", sep=""))
write.csv(ps, paste(tablename, ".csv", sep=""))

### CHANGE the order as desired for the graphic
ps <- read.csv(toString(paste(tablename, ".csv", sep="")))

# psextra <- read.csv(toString(paste("EiconM4TOP20WOut", ".csv", sep=""))) to add the top20 on a separate file
# ps <- merge(ps, psextra, all=TRUE)

ps <- read.csv("GO_MLSPAll_Groups_Negativo__Table.csv")

ps <- ps[seq(dim(ps)[1], 1),] # re-invert to the original
ps2 <- ps[,-1]
rownames(ps2) <- ps[,1]
ps <- ps2

# Eliminar el prefijo "GO_MLSP" de los nombres de las columnas
#colnames(ps) <- archivos
#colnames(ps) <- gsub("GO_MLSP", "", colnames(ps))

#ps <- ps[apply(ps, 1, function(x) sum(!is.na(x)) > 1), ]

# Actual Heatmap plot #
library(reshape)
library(lattice)
library(RColorBrewer)
colores <- c(rev(brewer.pal(9,"YlOrRd")[2:9]))
# colores <- c(brewer.pal(11,"RdYlBu")[1:11])

# colores2 <- c(rev(brewer.pal(9,"Blues")[2:5]))
paleta <- colorRampPalette(colores, space = "Lab")
# paleta2 <- colorRampPalette(colores, space = "Lab")
colores <- paleta(50)
interv <- seq(0, 0.01, by=.001)

pdfname <- toString(paste("heatmap_paper.pdf", sep=""))
pdf(pdfname, paper="a4r", width=8, height=8) # BP 11,6 CC 8,8
levelplot(t(ps), main="Positivo", xlab="", ylab="", col.regions=colores, cuts=100, at=interv, scales=list(x=list(rot=45, cex=0.7), y=list(rot=0, cex=0.7), draw=TRUE), colorkey=list(TRUE, space="right", contour=FALSE, pretty=TRUE), par.settings=list(axis.line=list(col=0)),
          panel=function(...) {
            panel.fill(col="gray96")
            panel.levelplot(...)
          })
dev.off()

png("heatmap_for_paper.png", units="px", width=5760, height=5760, res=600)
levelplot(t(ps), main="Positivo", xlab="", ylab="", col.regions=colores, cuts=100, at=interv, scales=list(x=list(rot=45, cex=1.5), y=list(rot=0, cex=1.5), draw=TRUE), colorkey=list(labels=list(cex=1.25), TRUE, space="right", contour=FALSE, pretty=TRUE), par.settings=list(axis.line=list(col=0)),
          panel=function(...) {
            panel.fill(col="gray96")
            panel.levelplot(...)
          })
dev.off()

library(viridis)
colvir <- viridis(100)

pdfname <- toString(paste("heatmap_paper_choosed_groups_LAST.pdf", sep=""))
pdf(pdfname, paper="a4r", width=8, height=8) # BP 11,6 CC 8,8
levelplot(t(ps), main="Positivo", xlab="", ylab="", col.regions=plasma(150), cuts=100, at=interv, scales=list(x=list(rot=45, cex=0.7), y=list(rot=0, cex=0.7), draw=TRUE), colorkey=list(TRUE, space="right", contour=FALSE, pretty=TRUE), par.settings=list(axis.line=list(col=0)),
          panel=function(...) {
            panel.fill(col="gray96")
            panel.levelplot(...)
          })
dev.off()

png("heatmap_for_paper_choosed_groups_LAST.png", units="px", width=3000, height=3000, res=350)
levelplot(t(ps), main="Positivo", xlab="", ylab="", col.regions=plasma(150), cuts=100, at=interv, scales=list(x=list(rot=45, cex=1), y=list(rot=0, cex=1), draw=TRUE), colorkey=list(TRUE, space="right", contour=FALSE, pretty=TRUE), par.settings=list(axis.line=list(col=0)),
          panel=function(...) {
            panel.fill(col="gray96")
            panel.levelplot(...)
          })
dev.off()
