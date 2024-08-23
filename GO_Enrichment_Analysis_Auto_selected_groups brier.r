#Script to measure if a list of Families is particularly enriched in genes pertaining to a Gene Onthology Term

####### ONLY VARIABLES THAT NEEDS TO BE CHANGED :  ######################
#### infile , GOanott,  phenotype, filenameout 
#### Those variables reduced to one so we also reduce the redundancy...
#########################################################################

###########################################
###########################################
###########################################
###########################################  POSITIVES POSITIVES POSITIVES POSITIVES POSITIVES POSITIVES
###########################################
###########################################
###########################################

rm(list=ls())
setwd("C:/Users/Estudiantes/Documents/Brier/Modificaciones/GO")
library(parallel)

cl <- makeCluster(detectCores()-1) 


###########*************
filename1<-toString(paste("PGLS.Estimated_longevity_days_to_normal.18.prueba."))  #Change PGLS results **********
filenameout <- gsub('^(?:[^ ]* ){1}',' ',filename1)#to exclude everything before the 1rd space and keep the name's part of interest to name the output files
filename2<-toString(paste(filename1, ".csv", sep=""))

infile<-read.csv(filename2,stringsAsFactors=FALSE)
head(infile)
#infile <-  subset(infile, select = -c(X.1)) #para quitar la columna extra de rownames numerados en la Perm

phenotype<-toString(paste("MLSP+Ei"))#Put Phenotype of interest###########*************


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

rm(list=setdiff(ls(), c("infile", "GOanott", "phenotype", "filenameout")))
#setwd("/Users/Administrator/Documents/PhD/Proyecto/GliaN_Paper_Obs/Pearson/PIC")
library(parallel)

cl <- makeCluster(detectCores()-1) 


filterbp<-toString(paste("benjamini BP50"))#Put filter and BP characteristics for name###########*************

outfile <- toString(paste("GO_",  phenotype, "_", filterbp, sep=""))
outfile1 <- toString(paste(outfile, ".csv", sep=""))

numreps<-1000



# mlsp_assoc <- read.csv("../Group5B_46GFs_MLSP.csv", stringsAsFactors=FALSE)
# head(mlsp_assoc)


mlsp_assoc <- infile[infile$p.adjusted.GFS.vs.Estimated_longevity_days_to_normal < 0.05 & infile$R.t.value.Estimated_longevity_days_to_normal > 0 ,]
head(mlsp_assoc)
dim(mlsp_assoc)

Familyset <- as.data.frame(infile[infile$X %in% mlsp_assoc$X, 1])

head(Familyset)
dim(Familyset)


names(Familyset)<- c("FamilyID")
dim(Familyset)


#Population of all gene families studied (Background)
dat<-as.data.frame(infile[,1])
names(dat)<- c("FamilyID")

GOanott<-read.csv("GO_Terms_per_family_Ens99BP_50.csv",header=TRUE,stringsAsFactors=FALSE) #Gene onthology annotations per family###########*************

GOsBPn<-unique(GOanott[,-1])
GOsBP<-GOsBPn[,1]

GOanott2<-matrix(NA,dim(GOanott)[1],2,dimnames = list(c(1:dim(GOanott)[1]),c("GO.ID", "FamilyID")))
GOanott2[,1]<-GOanott[,2]
GOanott2[,2]<-GOanott[,1]

counts<-vector("numeric")
annotation<-vector("numeric")
join<-merge(Familyset,GOanott2,by="FamilyID")

w <- dim(Familyset)[1]

#countssmpl<-matrix(NA,length(GOsBP),numreps)
zscore<-matrix(NA,length(GOsBP),1,dimnames=list(GOsBP,"Z-score"))
#annotationsmpl<-matrix(NA,w,numreps)
countssmpl2<-matrix(NA,length(GOsBP),numreps)
zscore2<-matrix(NA,length(GOsBP),1,dimnames=list(GOsBP,"Z-score"))


counts<-sapply(GOsBP, function(x) #Counts in Familyset per GO
  sum(join[,2]==x))

annotation<-apply(Familyset,1, function(x) #Annotation density of each family in Familyset
  sum(join[,1]==x))

counts2<-counts/mean(annotation) #Counts in Familyset per GO corrected against annotation density


clusterExport(cl, c("dat", "GOanott2", "GOsBP", "w"))
run3 <- function(c) {
  smpl <- as.data.frame(sample(dat[,1],w),stringsAsFactors=FALSE)
  names(smpl)<- c("FamilyID")
  joinsmpl<-merge(smpl,GOanott2,by="FamilyID")
  
  countssmpl<-sapply(GOsBP, function(x) #Cuentas por GO Term de genes en la smpl
    sum(joinsmpl[,2]==x))
  
  annotationsmpl<-apply(smpl,1, function(x) #Densidad de anotacion en reales
    sum(joinsmpl[,1]==x))
  
  
  return(list(countssmpl=countssmpl, annotationsmpl=annotationsmpl))
}

MCsmpl<-parSapply(cl,rep(1,numreps),run3)
stopCluster(cl)

countssmpl<-do.call(cbind,MCsmpl[1,])
annotationsmpl<-do.call(cbind,MCsmpl[2,])


means<-rowMeans(countssmpl)
meansannotation<-colMeans(annotationsmpl)
for(p in c(1:numreps))
  countssmpl2[,p]<-countssmpl[,p]/meansannotation[p]
means2<-rowMeans(countssmpl2)
SEMs<-apply(countssmpl,1,sd)
SEMs2<-apply(countssmpl2,1,sd)
zscore[,1]<-(counts-means)/SEMs
zscore2[,1]<-(counts2-means2)/SEMs2
pval<-pnorm(-abs(zscore))
adjpval<-p.adjust(pval,"fdr")
pval2<-pnorm(-abs(zscore2))
adjpval2<-p.adjust(pval2,"fdr")
numpval<-rowSums(countssmpl>=counts)/numreps
numpval2<-rowSums(countssmpl2>=counts2)/numreps
adjnumpval<-p.adjust(numpval,"fdr")
adjnumpval2<-p.adjust(numpval2,"fdr")

out1<-as.data.frame(cbind(GOsBPn[,1],GOsBPn[,2],counts,means,SEMs,zscore,pval, adjpval, numpval, adjnumpval),stringsAsFactors=FALSE)
#out2<-as.data.frame(cbind(GOsBPn[,1],GOsBPn[,2],counts2,means2,SEMs2,zscore2,pval2, adjpval2, numpval, adjnumpval),stringsAsFactors=FALSE)


names(out1)<-c("GOID","GO Term","Gene Families in the GO","Expected families per GO","Standard Deviation of the samples","Z-score","One tail p-value","FDR", "Numeric p-value", "Numeric FDR")
#names(out2)<-c("GOID","GO Term","Gene Families in the GO","Expected families per GO","Standard Deviation of the samples","Z-score","One tail p-value","FDR", "Numeric p-value", "Numeric FDR")


#write.csv(out1, outfile, row.names=F)
write.csv(out1, outfile1, row.names=F)

alarm()


