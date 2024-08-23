## for pre lm analyses
library("logr")
library("dplyr")
library("ggpubr")

## for pgls
library("ape")
library("nlme")
library("parallel")

## for the plot 
library("reshape")
library("ggplot2")

print("change always the ulimit of the computer $ulimit -s 21000 !!!!!!!!!!!!")
rm(list=ls())

#*************   Change parameters!!!!!!!   ********###
trait<- "Estimated_longevity_days_to_normal"
Correcting<- "benjamini"
#*************   ************************** ********###

tmp <- file.path(getwd(), paste("GenefamilySize.",trait,".associated.genes.", sep = "")) #,Correcting,".correction.log"
lf <- log_open(tmp, traceback = F)

setwd("C:/Users/Johana Banana/Documents/fantom/Modifiaciones/LM")
#setwd("/home/jvar2023/Bases")

##### First we need to do the genome completion filtering. ####

## upload all the orthologs (gene families) per spp list. 

#*************   Change parameters!!!!!!!   ********###
list_orthogroups_counts<-read.csv("Normalizacion_Fantom_TMMVERDADERO.csv")
traits<-read.csv("Long_days_Fantom_exacto.csv")
#*************   ************************** ********###

# Check if the phenotype interested is normally distributed. If not get the LOG values
#traits$clog.longevity<-log10(traits$Maximum.longevity..yrs.)
traits.filtered <-subset(traits,traits$fantom.data%in% colnames(list_orthogroups_counts))



#Remove raits without phenotypic or data of interest
traits.filtered <- traits[!is.na(traits[[trait]]),]
traits.filtered <- traits.filtered[!is.na(traits$fantom.data),]


# Obtener los nombres de las columnas que están en ambas bases de datos

#*gene_numbers<<-list_orthogroups_counts...
#traits.filtered<<-traits[!is.na(traits[,phenotypes[i]]),]

#list_orthogroups_counts <- list_orthogroups_counts[1:8,] #Temporal, solo deben haber nombres unicos (sin repetirse nombres de filas). 



row.names(list_orthogroups_counts)<-list_orthogroups_counts$ENSG#Assign gene names as row names

gene.numbers.filtered<<-list_orthogroups_counts[colnames(list_orthogroups_counts) %in% traits.filtered$fantom.data]

traits.filtered$fantom.data
colnames(gene.numbers.filtered)

#*************   IT IS IMPORTANT TO HAVE THE SAME SPP ORDER BETWEEN TRAITS AND GENE NUMBERSS!!!!!!!   ********###
gene.numbers.filtered <- gene.numbers.filtered[,traits.filtered[,4]]
#*************   ************************** ********###

traits.filtered$fantom.data
colnames(gene.numbers.filtered)


# total counts per tissue
sumgenes<-rowSums(t(gene.numbers.filtered))
sumgenes

#?summary.lm

###############################################
# Open function for more detail of the process#
###############################################

################
##### LM #####
### function ###
################


LM.GF.size<- function(traits, pheno.col.nums, spp.col.num, gene.numbers, No.variables, where.Save.it){
  phenotypes<<- colnames(x = traits)[c(spp.col.num, pheno.col.nums)]
  log_print("#########################", hide_notes = T)
  log_print("##### LM LM LM ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  if (No.variables == 1){
    ###################################
    log_print("#### LM ONE Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,1,FUN = function(x){
      
      traits2<<-cbind(traits.filtered[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      
      
      ## creating the formula for the model for 1 variable
      formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
      lmModel <- try(lm(formilin, data = traits2)) #method = "ML",
      # pglsModel<-try(gls(model = formilin, correlation = corBrownian(phy = tree2, form = formilin2), method = "ML", data = traits2))
      
      if (inherits(lmModel, "try-error")) 
        return(c(NA,NA))
      else
        return(c(summary(lmModel)$adj.r.squared, summary(lmModel)$r.squared, summary(lmModel)$coefficients[2,4], summary(lmModel)$coefficients[2,"t value"]))
    })))
    #colnames(out)<<-c("adj_r_squared", "r_squared", "p-value",paste
    colnames(out)<<-c(paste("adj_r_squared", phenotypes[2], sep=" "),
                      paste("r_squared", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "), paste("Coef tval", phenotypes[2]))
    
  } 
  else if (No.variables  > 2) {
    ###################################
    log_print("#### LM TWO Variable ######", hide_notes = T)
    ###################################
    out<<-as.data.frame(t(apply(X = gene.numbers,MARGIN = 1,FUN = function(x){
      pheno1<<-phenotypes[3]
      pheno2<<-phenotypes[2]
      traits2<<-cbind(traits[,c(phenotypes[1], pheno2, pheno1)],as.numeric(x)) #add phenotypes as needed
      # traits2<<-cbind(traits.filtered[,c(phenotypes[1], pheno2, pheno1)],t(gene.numbers.filtered[,1])) #add phenotypes as needed
      # traits2<<-cbind(traits[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
      names(traits2)[dim(traits2)[2]]<<-"GFS"
      #next bit removes Genus.Species with NAs in events from tree, only necessary for losses, deletion, not for GFS
      
      ## creating the formula for the model for 2 variables
      # formilin<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
      formilin<<-as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
      #log_print((paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ ")))
      # formilin2<-as.formula(paste("~",colnames(traits)[spp.col.num], sep = ""))
      formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
      #log_print((paste("~",colnames(traits)[spp.col.num], sep = "")))
      # pglsModel<-try(gls(model = formilin, correlation = corBrownian(phy = tree, form = formilin2), method = "ML", data = traits2))
      lmModel <- try(lm(model = formilin, data = traits2))  # method = "ML",
      
      if (inherits(lmModel, "try-error"))
        return(c(NA,NA))
      else
        return(c(summary(lmModel)$adj.r.squared, summary(lmModel)$r.squared, summary(lmModel)$coefficients[2,4], summary(lmModel)$coefficients[2,"t value"]))
    }))) 
    colnames(out)<<-c(paste("adj_r_squared", phenotypes[2], sep=" "),
                      paste("r_squared", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "), paste("Coef tval", phenotypes[2]))
    
    
  } ######CHECK if the p-value is the value we extract here summary(lmModel)$coefficients[2,4], for ONE variable is OK#####################
  else if (No.variables > 2){
    ###################################
    ## LM more than two Variables ###
    ###################################
    log_print("At the moment the function only works with upto 2 variables, sorry :P", hide_notes = T)
  }
  for (i in 1:No.variables) {
    
    ###########################
    ### R^2 calculation ###
    ###########################
    log_print("Calc. Rs", hide_notes = T)
    sps<-nrow(traits2)
    DeFr <- sps - (1 + No.variables)
    Coef.tval1<-toString(paste("Coef tval", phenotypes[i+1], sep=" "))
    R1 <- toString(paste("R t value", phenotypes[i+1], sep=" "))
    R_t_value1 <- out[[Coef.tval1]] / (sqrt(((out[[Coef.tval1]])^2) + DeFr)) # tvalue / square root of(tvalue^2 + DF))
    out[[R1]] <<- R_t_value1
    ### benjamini correction
    log_print("Calc. benjamini correction", hide_notes = T)
    out[[paste("p.adjusted GFS vs.", phenotypes[i+1], sep = "")]] <<- p.adjust(out[,paste("p_value", phenotypes[i+1], sep=" ")], method = "fdr")
  }
  
  out <- out[complete.cases(out),]
  nameout<- paste("LM", paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits2)[1], "Spp", dim(out)[1], "GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  log_print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  log_print(paste("object created:", paste( "out", sep = "")), hide_notes = T)
  write.csv(out, paste(where.Save.it, nameout, sep = "/"))
}

### LM function
 #dirSave<-"/Volumes/PortableSSD/Mammalian_orthologs_project/PLGS"
#*************   Change parameters!!!!!!!   ********###
#dirSave<-"/home/jvar2023/Bases"
dirSave<- "C:/Users/Johana Banana/Documents/fantom/Modifiaciones/LM"



#*************   ************************** ********###

log_print(paste("You have ",dim(gene.numbers.filtered)[1], " gene families to input the LM", sep = ""), hide_notes = T)


#*************   Change parameters!!!!!!! #CHECK that column numbers are the ones we want (in pheno.col.nums and  spp.col.num)
#*  ********###
LM.GF.size(traits = traits, pheno.col.nums = c(3), No.variables =1, 
           spp.col.num = 4, gene.numbers = gene.numbers.filtered, 
           where.Save.it = dirSave)

#*************   ************************** ********###



##############https://www.statology.org/interpret-prt-regression-output-r/


####At the very bottom of the output we can see the overall p-value 
##for the regression model.
####If we would like to only extract this p-value from the model, 
#we can define a custom function to do so:
#########----
#define function to extract overall p-value of model
#overall_p <- function(lmModel) {
#  f <- summary(lmModel)$fstatistic
#  p <- pf(f[1],f[2],f[3],lower.tail=F)
#  attributes(p) <- NULL
#  return(p)
#}

#extract overall p-value of model
#overall_p(lmModel)

#check
#summary(lmModel)

#lo mismo que esto:
##Extract p-value, located in $coefficients under Pr(>|t|) column
#summarygls$coefficients[2,4] ##Pr(>|t|) column represents the p-value associated with the value in the t value column.

#*************   ************************** ********###

#closeAllConnections()

colnames(out)

#log_close()
writeLines(readLines(lf))

###################################
#### PERMS PERMS PERMS PERMS ######
###################################
LM.GF.size.perms<-function(traits.cl, pheno.col.nums, cores, spp.col.num, gene.numbers.cl, where.Save.it, No.perms){
  log_print("###########################", hide_notes = T)
  log_print("##### Perms perms LM ####", hide_notes = T)
  log_print("###########################", hide_notes = T)
  log_print("calc. cores", hide_notes = T)
  nucleos<-detectCores()
  nucleos<-nucleos + cores - nucleos
  cl <- makeCluster(nucleos) 
  log_print(paste(nucleos + cores - nucleos," cores detected and ready to roll"), hide_notes = T)
  
  gene.numbers.cl<- gene.numbers.cl
  traits.cl<- traits.cl
  #tree.cl<- tree.cl
  phenotypes<<- colnames(x = traits.cl)[c(spp.col.num, pheno.col.nums)]
  traits.cl<<-traits.cl
  No.perms <- No.perms
  
  ##### this way  clusterExports stops using .GlobalEnv and uses the eviropnment inside the function. 
  clusterExport(cl = cl , varlist = c("gene.numbers.cl", "traits.cl", "lm", "phenotypes", "No.perms"), envir = environment())
  
  log_print("Starting perms LM", hide_notes = T)
  
  if (length(pheno.col.nums) == 2) {
    perms<<- parLapply(cl, 1:No.perms, function(i,...){
      traits.cl<-traits.cl
      y<-sample(c(1:length(gene.numbers.cl)))
      gene_numbers2 <- gene.numbers.cl[,y]
      outperm1<-apply(gene_numbers2, MARGIN = 1,function(x){
        traits2<<-cbind(traits.cl[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
        #traits2<<-as.data.frame(cbind(traits.cl[,c("Spp.Name" , "Av.Body.mass" , "SSD")],x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c(cat(paste0(dQuote(phenotypes[1:length(phenotypes)], F), collapse=" , ")))],x)) #add phenotypes as needed
        # The dQuote adds double quotes, the paste0 inserts the commas and cat shows the result without escaping special characters.
        #list2env(traits2, envir = .GlobalEnv)## Assign them to the global environment
        names(traits2)[dim(traits2)[2]]<-"GFS"
        
        # formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        formilin<<-as.formula(paste("GFS", paste(phenotypes[c(3,2)], collapse = " + "), sep = " ~ "))
        # formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        # formilin2<<-as.formula(paste("~",colnames(traits.cl)[spp.col.num], sep = ""))
        formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
        
        lmModel <- try(lm(formilin, method = "ML", data = traits2)) #method = "ML",
        
        if (inherits(lmModel, "try-error")) 
          return(c(NA,NA))
        else
          return(c(summary(lmModel)$adj.r.squared, summary(lmModel)$r.squared, summary(lmModel)$coefficients[2,4],summary(lmModel)$coefficients[2,"t value"]))
      })
      outperm1<-as.data.frame(t(outperm1))
      colnames(outperm1)<-c(paste("adj_r_squared", phenotypes[2], sep=" "),
                            paste("r_squared", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "), paste("adj_r_squared", phenotypes[3], sep=" "),
                            paste("r_squared", phenotypes[3], sep=" "), paste("p_value", phenotypes[3], sep=" "),paste("Coef tval", phenotypes[2]),paste("Coef tval", phenotypes[3]))
    
      lmList<- list(outperm1)
      return(outperm1)
      
    })
    
  } else if (length(pheno.col.nums) == 1) {
    perms<<- parLapply(cl, 1:No.perms, function(i,...){
      traits.cl<-traits.cl
      y<-sample(c(1:length(gene.numbers.cl)))
      gene_numbers2 <- gene.numbers.cl[,y]
      outperm1<-apply(gene_numbers2, MARGIN = 1,function(x){
        traits2<<-cbind(traits.cl[,c(phenotypes)],as.numeric(x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c("Spp.Name" , "SSD" , "Av.Body.mass")],x)) #add phenotypes as needed
        # traits2<<-as.data.frame(cbind(traits.cl[,c(cat(paste0(dQuote(phenotypes[1:length(phenotypes)], F), collapse=" , ")))],x)) #add phenotypes as needed
        # The dQuote adds double quotes, the paste0 inserts the commas and cat shows the result without escaping special characters.
        #list2env(traits2, envir = .GlobalEnv)## Assign them to the global environment
        names(traits2)[dim(traits2)[2]]<-"GFS"
        
        formilin<<-as.formula(paste("GFS", paste(phenotypes[2:length(phenotypes)], collapse = " + "), sep = " ~ "))
        # formilin2<<-as.formula(paste("~",colnames(traits.cl)[spp.col.num], sep = ""))
        formilin2<<-as.formula(paste("~",phenotypes[1], sep = ""))
        
        
        # pglsModel<-pgls(formilin, correlation = corBrownian(1, phy = tree.cl, form = formilin2), method = "ML", data = traits2) #add phenotypes as needed
        # pglsModel<-gls(formilin, correlation = corBrownian(1, phy = tree.cl, form = formilin2), method = "ML", data = traits2, ) #add phenotypes as needed
        
        #add phenotypes as needed
        lmModel <- try(lm(formilin, method = "ML", data = traits2))  # method = "ML",
        
        if (inherits(lmModel, "try-error"))
          return(c(NA,NA))
        else
          return(c(summary(lmModel)$adj.r.squared, summary(lmModel)$r.squared, summary(lmModel)$coefficients[2,4],summary(lmModel)$coefficients[2,"t value"]))
      })
      outperm1<-as.data.frame(t(outperm1))
      colnames(outperm1)<-c(paste("adj_r_squared", phenotypes[2], sep=" "),
                            paste("r_squared", phenotypes[2], sep=" "), paste("p_value", phenotypes[2], sep=" "),paste("Coef tval", phenotypes[2]))
      
      
      lmList<- list(outperm1)
      return(outperm1)
      
    })
  }
  
  log_print("Perms LM done", hide_notes = T)
  
  perms<<- as.data.frame(do.call(rbind, perms), rownames=TRUE)
  log_print(colnames(perms))
  
  for (i in 1:length(pheno.col.nums)) {
    
    log_print("creating variables for stats", hide_notes = T)
    
    #traits2<-as.data.frame(cbind(traits.cl[,c(phenotypes[1:length(phenotypes)])],x)) #add phenotypes as needed
    
    # perms<- bind_rows(perms,.id = NULL)
    sps<-nrow(traits.cl)
    DeFr<- nrow(traits.cl) - (1 + (length(phenotypes)-1)) ## the -1 removes the col.name Spp.Name that its not a phenotype.
    
    head(traits.cl)
    #### the 62 are the numbers of spp minus 2 = degrees of freedom
    log_print("tval and R values", hide_notes = T)
    
    Coef.tval1<<-toString(paste("Coef tval", phenotypes[i+1], sep=" "))
    
    head(perms)
    colnames(perms[[1]])
    rownames(perms)
    R_t_value <<- ((perms[,Coef.tval1]))/sqrt((((perms[,Coef.tval1])^2) + DeFr )) # square root of (tvalue^2 / (tvalue^2 + DF))
    log_print(head(R_t_value), hide_notes = T)
    
    R1 <- paste("R t value", phenotypes[i+1], sep=" ")
    perms[R1]<<- R_t_value
    log_print(" benjamini correction", hide_notes = T)
    
    ## benjamini correction
    perms[paste("p.adjusted Benjamini GFS vs ", phenotypes[i+1], sep = "")] <<- p.adjust(perms[,paste("p_value", phenotypes[i+1], sep=" ")], method = "fdr")
    #permsRs<<-permsRs
  }
  log_print("Creating output file", hide_notes = T)
  nameout<- paste("LM_perms", paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits.cl)[1], "Spp", dim(perms)[1], "Perm_GeneFams", "Rs_benjamini", Sys.Date(), "csv", sep = ".")
  log_print(paste("Location and name of your file:", paste(where.Save.it, nameout, sep = "/")), hide_notes = T)
  write.csv(perms, paste(where.Save.it, nameout, sep = "/"))
  perms<<- perms
  
  closeAllConnections()
}

# works with 2 phenos only at the moment
#*************   Change parameters!!!!!!!   ********###
LM.GF.size.perms(traits.cl = traits.filtered, cores = 20, pheno.col.nums = 3, spp.col.num = 4, 
                   gene.numbers.cl = gene.numbers.filtered, where.Save.it = dirSave, No.perms = 1000)
#*************   ************************** ********###

colnames(perms)

head(perms,3)

#########################
##### plot plot plot ####
#########################
DensityPlots<- function(out, outcol ,perms , permscol, traits.cl, spp.col.num, pheno.col.nums, No.perms, where.Save.it, smoothing, Opacity, xaxisMin, XaxisMax, Phenos){
  
  log_print("#########################", hide_notes = T)
  log_print("##### plot plot plot ####", hide_notes = T)
  log_print("#########################", hide_notes = T)
  
  log_print("It is important to set genefamilies as rownames for out and perms objects", hide_notes = T)  
  
  phenotypes<<- colnames(x = traits.cl)[c(spp.col.num, pheno.col.nums)]
  log_print("Phenos in...", hide_notes = T)
  
  # R1<-paste("R t value", phenotypes[1+1], sep=" ")
  #R1<-paste("R.t.value", phenotypes[1+1], sep=".")
  R1<-paste("R_t_value", phenotypes[1+1], sep=".") # kililis
  # P1<-paste("R t value", phenotypes[1+1], sep=" ")
  P1<-paste("R.t.value", phenotypes[1+1], sep=".")
  
  log_print("Change R1 if you there is an issue here", hide_notes = T)
  log_print(R1, hide_notes = T)
  log_print(P1, hide_notes = T)
  
  # C1<-data.frame(rownames(out),out[R1])
  # C2<-data.frame(rownames(perms),perms[,P1])
  
  C1<-data.frame(rownames(out),out[,outcol])
  C2<-data.frame(rownames(perms),perms[,permscol])
  
  dat<-merge(x = C1, y = C2, by.x = "rownames.out.", by.y = "rownames.perms.", all.y = T)
  dat$rownames.out.<-NULL
  
  #dat<-data.frame(perms[P1])
  #dat<-data.frame(out[R1])
  
  log_print("Setting columns [][]", hide_notes = T)
  #for two variables pĺot
  colnames(dat) <- c(paste("Real R²", phenotypes[1+1], sep=" "),paste("Perms R²", phenotypes[1+1], sep=" "))
  #colnames(dat) <- c(paste("Real R²", phenotypes[1+1], sep=" "))
  #colnames(dat) <- c(paste("Real R²", phenotypes[1+1], sep=" "),paste("Perms R²", phenotypes[1+1], sep=" "))
  #colnames(dat) <- c(paste("Perms R²", phenotypes[1+1], sep=" "))
  dat2<-melt(dat)
  
  # dat<-cbind(out[[1]][,6],reales[,6])
  # colnames(dat) <- c("perm", "real")
  # dat2<-melt(analisis$R_F_value.L_Ei)
  log_print("Name of the plot", hide_notes = T)
  #nameout.plot<- paste("PGLS_perms", "1000",paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits.filtered)[1], "Spp", dim(perms)[1], "Perm_GeneFams", "Rs_benjamini", Sys.Date(), "pdf", sep = ".")
  nameout.plot<- paste("LM_perms", Phenos,"Phenotype",No.perms,"perms",paste(phenotypes[2:length(phenotypes)], collapse  = "_"), dim(traits.cl)[1], "Spp", dim(perms)[1], "Perm_GeneFams", "Rs_benjamini", Sys.Date(), "pdf", sep = ".")
  log_print(nameout.plot, hide_notes = T)
  pdf(paste(where.Save.it, nameout.plot, sep = "/")) ### not yet possible to save it from
  #pdf(paste(dirSave, nameout.plot, sep = "/")) ### not yet possible to save it from
  
  #xaxisMin<-(max(out[R1])+max(perms[P1])/2)+((max(out[R1])+max(perms[P1])/2)*0.10)
  #XaxisMax<-(max(out[R1])+max(perms[P1])/2)+((max(out[R1])+max(perms[P1])/2)*0.10)
  
  print("Building plot")
  plot1<-ggplot(dat2, aes(x = value, fill = variable)) + geom_density(alpha= Opacity, na.rm = T, adjust = smoothing, ) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) + xlim(xaxisMin, XaxisMax)
  
  # plot1<-ggplot(dat2, aes(x=value, fill = variable)) + stat_function(fun = dnorm, geom = "density", xlim = c(-5, 5), fill = rgb(0, 0, 1, 0.1)) + 
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
  #         axis.line = element_line(colour = "black"))
  
  
  
  print("Saving plot")
  print(plot1)
  dev.off()
}

#out<-read.csv("PGLS.cEi.46.Spp.4121.GeneFams.Rs_benjamini.2023-02-10.csv", row.names = 1)

#*************   Change parameters!!!!!!!   ********###
DensityPlots(out = out, outcol = "R t value Estimated_longevity_days_to_normal", perms = perms, permscol= "R t value Estimated_longevity_days_to_normal", 
             traits.cl = traits.filtered, spp.col.num = 4, pheno.col.nums = 3, 
             No.perms = 1000, Opacity = 0.1 , smoothing = 1.5, xaxisMin= -0.8, XaxisMax = 0.8, 
             Phenos = "Estimated_longevity_days_to_normal", where.Save.it = dirSave)
#*************   ************************** ********###

closeAllConnections()

colnames(out)

log_close()
writeLines(readLines(lf))