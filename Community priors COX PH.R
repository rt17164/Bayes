########################################################
# Community of priors (Cox PH model)
########################################################
# 
# started 30/05/2018
# 
# 
#
########################################################

# packages
library(foreign)
library(coda)
library(crayon)
library(matrixcalc)
library(survival)

############
work = T ; if (work==T) {computerwd <- "C:/Users/rt17164/Google Drive/"} else {computerwd <- "/Users/russt/Google Drive/"}
setwd(paste(computerwd, "Bristol/NIHR RMF/GAP Bayesian/Simulations/BUGS input", sep = ""))
start <- proc.time()

# Trial Data
source(paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Simulations CHARM v0.1.R", sep=""))
"C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Simulations CHARM v0.1.R"
"/Users/russt/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Simulations CHARM v0.1.R"

# find best dataset to mimic CHARM trials
results1 <- results[["n = 2000, baseline = Exp, model = Cox"]]
best <- results1[results1[,"logHRdif"]==min(results1$logHRdif),]$dataset
data1 <- data[["n = 2000, baseline = Exp"]][[best]]

# gen Subgroup variable
ADDED <- ifelse(data1$ADDED==1,1,0)
ALT <- ifelse(data1$ALT==1,2,0)
PRES <- ifelse(data1$PRES==1,3,0)
data1$Subg <- ADDED+ALT+PRES

# gen Subgroup specific treatment effect (4 seperate treatment effects: in placebo group, ADDED, ALT and PRES)
trtbase <- ifelse(data1$trt==0,1,0)
trtADDED <- ifelse(data1$trt==1 &  data1$ADDED==1,2,0)
trtALT <- ifelse(data1$trt==1 &  data1$ALT==1,3,0)
trtPRES <- ifelse(data1$trt==1 &  data1$PRES==1,4,0)
data1$trt1 <- trtbase+ trtADDED + trtALT + trtPRES
prop.table(table(data1$trt1))

## Modeling Cox PH
data1[1:13,]
write.dta(data1, file =paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Simulations/data.dta", sep=""))
fit <- survival::coxph(Surv(t, d) ~  factor(Subg)+factor(trt1) , data = data1)
# HR <- list(HRADDED= as.numeric(exp(svycontrast(fit,c("trt"=1 ))[1])) ,                     # HR TRT Thor VS Pla Thor
#            HRALT= as.numeric(exp(svycontrast(fit,c("trt"=1, "trt:ALT"=1 ))[1])) ,          # HR TRT Abdo VS Pla Abdo
#            HRPRES= as.numeric(exp(svycontrast(fit,c("trt"=1, "trt:PRES"=1 ))[1])))         # HR TRT Card VS Pla Card
# HR


#######################################################################################

# Data
# Store values of model
DataInfo <- list()
DataInfo$N <- nrow(data1)
DataInfo$T <- length(unique(data1$t))-1
DataInfo$eps <- 1E-10
DataInfo$obs.t <- data1$t
DataInfo$d <- data1$d
DataInfo$trt1 <- data1$trt1
DataInfo$t <- sort(unique(data1$t))
DataInfo$ntrt <- 4
DataInfo$nSubg <- 3
DataInfo$nArms <- 2 *  DataInfo$nSubg
DataInfo$Subg <- data1$Subg


# Analysis script
source(paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Parallelise OpenBUGS v0.3.R", sep=""))
"C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Parallelise OpenBUGS v0.3.R"
post <- list()
 for (expertnb in 1:1){
   # summvague <- summ
   vag <- RunParallelBUGS(prior = "vague",expertnb = expertnb,  DataInfo= DataInfo, work=work)
   ent <- RunParallelBUGS(prior = "enthusiastic",expertnb = expertnb,  DataInfo= DataInfo, work=work)
   scep <- RunParallelBUGS(prior = "sceptical",expertnb = expertnb,  DataInfo= DataInfo, work=work)
   int <- RunParallelBUGS(prior = "interaction",expertnb = expertnb,  DataInfo= DataInfo, work=work)
   int2 <- RunParallelBUGS(prior = "interaction2",expertnb = expertnb,  DataInfo= DataInfo, work=work)

   # Posterior plot
   source(paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Figure 3 & 4.R", sep=""))
   "C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Figure 3 & 4.R"
   post[[expertnb]] <-  posterior_plot(expertnb = expertnb)
}
 
# 
# # Combined
# vag <- RunParallelBUGS(prior = "vague",combined=T,  DataInfo= DataInfo)
# entComb <- RunParallelBUGS(prior = "enthusiastic",combined=T,  DataInfo= DataInfo)
# scepComb <- RunParallelBUGS(prior = "sceptical",combined=T,  DataInfo= DataInfo)
# intComb <- RunParallelBUGS(prior = "interaction",combined=T,  DataInfo= DataInfo)
# int2Comb <- RunParallelBUGS(prior = "interaction2",combined=T,  DataInfo= DataInfo)
# 
# # Alternative Prior constructions
# intCombMarg <- RunParallelBUGS(prior = "interaction",combined=T,  DataInfo= DataInfo, marginal = T)
# intCombMean <- RunParallelBUGS(prior = "interaction",combined=T,  DataInfo= DataInfo, MeanSens = T)
# intCombCov1 <- RunParallelBUGS(prior = "interaction",combined=T,  DataInfo= DataInfo, CovSens = 1)
# intCombCov2 <- RunParallelBUGS(prior = "interaction",combined=T,  DataInfo= DataInfo, CovSens = 2)
# 
# 
# 
# # Posterior plot 
# source(paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Figure 3 & 4.R", sep=""))
# "C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Figure 3 & 4.R"
# 
# postComb <- posterior_plot(combined = T)
# write.table(postComb, file = paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS output/entExpertComb.txt", sep = ""))

elapsedT <- (proc.time()-start)
cat("\n","The analysis took ",
    floor(elapsedT[3]/3600), "hours," ,
    floor((elapsedT[3]/3600-floor(elapsedT[3]/3600))*60), "minutes and",
    round(elapsedT[3] - floor(elapsedT[3]/60)*60), "seconds")

#######################################################################################