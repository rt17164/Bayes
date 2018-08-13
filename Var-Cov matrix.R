########################################################
# derive Var-Cov matrix elements from normal approximation
########################################################
# 
# started 14/03/2018
# 
# 
#
########################################################
VARCOV <- function(expertnb, marginal=F, CovSens=F, MeanSens=F, work=T, combined=F ){
  
# difference in elicitation
d=log(log(1-0.12)/log(1-0.18))
  
if (work==T) {computerwd <- "C:/Users/rt17164/Google Drive/"} else {computerwd <- "/Users/russt/Google Drive/"}
if (combined==T){
expert_mode1 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG1","/expertcombUnconditionalG1", ".csv",sep = ""), sep = ",",header = T)
expert_mode2 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG2","/expertcombUnconditionalG2", ".csv",sep = ""), sep = ",",header = T)
expert_mode3 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G2","/expertcombConditional1G2", ".csv",sep = ""), sep = ",",header = T)
expert_mode4 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G2","/expertcombConditional2G2", ".csv",sep = ""), sep = ",",header = T)
expert_mode5 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG3","/expertcombUnconditionalG3", ".csv",sep = ""), sep = ",",header = T)
expert_mode6 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G3","/expertcombConditional1G3", ".csv",sep = ""), sep = ",",header = T)
expert_mode7 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G3","/expertcombConditional2G3", ".csv",sep = ""), sep = ",",header = T)
} else {
expert_mode1 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG1","/expert", expertnb, "UnconditionalG1", ".csv",sep = ""), sep = ",",header = T)
expert_mode2 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG2","/expert", expertnb, "UnconditionalG2", ".csv",sep = ""), sep = ",",header = T)
expert_mode3 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G2","/expert", expertnb, "Conditional1G2", ".csv",sep = ""), sep = ",",header = T)
expert_mode4 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G2","/expert", expertnb, "Conditional2G2", ".csv",sep = ""), sep = ",",header = T)
expert_mode5 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG3","/expert", expertnb, "UnconditionalG3", ".csv",sep = ""), sep = ",",header = T)
expert_mode6 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G3","/expert", expertnb, "Conditional1G3", ".csv",sep = ""), sep = ",",header = T)
expert_mode7 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G3","/expert", expertnb, "Conditional2G3", ".csv",sep = ""), sep = ",",header = T)
}


# logHR mean
p1 <- expert_mode1[1,"LHR.mean"]
p2 <- expert_mode2[1,"LHR.mean"]
p3 <- expert_mode5[1,"LHR.mean"]
p20 <- expert_mode3[1,"LHR.mean"]
p2d <- expert_mode4[1,"LHR.mean"]
p30 <- expert_mode6[1,"LHR.mean"]
p3d  <- expert_mode7[1,"LHR.mean"]

# logHR var
s1sq <- expert_mode1[1,"LHR.sd"]^2
s2sq <- expert_mode2[1,"LHR.sd"]^2
s3sq <- expert_mode5[1,"LHR.sd"]^2
s20sq <- expert_mode3[1,"LHR.sd"]^2
s2dsq <- expert_mode4[1,"LHR.sd"]^2
s30sq <- expert_mode6[1,"LHR.sd"]^2
s3dsq <- expert_mode7[1,"LHR.sd"]^2


# step 1
V11 <-  s1sq

# step 2

b12 <- (p2d-p20)/d
b21 <- b12
V12 <- V21 <- b21*V11

b13 <- (p3d-p30)/d
b31 <- b13
V13 <- V31 <- b31*V11

### Sensitivity ###
if (MeanSens==T){
  # prespecified conditional mean method
  p2 <- p20+b21*p1
  p2 <- p2d+b21*p1}

# step 3
if (marginal==T) {
# (a) marginal method
  V22 <- s2sq
  V33 <- s3sq
} else {
# (b) conditional method
  varCond1G2 <- (s20sq+s2dsq)/2  # average of conditional variances in group 2 given group 1
  varCond1G3 <- (s30sq+s3dsq)/2  # average of conditional variances in group 3 given group 1
  V22 <-  varCond1G2 + b12^2*V11
  V33 <-  varCond1G3 + b13^2*V11
}

# step 4
V23 <- V32 <- V12*V13/V11

if(CovSens==1) {
  # covariance between theta2 and theta3 <- covariance between theta1 and theta2
  V23 <- V12
}
if(CovSens==2) {
  # covariance between theta2 and theta3 <- covariance between theta1 and theta3
  V23 <- V13
}



vect <- c(p1,p2,p3)
mat <- matrix(c(V11,V12,V13,V12, V22, V23,V13, V23,V33), nrow = 3, byrow = F)

VARCOV <- list(p1=p1, p2=p2, p3=p3, V11=V11, V22=V22, V33=V33,V12=V12,V13=V13,V23=V23)

return(VARCOV)
}


