########################################################
# Sceptical prior construction
########################################################
#
# started 16/03/2018
#
# Normal approximation to the elicited expert opinions
# 
#
########################################################
sceptical <- function(expertnb, combined=F, nSubg=3,marginal=F, CovSens=F, MeanSens=F,
                      Arm1trt = Arm1trt, Arm1base = Arm1base, Arm2trt = Arm2trt, Arm2base = Arm2base, Arm3trt = Arm3trt, Arm3base = Arm3base) {
  source('C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Var-Cov matrix.R')
  if (combined==T){ a <- VARCOV(combined=T, marginal=marginal, CovSens=CovSens, MeanSens=MeanSens)
  } else {a <- VARCOV(expertnb = expertnb, marginal=marginal, CovSens=CovSens, MeanSens=MeanSens)}
  
  sceptical <- paste("list(mutheta=c(",paste(c(rep("0,", nSubg-1),0), collapse =  ""), "), ",
                        "diagv=c(" ,paste(a$V11,a$V22,a$V33,sep = ", "), "), ",
                        "v12=",a$V12, ", v13=",a$V13,", v23=",a$V23 , ",
                        nArms = ",2*nSubg,", ntrt=",4,", nSubg=", nSubg,
                        ",
                 Subg=c(", paste(c(rep("1,", nSubg-1)), collapse =  "") , paste(c(rep("2,", nSubg-1)), collapse =  ""),paste(c(rep("3,", nSubg-2),3), collapse =  ""),")",
                        ",
                 trt=c(2,1,3,1,4,1),  # trt[1] Placebo trt[2] ADDED trt[3] ALT trt[4] PRESERVED",
                        "
                 r=c(", paste(Arm1trt[1],Arm1base[1],   Arm2trt[1],Arm2base[1], Arm3trt[1],Arm3base[1],  sep = ","), "), # number of deaths or first hosp admissions",
                        "
                 n=c(" ,paste(Arm1trt[2],Arm1base[2],   Arm2trt[2],Arm2base[2], Arm3trt[2],Arm3base[2], sep = ","), ") # number of person at risk","
  )",
                        sep = ""    
  )
  write.table(sceptical,quote = F,row.names  = F,col.names  = F,
              file="C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/model 5 - CHARM with priors/sceptical/data.txt")
  
  # return(cat(sceptical, is.positive.definite(matrix(c(a$V11,a$V12,a$V13,
  #                                                     a$V12,a$V22,a$V23,
  #                                                     a$V13,a$V23,a$V33), byrow = T, nrow = 3))))
  
}

