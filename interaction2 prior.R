########################################################
# Interraction 2 prior construction 
# (only using interraction with no main trt effect)
########################################################
# 
# started 20/03/2018
#
########################################################
interaction2<- function(expertnb,combined=F, nSubg=3,marginal=F, CovSens=F, MeanSens=F,
                        Arm1trt = Arm1trt, Arm1base = Arm1base, Arm2trt = Arm2trt, Arm2base = Arm2base, Arm3trt = Arm3trt, Arm3base = Arm3base) {
  source('C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Var-Cov matrix.R')
  if (combined==T){ a <- VARCOV(combined=T,marginal=marginal, CovSens=CovSens, MeanSens=MeanSens)
  } else {a <- VARCOV(expertnb = expertnb,marginal=marginal, CovSens=CovSens, MeanSens=MeanSens)}
  C <- matrix(c( 1/3 , 1/3 , 1/3 ,
                 -1  ,  1  ,  0  ,
                 -1  ,  0  ,  1  ), byrow = T, nrow = 3)
  
  V <- matrix(c(a$V11,a$V12,a$V13,
                a$V12,a$V22,a$V23,
                a$V13,a$V23,a$V33), byrow = T, nrow = 3)
  
  W <- C %*% V %*% t(C)
  
  W_star <-  matrix(c( 1000 ,    0    ,     0    ,
                       0  , W[2,2]  ,  W[2,3]  ,
                       0  , W[3,2]  ,  W[3,3]  ), byrow = T, nrow = 3)
  
  V_star <- solve(C) %*% W_star %*% solve(t(C))
  
  # mean_vector <- c(a$p1, a$p2, a$p3)
  # m4 <- as.list(C %*% mean_vector) ; names(m4) <- c("p1","p2","p3")
  
  interaction2 <- paste("list(mutheta=c(",paste(c(rep("0,", nSubg-1),0), collapse =  ""), "), ",
                        "diagv=c(" ,paste(V_star[1,1],V_star[2,2],V_star[3,3],sep = ", "), "), ",
                        "v12=",V_star[1,2], ", v13=",V_star[1,3],", v23=",V_star[2,3] , ",
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
  write.table(interaction2,quote = F,row.names  = F,col.names  = F,
              file="C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/model 5 - CHARM with priors/interaction2/data.txt")
  # return(cat(interaction2, is.positive.definite(matrix(c(a$V11,a$V12,a$V13,
  #                                                        a$V12,a$V22,a$V23,
  #                                                        a$V13,a$V23,a$V33), byrow = T, nrow = 3))))
}

