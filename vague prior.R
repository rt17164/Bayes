########################################################
# Vague prior construction
########################################################
# 
# started 15/03/2018
#
########################################################
vague <- function(nSubg,Arm1trt = Arm1trt, Arm1base = Arm1base, Arm2trt = Arm2trt, Arm2base = Arm2base, Arm3trt = Arm3trt, Arm3base = Arm3base) {
vague <- paste("list(mutheta=c(",paste(c(rep("0,", nSubg-1),0), collapse =  "") , "), ",
      "diagv=c(" ,paste(c(rep("200,", nSubg-1),200), collapse =  ""), "), ",
      "v12=",0, ", v13=",0,", v23=",0 , ",
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
write.table(vague,quote = F,row.names  = F,col.names  = F,
            file="C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/model 5 - CHARM with priors/vague/data.txt")
# return(cat(vague))

}
  