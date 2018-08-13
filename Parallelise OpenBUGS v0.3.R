########################################################
# Parallelise OpenBUGS using Cox PH model
########################################################
# 
# started 11/06/2018
# 
# 2. with alpha which is the subgroup specific effect on survival (Bayes Cox PH model)
# 3. different priors (ent, scep, int, int2)
# 
########################################################
# prior = "vague";  DataInfo= DataInfo; combined=F; marginal=F; CovSens=F; MeanSens=F; expertnb=1
RunParallelBUGS <- function(prior, expertnb, combined=F, marginal=F, CovSens=F, MeanSens=F, DataInfo, work){
  ptm <- proc.time()
  if (work==F){setwd("~/")}
  # setwd(paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Simulations/BUGS input", sep=""))
  
  # 1. loading the libraries for parallel processing
  library(snowfall)
  
  # 2. setting the number of CPUs to be 2 / 4 / 6
  sfInit(parallel=TRUE, cpus=6)
  
  # 3. and assigning the R2OpenBUGS library to each CPU
  sfLibrary(R2OpenBUGS)
  
  # 4. generation DataInfo file
  if (prior=="vague"){DataInfo$mutheta <- c(rep(0,DataInfo$nSubg)) ; 
  DataInfo$diagv <- c(1000,1000,1000); DataInfo$V12=0; DataInfo$V13=0; DataInfo$V23=0}
  
  source(paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Var-Cov matrix.R", sep=""))
  "C:/Users/rt17164/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Var-Cov matrix.R"
  "/Users/russt/Google Drive/Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/Var-Cov matrix.R"
  if (combined==T){ a <- VARCOV(combined=T,marginal=marginal, CovSens=CovSens, MeanSens=MeanSens, work=work)} else {
    a <- VARCOV(expertnb = expertnb,marginal=marginal, CovSens=CovSens, MeanSens=MeanSens,work=work)}
  
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
  
  
  if (prior=="enthusiastic"){ DataInfo$mutheta <- c(a$p1,a$p2,a$p3) ; 
  DataInfo$diagv <- c(a$V11,a$V22,a$V33); DataInfo$V12=a$V12; DataInfo$V13=a$V13; DataInfo$V23=a$V23}
  if (prior=="sceptical"){DataInfo$mutheta <- c(rep(0,DataInfo$nSubg)) ; 
  DataInfo$diagv <- c(a$V11,a$V22,a$V33); DataInfo$V12=a$V12; DataInfo$V13=a$V13; DataInfo$V23=a$V23 }
  if (prior=="interaction"){
    DataInfo$mutheta=c(a$p1,a$p2,a$p3); 
    DataInfo$diagv=c(V_star[1,1],V_star[2,2],V_star[3,3]);
    DataInfo$V12=V_star[1,2];
    DataInfo$V13=V_star[1,3];
    DataInfo$V23=V_star[2,3] }
  if (prior=="interaction2"){  
    DataInfo$mutheta=c(rep(0,DataInfo$nSubg)); 
    DataInfo$diagv=c(V_star[1,1],V_star[2,2],V_star[3,3]);
    DataInfo$V12=V_star[1,2];
    DataInfo$V13=V_star[1,3];
    DataInfo$V23=V_star[2,3] }
  
  BUGS.data  <- list(N=DataInfo$N, T=DataInfo$T, eps=DataInfo$eps, ntrt=DataInfo$ntrt, 
                     nSubg=DataInfo$nSubg, mutheta=DataInfo$mutheta, diagv=DataInfo$diagv , 
                     V12=DataInfo$V12, V13=DataInfo$V13, V23=DataInfo$V23, t=DataInfo$t,
                     obs.t=DataInfo$obs.t, d=DataInfo$d, trt1=DataInfo$trt1, Subg=DataInfo$Subg)
  
  bugs.data(BUGS.data, digits = 3,dir =paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/Simulations/BUGS input", sep = ""), data.file = "data.txt")
  
  # 5. creating separate directory for each CPU process
  folder1 <- paste(getwd(), "/chain1", sep="")
  folder2 <- paste(getwd(), "/chain2", sep="")
  dir.create(folder1); dir.create(folder2);
  
  # 6. sinking the model into a file in each directory
  for (folder in c(folder1, folder2)){
    sink(paste(folder, "/model.txt", sep=""))
    cat("
        model
        {
        # Set up data
        for(i in 1:N) {
        for(j in 1:T) {
        # risk set = 1 if obs.t >= t
        Y[i,j] <- step(obs.t[i] - t[j] + eps)
        # counting process jump = 1 if obs.t in [ t[j], t[j+1] )
        #                      i.e. if t[j] <= obs.t < t[j+1]
        dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * d[i]
        }
        }
        
        # Priors
        alpha[1] <- 0
        theta[1] <- 0.00000E+00
        theta[2:ntrt] ~ dmnorm(mutheta[1:nSubg], prectheta[1:nSubg, 1:nSubg])
        for (i in 1:nSubg) {
        V[i, i] <- diagv[i]
        }
        V[1, 2] <- V12
        V[1, 3] <- V13
        V[2, 3] <- V23
        for (i in 2:nSubg) {
        for (j in 1:(i - 1)) {
        V[i, j] <- V[j, i]
        }
        }
        prectheta[1:nSubg, 1:nSubg] <- inverse(V[1:nSubg, 1:nSubg])
        
        for (i in 2:nSubg){alpha[i]~dnorm(0,0.0001)}
        
        # Likelihood 
        for(j in 1:T) {
        for(i in 1:N) {
        dN[i, j]   ~ dpois(Idt[i, j])              # likelihood
        Idt[i, j] <- Y[i, j] * exp(alpha[Subg[i]] + theta[trt1[i]]) * dL0[j] 	# Intensity 
        }     
        dL0[j] ~ dgamma(mu[j], c)
        mu[j] <- dL0.star[j] * c    # prior mean hazard
        
        # Survivor function = exp(-Integral{l0(u)du})^exp(theta*z)    
        # S.treat[j] <- pow(exp(-sum(dL0[1 : j])), exp(theta * -0.5));
        # S.placebo[j] <- pow(exp(-sum(dL0[1 : j])), exp(theta * 0.5));	
        }
        c <- 0.001
        r <- 0.1
        for (j in 1:T) {
        dL0.star[j] <- r * (t[j + 1] - t[j])
        }
        
        # Derivation
        for (iv in 1:ntrt) {
        # Calculate HRs
        hr[iv] <- exp(theta[iv])
        # Calculate Bayesian p-values - Probability of benefit p(theta)<0
        p0[iv] <- step(0-theta[iv])
        }
        
        # Dummy variables
        hr.ADDED <- hr[2]
        hr.ALT <- hr[3]
        hr.PRES <- hr[4]
        p.ADDED <- p0[2]
        p.ALT <- p0[3]
        p.PRES <- p0[4]
        
        # Ranking plot
        for (v in 1:ntrt) {
        for (vi in 1:ntrt) {
        rk[v,vi] <- equals(ranked(theta[],vi),theta[v])
        }
        }

        } # model end")
sink()
}
  
  # 7. defining the function that will run MCMC on each CPU
  parallel.bugs <- function(chain, BUGS.data, pars, DataInfo, work){
    # 7a. defining directory for each CPU
    sub.folder <- paste(getwd(),"/chain", chain, sep="")
    
    # 7b. specifying the initial MCMC values
    inits <- function()
    {
      list(theta = c(NA,rep(rbinom(size=1,prob = 0.5,n=1)*0.5 ,3 )),
           alpha = c(NA,rep(runif(1,-0.5,0.5),DataInfo$nSubg-1)),
           dL0 = c(rep(round(runif(1,2,n = 1)),DataInfo$T)))
    }
    
    # 7c. calling OpenBugs
    # (you may need to change the OpenBUGS.pgm directory)
    if (work){bugs(data=BUGS.data, inits=inits,
                   parameters.to.save=pars, model.file="model.txt",
                   n.chains=1, n.iter=6000, n.burnin=2500, n.thin=1,
                   debug=F, DIC=F, codaPkg=T,
                   working.directory = sub.folder)} else {
              setwd("/Users/russt")
              bugs(data=BUGS.data, inits=inits,
                   parameters.to.save=pars, model.file="model.txt",
                   n.chains=1, n.iter=6000, n.burnin=2500, n.thin=1,
                   debug=F, DIC=F, codaPkg=T,
                   working.directory = sub.folder,
                   WINE="/Users/russt/Downloads/usr/bin/wine",
                   WINEPATH ="/Users/russt/Downloads/usr/bin/winepath",
                   OpenBUGS.pgm="/Users/russt/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe",
                   useWINE = T)
              }
  }
  
  # 8. setting the parameters to be monitored
  pars<- c("hr.ADDED", "hr.ALT", "hr.PRES", "rk", "p.ADDED", "p.ALT", "p.PRES")
  
  
  # 9. calling the sfLapply function that will run
  # parallel.bugs on each of the 3 CPUs
  sfLapply(1:2, fun=parallel.bugs, BUGS.data=BUGS.data, pars=pars, DataInfo=DataInfo, work=work)
  sfStop()
  # 10. locating position of each CODA chain
  chain1 <- paste(folder1, "/CODAchain1.txt", sep="")
  chain2 <- paste(folder2, "/CODAchain1.txt", sep="")
  
  # 11. and, finally, getting the results
  line.coda <- read.bugs(c(chain1, chain2))
  line.coda.hr <- mcmc.list(line.coda[[1]][,c("hr.ADDED", "hr.ALT", "hr.PRES")],
                            line.coda[[2]][,c("hr.ADDED", "hr.ALT", "hr.PRES")])
  summ <- summary(line.coda, digits=3)
  stats <- as.data.frame(summ$statistics)
  quantiles <- as.data.frame(summ$quantiles)
  hr.sum <- matrix(as.numeric(c(format(round(summ$quantiles[c("hr.ADDED", "hr.ALT", "hr.PRES"),c("2.5%", "50%" , "97.5%" )],3),nsmall=2),
                                format(round(summ$statistics[c("p.ADDED", "p.ALT", "p.PRES"),c("Mean" )],3),nsmall = 3))),nrow = 3,byrow = F)
  
  rownames(hr.sum) <- c("hr.ADDED", "hr.ALT", "hr.PRES" )
  colnames(hr.sum) <- c("2.5%", "50%" , "97.5%", "P(benefit)" )
  hr.sum
  
  # reduce line.coda to hr
  head(line.coda[[1]][,grepl("rk", dimnames(line.coda[[1]])) ] )
  
  HPDinterval <- HPDinterval(mcmc.list(line.coda[[1]][,c("hr.ADDED", "hr.ALT", "hr.PRES")],
                                       line.coda[[2]][,c("hr.ADDED", "hr.ALT", "hr.PRES")])) # results of both chains
  plot <- plot(line.coda.hr)
  gelman.diag(line.coda.hr )
  gelman.plot(line.coda.hr)
              
              
  # # ranking plot
  test <- matrix(stats[grepl("rk", rownames(stats)),"Mean"],nrow = 4, byrow = F)
  colnames(test) <- c("Placebo","CHARM Added","CHARM Alternative", "CHARM Preserved")
  rownames(test) <- c(paste("Ranking 1", "Most effective", sep = "\n"), "Ranking2", "Ranking3", paste("Ranking 4", "Least effective", sep = "\n"))
  source(paste(computerwd, "Bristol/NIHR RMF/GAP Bayesian/Histograms.Shiny/ggcorrplot_adapt.R" , sep = ""))
  rkplot <- ggcorrplot_adapt(test, lab = T, lab_col = 1)

  elapsedT <- (proc.time()-ptm)
  if (combined==F) {cat("\n","Running BUGS model using parallel using", prior,"prior from expert", expertnb, "took",
                        floor(elapsedT[3]/60), " minutes and ", round(elapsedT[3] - floor(elapsedT[3]/60)*60), "seconds","\n")} else {
                          cat("\n","Running BUGS model using parallel using", prior,"prior from combined experts took", 
                              floor(elapsedT[3]/60), " minutes and ", round(elapsedT[3] - floor(elapsedT[3]/60)*60), "seconds","\n")
                        }
  return(list(summary=hr.sum, BUGSchains=line.coda, RankingPlot=rkplot))
}


