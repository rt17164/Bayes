posterior_plot <-  function(expertnb=expertnb, combined=F) { 
  if (combined==T){
    expert_mode1 <<- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG1","/expertcombUnconditionalG1", ".csv",sep = ""), sep = ",",header = T)
    expert_mode2 <<- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG2","/expertcombUnconditionalG2", ".csv",sep = ""), sep = ",",header = T)
    expert_mode3 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G2","/expertcombConditional1G2", ".csv",sep = ""), sep = ",",header = T)
    expert_mode4 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G2","/expertcombConditional2G2", ".csv",sep = ""), sep = ",",header = T)
    expert_mode5 <<- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG3","/expertcombUnconditionalG3", ".csv",sep = ""), sep = ",",header = T)
    expert_mode6 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G3","/expertcombConditional1G3", ".csv",sep = ""), sep = ",",header = T)
    expert_mode7 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G3","/expertcombConditional2G3", ".csv",sep = ""), sep = ",",header = T)
  } else {
    expert_mode1 <<- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG1","/expert", expertnb, "UnconditionalG1", ".csv",sep = ""), sep = ",",header = T)
    expert_mode2 <<- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG2","/expert", expertnb, "UnconditionalG2", ".csv",sep = ""), sep = ",",header = T)
    expert_mode3 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G2","/expert", expertnb, "Conditional1G2", ".csv",sep = ""), sep = ",",header = T)
    expert_mode4 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G2","/expert", expertnb, "Conditional2G2", ".csv",sep = ""), sep = ",",header = T)
    expert_mode5 <<- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/UnconditionalG3","/expert", expertnb, "UnconditionalG3", ".csv",sep = ""), sep = ",",header = T)
    expert_mode6 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional1G3","/expert", expertnb, "Conditional1G3", ".csv",sep = ""), sep = ",",header = T)
    expert_mode7 <- read.csv(file=paste(computerwd,"Bristol/NIHR RMF/GAP Bayesian/GAP Subgroup/BUGS input/Elicitations/By mode/Conditional2G3","/expert", expertnb, "Conditional2G3", ".csv",sep = ""), sep = ",",header = T)
  }
  x <- c(1,2,3)
  y <- c(seq(1,6,1))
  y1 <- vague <- c(vag$summary[4],vag$summary[5],vag$summary[6])

  if (combined){
  y2 <- enthusiastic <- c(entComb$summary[4],entComb$summary[5],entComb$summary[6])
  y3<- sceptical <- c(scepComb$summary[4],scepComb$summary[5],scepComb$summary[6])
  y4 <- interaction <- c(intComb$summary[4],intComb$summary[5],intComb$summary[6])
  y5 <- interaction2 <- c(int2Comb$summary[4],int2Comb$summary[5],int2Comb$summary[6])
  } else {   y2 <- enthusiastic <- c(ent$summary[4],ent$summary[5],ent$summary[6])
  y3<- sceptical <- c(scep$summary[4],scep$summary[5],scep$summary[6])
  y4 <- interaction <- c(int$summary[4],int$summary[5],int$summary[6])
  y5 <- interaction2 <- c(int2$summary[4],int2$summary[5],int2$summary[6])}
  
  x0 <- c("ADDED", "ALT", "PRES")
  
  priormean <- as.numeric(format(round(exp(c(expert_mode1[1,c("LHR.mean")],
                                             expert_mode2[1,c("LHR.mean")],
                                             expert_mode5[1,c("LHR.mean")])),3),nsmall = 2))
  CI_prior_L <- as.numeric(format(round(exp(c(expert_mode1[1,c("LHR.mean")]-1.96*expert_mode1[1,c("LHR.sd")],
                            expert_mode2[1,c("LHR.mean")]-1.96*expert_mode2[1,c("LHR.sd")],
                            expert_mode5[1,c("LHR.mean")]-1.96*expert_mode5[1,c("LHR.sd")])),3),nsmall = 2))
  CI_prior_U <- as.numeric(format(round(exp(c(expert_mode1[1,c("LHR.mean")]+1.96*expert_mode1[1,c("LHR.sd")],
                            expert_mode2[1,c("LHR.mean")]+1.96*expert_mode2[1,c("LHR.sd")],
                            expert_mode5[1,c("LHR.mean")]+1.96*expert_mode5[1,c("LHR.sd")])),3),nsmall = 2))
  Pbenefit <- as.numeric(format(round(1-pnorm(c(expert_mode1[1,c("LHR.mean")]/expert_mode1[1,c("LHR.sd")],
                expert_mode2[1,c("LHR.mean")]/expert_mode2[1,c("LHR.sd")],
                expert_mode5[1,c("LHR.mean")]/expert_mode5[1,c("LHR.sd")])),3),nsmall = 3))

  Prior <- c(rep(1:6, each=3))
  Subgroup <- c(rep(1:3, times=6))
  HR <- c(priormean, vague, enthusiastic,sceptical, interaction,interaction2)
  
  df <- data.frame(Prior, Subgroup, HR)
  
  nPriors <- max(df$Prior)
  
  xrange <- range(df$Subgroup)
  if (combined){ yrange <- c(0.75,0.90) } else{ yrange <- c(0.70,0.90) } 
  
  # set up the plot 
  plot(xrange, yrange , type="n", xlab="",
       ylab="Hazard ratio", axes = F ) 
  
  colors <- c("blue3", "brown4",  "palegreen4", "goldenrod2", "lightcyan4","indianred3") 
  linetype <- c(3,5,1,1,1,1) 
  plotchar <- c(1,1,1,5,2,7)
  for (i in c(seq(yrange[1],yrange[2],0.05))) {abline(a=i, b=0, col="gainsboro") }
  
  # add lines 
  
  for (i in 1:nPriors) { 
    priortype <- subset(df, Prior==i) 
    lines(priortype$Subgroup, priortype$HR, type="b", lwd=2,
          lty=linetype[i], col=colors[i], pch=plotchar[i] )
  } 
  
  axis(2,at=c(seq(0.7,1.0,0.05)))
  axis(1,at=c(seq(1,3,1)), labels=c("ADDED", "ALTERNATIVE", "PRESERVED"))
  
  # add a title and subtitle 
  if (combined) {title("Posterior median treatment effect on primary outcome in each
CHARM trial using the pooled experts' prior in different ways")
} else {title(paste("Posterior median treatment effect on primary outcome in each
CHARM trial using the prior from expert ", expertnb, " in different ways", sep=""))}

  
  # add a legend 
  legend("bottomleft", 1:nPriors, cex=0.75, col=colors,
         pch=plotchar, lty=linetype, title="Prior", 
         legend = c("prior", "vague", "enthusiastic","sceptical", "interaction","interaction-variance"),
         lwd=1.5, inset=c(0.01,0.01)
  )
  
  prior=matrix(c(CI_prior_L,priormean,CI_prior_U,Pbenefit),nrow = 3, byrow=F) ; 
  colnames(prior) <- c("2.5%", "50%","97.5%", "P(benefit)") ; 
  rownames(prior) <- c("hr.ADDED", "hr.ALT", "hr.PRES")
  
  
  if (combined){ return(list(prior=prior, vague=vag$summary, enthusiastic=entComb$summary, sceptical=scepComb$summary, interaction=intComb$summary, interaction2=int2Comb$summary,
                            Marginal=intCombMarg$summary, ConditionalM=intCombMean$summary, Cov1=intCombCov1$summary, Cov2=intCombCov2$summary))
  } else { return(list(prior=prior, vague=vag$summary, enthusiastic=ent$summary, sceptical=scep$summary, interaction=int$summary, interaction2=int2$summary ))}
  
}
