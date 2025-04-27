tdcm.mean.no.compare <- function(betai, treat,nutnet = F, prefixi,scale.num = FALSE, rand = 1000, alpha = TRUE)
  {
  if (nutnet == F){
  betai <- as.data.frame(betai)
  library(lme4)
  library(car)
  tdc.lmm <- function(betai, treat,scale.num = FALSE, prefixi = NULL) {
    warming <- treat$treat
    warm.lev <- unique(warming)
    out <- list()
    for (j in 1:length(warm.lev))
    {
      idj <- which(warming == warm.lev[j])
      cij.use <- betai[idj, 1]
      year.use <- treat$year1[idj]
      block.use <- treat$block[idj]
      if (scale.num){
        cij.use <- scale(cij.use)
        year.use <- scale(year.use)
      }
      lmij <- lmer(cij.use ~ year.use + ((1+year.use)|block.use))
      lmijsm <- summary(lmij)
      AIC1 <- AIC(lmij)
      r2ij <- rsquared.glmm(lmij)
      lmijCS=Anova(lmij,type = "II")
      out[[j]] <- c(
        slope.fix = lmijsm$coefficients[2, 1], slope.se = lmijsm$coefficients[2, 2],
        R2M = r2ij$Marginal,
        R2C = r2ij$Conditional,
        AIC1 = AIC1,
        AIC2=r2ij$AIC,
        P.typeII=lmijCS[[3]],
        Chisq=lmijCS[[1]]
      )
    }
    outs <- Reduce(cbind, out)
    colnames(outs) <- warm.lev
    outs
  }

  tdci <- tdc.lmm(betai = betai, treat = treat,scale.num = scale.num, prefixi = prefixi)
  if (alpha == TRUE){
    output <- tdci
    output <- as.data.frame(output)
    output
  } else {
    r2.obs=as.vector(tdci[3:4,])
    aic.obs=as.vector(tdci[5:6,])
    
    # randomize time points and other the same as observed.
    year.lev=unique(treat$year1)
    year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
    trace.seq=seq(from=1,to=rand,by = 100)
    ind.rand=lapply(1:nrow(year.perm),
                    function(k)
                    {
                      if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                      out=list()
                      idi=year.perm[k,]
                      perm.treat=treat
                      perm.treat[,"year1"]=year.lev[idi[match(treat$year1,year.lev)]]
                      tdcr=tdc.lmm(betai = betai, treat = perm.treat,scale.num = scale.num)
                      out$r2=as.vector(tdcr[3:4,])
                      out$aic=as.vector(tdcr[5:6,])
                      out$ds=(-tdcr[1,1])-(-tdcr[1,2])
                      out
                    })
    r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
    aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
    EPS <- sqrt(.Machine$double.eps)
    p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
    p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
    
    p.values=rbind(matrix(p.r2,2,length(unique(treat$treat))),matrix(p.aic,2,length(unique(treat$treat))))
    rownames(p.values)=c("P.R2M","P.R2C","P.AIC1","P.AIC2")
    output=rbind(tdci,p.values)
    output = as.data.frame(output)
    output
  }
  } else {
      betai <- as.data.frame(betai)
      library(lme4)
      library(car)
      tdc.lmm <- function(betai, treat,scale.num = FALSE, prefixi = NULL) {
        warming <- treat$Treatment
        warm.lev <- unique(warming)
        out <- list()
        for (j in 1:length(warm.lev))
        {
          idj <- which(warming == warm.lev[j])
          cij.use <- betai[idj, 1]
          year.use <- treat$year1[idj]
          block.use <- treat$Block[idj]
          site.use <- treat$Site.code[idj]
          if (scale.num){
            cij.use <- scale(cij.use)
            year.use <- scale(year.use)
          }
          lmij <- lmer(cij.use ~ year.use + ((1+year.use)|block.use)+ ((1+year.use)|site.use))
          lmijsm <- summary(lmij)
          AIC1 <- AIC(lmij)
          r2ij <- rsquared.glmm(lmij)
          lmijCS=Anova(lmij,type = "II")
          out[[j]] <- c(
            slope.fix = lmijsm$coefficients[2, 1], slope.se = lmijsm$coefficients[2, 2],
            R2M = r2ij$Marginal,
            R2C = r2ij$Conditional,
            AIC1 = AIC1,
            AIC2=r2ij$AIC,
            P.typeII=lmijCS[[3]],
            Chisq=lmijCS[[1]]
          )
        }
        outs <- Reduce(cbind, out)
        colnames(outs) <- warm.lev
        outs
      }
      
      tdci <- tdc.lmm(betai = betai, treat = treat,scale.num = scale.num, prefixi = prefixi)
      if (alpha == TRUE){
        output <- tdci
        output <- as.data.frame(output)
        output
      } else {
        r2.obs=as.vector(tdci[3:4,])
        aic.obs=as.vector(tdci[5:6,])
        
        # randomize time points and other the same as observed.
        year.lev=unique(treat$year1)
        year.perm=vegan:::getPermuteMatrix(rand,length(year.lev))
        trace.seq=seq(from=1,to=rand,by = 100)
        ind.rand=lapply(1:nrow(year.perm),
                        function(k)
                        {
                          if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                          out=list()
                          idi=year.perm[k,]
                          perm.treat=treat
                          perm.treat[,"year1"]=year.lev[idi[match(treat$year1,year.lev)]]
                          tdcr=tdc.lmm(betai = betai, treat = perm.treat,scale.num = scale.num)
                          out$r2=as.vector(tdcr[3:4,])
                          out$aic=as.vector(tdcr[5:6,])
                          out$ds=(-tdcr[1,1])-(-tdcr[1,2])
                          out
                        })
        r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
        aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
        EPS <- sqrt(.Machine$double.eps)
        p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
        p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
        
        p.values=rbind(matrix(p.r2,2,length(unique(treat$treat))),matrix(p.aic,2,length(unique(treat$treat))))
        rownames(p.values)=c("P.R2M","P.R2C","P.AIC1","P.AIC2")
        output=rbind(tdci,p.values)
        output = as.data.frame(output)
        output
    }
  }
}
