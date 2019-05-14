
# Routine for estimating marginal effects on #
# strategic interactions #####################
##############################################

me_on_interaction = function(model,X){
  pars    = model$pars
  mod_se  = model$sum$std.error
  mar_effs = list()
  for(j in colnames(X)){
    X_means = matrix(rep(rbind(apply(X,2,mean)),
                         len=100*ncol(X)),
                     ncol=ncol(X),byrow = T)
    colnames(X_means) = colnames(X)
    X_means[,j] = seq(min(X[,j]),max(X[,j]),len=100)
    delta = mod.strat$pars[ncol(X)+1]/
      (exp(X_means%*%pars[1:ncol(X)])+1)
    delta.b = matrix(0,ncol=10000,nrow=100)
    for(i in 1:10000){
      err <- c()
      for(k in 1:(ncol(X)+1)){
        err[k] = rnorm(n=1,sd=mod_se[k])
      }
      delta.b[,i] = (pars[ncol(X)+1]+err[ncol(X)+1])/
        (exp(X_means%*%(pars[1:ncol(X)]+err[1:ncol(X)]))+1)
    }
    mar_effs[[j]] = data.frame(delta = delta,
                               se = apply(delta.b,1,sd),
                               var = X_means[,j],
                               term = j)
  }
  out = do.call(rbind,mar_effs)
  rownames(out) = NULL
  return(out)
}
