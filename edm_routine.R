
# Routine for Estimating EDM #
##############################


edm_fun = function(X,y,R,A,b){
  beta = exp(X%*%b[1:ncol(X)])
  gamma = (b[ncol(X)+1])
  delta1 = beta/(beta + 1)
  delta2 = gamma/(beta + 1)
  y_hat = delta1*R - delta2*A
  rss   = sum((y_hat-y)^2)
  return(rss)
}
edm_pred = function(X,y,R,A,b,CI=FALSE,alpha=0.95){
  beta = exp(X%*%b[1:ncol(X)])
  gamma = (b[ncol(X)+1])
  delta1 = beta/(beta + 1)
  delta2 = gamma/(beta + 1)
  y_hat = delta1*R - delta2*A
  if(CI==FALSE){
    return(y_hat)
  } else {
    se.fit = sqrt(rowSums((cbind(X,A)%*%mod.strat$vcov) * 
                            cbind(X,A)))
    Qt = c(-1,1)*qt((1 - alpha)/2,
                    (nrow(X)+1)-(ncol(X)+1),
                    lower.tail=F)
    out = data.frame(prediction=y_hat,
                     lo = y_hat + Qt[1]*se.fit,
                     hi = y_hat + Qt[2]*se.fit)
    return(out)
  }
}
edm = function(X,y,R,A){
  out = optim(edm_fun,
              par=c(rep(0,len=1+ncol(X))),
              hessian = TRUE, method = "BFGS", 
              control = list(REPORT = 10, 
                             trace = 1, 
                             maxit = 50000),
              X=X,y=y,R=R,A=A)
  vcov = try(as.matrix(solve(out$hessian, 
                             tol=1e-24)), T)
  sum = data.frame(term = c(colnames(X),"gamma"),
                   estimate=out$par,
                   std.error = sqrt(diag(vcov)),
                   statistic = out$par/
                     sqrt(diag(vcov)),
                   p.value = round(2*pnorm(abs(out$par/
                                                 sqrt(diag(vcov))),
                                           lower.tail=FALSE),4))
  fit = edm_pred(X=X,R=R,A=A,
                    b=out$par)
  return(list(pars=out$par,
              vcov=vcov,
              rss=out$value,
              fit=fit,
              sum=sum))
}

