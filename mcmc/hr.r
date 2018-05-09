
# intersects line x+t*v with G*x<h polytope
# returns: c(X1,minT), where X1 is the intersection points, and minT is the T values at intersection point..
intersectLineWithPlanes <- function(G,h,X,V){
  
  #hGX= ( G %*% X )  - h 
  hGX = sweep( -G %*% X , 2 , h, '+',check.margin = FALSE)
  
  T = hGX / G %*% V
  
  T[T<0] = Inf
  
  minT = apply(T,2,min)
  
  X1 = sweep( X , 2 , V %*% minT , '+',check.margin = FALSE) 
  
  #ret = c(X1,minT)
  ret = vector('list',2);
  ret[[1]] = X1
  ret[[2]] = minT
  return(ret);
}


# compute hr step
hrstep_unif <- function (G,h,X) {
  
  d=ncol(G)
  m=nrow(G)
  n=ncol(X);
  
  # compute random directions:
  vr = matrix(rnorm(n*d),nrow=d)
  vr = vr / apply(vr,2,sum)
  
  # compute intersections:
  tbr = intersectLineWithPlanes(G,h,X,vr)
  tar = intersectLineWithPlanes(G,h,X,-vr)
  ta = -tar[[2]]
  tb = tbr[[2]]
  
  # compute random new points:
  r01 = matrix(runif(n))
  t01 = ta + (tb-ta) * t(r01) 
  
  Xnext = X + sweep( vr , 2 , t01 , '*')# ( t01 * vr ) 
  return(Xnext)
}

# performs hr step, using a trivial method to estimate 1d density
# ns2d is number of samples for 1d density estim
hrstep_01 <- function (G,h,f,X,ns1d) {
  d=ncol(G)
  m=nrow(G)
  n=ncol(X);
  
  # compute random directions:
  vr = matrix(rnorm(n*d),nrow=d)
  vr = vr / repmat(colSums(vr^2),d,1)
  
  # compute intersections:
  tbr = intersectLineWithPlanes(G,h,X,vr)
  tar = intersectLineWithPlanes(G,h,X,-vr)
  ixa = tar[[1]]
  ixb = tbr[[1]]
  
  # compute random new points:
  uu = runif(n)
  #loop over samples and process each separately..
  for(zi in 1:n){
    #xx = sweep( ixa[,zi,drop=FALSE] , 2 , sweep( (ixb[,zi,drop=FALSE]-ixa[,zi,drop=FALSE]),2, seq(0,1,length=(ns1d+2)) ,'*' ) , '+')
    #xx = bsxfun( '+' , ixa[,zi,drop=FALSE] , bsxfun('*',(ixb[,zi,drop=FALSE]-ixa[,zi,drop=FALSE]),seq(0,1,length=(ns1d+2))) )
    #xx = outer( (ixb[,zi,drop=FALSE]-ixa[,zi,drop=FALSE]),seq(0,1,length=(ns1d+2)) )
    xx = repmat( ixa[,zi,drop=FALSE] , 1 , ns1d+2 ) +  ( (ixb[,zi,drop=FALSE]-ixa[,zi,drop=FALSE]) %*% seq(0,1,length=(ns1d+2)) )
    ff   = apply(xx[,-c(1,ns1d+2)],2,f)
    ffcs = cumsum(ff);
    ffcs = ffcs / tail(ffcs,1)
    
    # now draw unif. rd to get next point:
    u  = uu[zi];
    # find first entry larger than u
    ux = min(which(ffcs>u))
    # par(mfrow=c(2,1)); plot(ff); plot(ffcs)
    # pos now is: center of xx[ux] andxx[ux+1]
    xn = 0.5*(xx[,ux] + xx[,ux+1]) 
    
    X[,zi] = xn
  }
  return(X)
}

# compute hr steps
# nsteps
hr_unif <- function(G,h,X,nsteps) {
  
  cs = 0
  XX = vector('list',length(nsteps))
  
  ii = 1;
  for(nsi in nsteps){
    stepsi = nsi-cs;
    for( zi in 1:stepsi ){
      X = hrstep_unif(G,h,X)
    }
    XX[[ii]] = X;
    ii = ii+1
    cs = nsi;
  }
  return(XX);
}

hr <- function(G,h,f,X,nsteps,ns1d) {
  
  cs = 0
  XX = vector('list',length(nsteps))
  
  ii = 1;
  for(nsi in nsteps){
    stepsi = nsi-cs;
    for( zi in 1:stepsi ){
      X = hrstep_01(G,h,f,X,ns1d)
    }
    XX[[ii]] = X;
    ii = ii+1
    cs = nsi;
  }
  return(XX);
}




