library(gss)


# uilities
centering <- function(x){
  return(x-mean(x))
}

weight.nodes <- function(pt,xdomain, scale = 1){
  qpt <- (pt-xdomain[1])/(xdomain[2]-xdomain[1])
  qpt.minus <- c(max(0,min(qpt)-1e-6),qpt[1:(length(qpt)-1)])
  wt = scale*(qpt-qpt.minus)/sum(qpt-qpt.minus)
  return(wt)
}

partK <- function(x1,x2)
{# x1 and x2 vectors of same length
  # return value = k2(x1)*k2(x2)-k4(abs(x1 - x2)).
  val <- ((x1-0.5)^2-1/12)/2*((x2-0.5)^2-1/12)/2
  wk <- abs(x1-x2)
  val <- val-((wk-0.5)^4-(wk-0.5)^2/2+7/240)/24
}

# functional concurrent model with functional response.
# adapt to irregular sampling
# assume the domains of x and y to be [0,1].
# Inputs:
#  yqlist - list of vectors with yqlist[[i]][j]=y_i(t_j), t_j are quadrature nodes.
#  xqlist - list of vectors with xqlist[[i]][j]=x_i(t_j), t_j are quadrature nodes.
#  quad - quadrature to be used for integration (not necessary in some situations)
#  method - similar to method in ssanova.
#  alpha - alpha used in the GCV score for smoothing parameter selection
#  id.basis - index vector for selected "knots" (x's).
#  xdomain - default is [0,1].
flr2d <- function (yqlist,xqlist,quad,method="v",alpha=1.4,id.basis=NULL,xdomain=c(0,1))
{
  # need to center x and y first
  #xqlist.center <- lapply(xqlist,centering)
  #yqlist.center <- lapply(yqlist,centering)
  #xmean = lapply(xqlist,mean)
  #ymean = lapply(yqlist,mean)
  
  xmat <- matrix(unlist(xqlist),ncol = length(xqlist))
  ymat <- matrix(unlist(yqlist),ncol = length(xqlist))
  xmat <- xmat - apply(xmat, 1, mean)
  ymat <- ymat - apply(ymat, 1, mean)
  xmat <- as.vector(xmat)
  ymat <- as.vector(ymat)
  
  n = length(xmat)
  if (is.null(id.basis)) id.basis <- 1:n
  nbasis <- length(id.basis)
  # qpt: design points 
  # wvec: weights for design points
  # Wxmat
  # Wymat
  qpt <- (quad$pt-xdomain[1])/(xdomain[2]-xdomain[1]) 
  wvec <- sqrt(quad$wt) # weights for design points
  Wxmat <- wvec*xmat
  Wymat <- wvec*ymat
  
  # need to delete the rows with x = 0
  #select.indx = which(Wxmat !=0)
  #if (length(select.indx) ==0) stop("values of the covariate are all zero!")
  #x <- Wxmat[select.indx]
  #y <- Wymat[select.indx]
  x <- Wxmat
  y <- Wymat
  Kmat <- outer(qpt,qpt[id.basis],partK)
  
  # s matrix Unpenalized terms evaluated at data points. X*S
  # r matrix Basis of penalized terms evaluated at data points. X*Q
  # q matrix Penalty matrix. Q
  s = cbind(1,qpt)*x
  q = Kmat[id.basis,]
  r = Kmat*x
  
  wt <- NULL
  if (qr(s)$rank < dim(s)[2])
    stop("gss error in fcl: fixed effects are linearly dependent")
  z <- sspreg1(s,r,q,y,wt,method,alpha,varht=1,random=NULL)
  c(list(id.basis=id.basis,alpha=alpha,xdomain=xdomain,
         x=x,y=y,Wxmat =Wxmat,Wymat= Wymat,quad=quad),z)
}

sspreg1 <- function(s,r,q,y,wt,method,alpha,varht,random)
{
  qr.trace <- FALSE
  if ((alpha<0)&(method%in%c("u","v"))) qr.trace <- TRUE
  alpha <- abs(alpha)
  ## get dimensions
  nobs <- nrow(r)
  nxi <- ncol(r)
  if (!is.null(s)) {
    if (is.vector(s)) nnull <- 1
    else nnull <- ncol(s)
  }
  else nnull <- 0
  if (!is.null(random)) nz <- ncol(as.matrix(random$z))
  else nz <- 0
  nxiz <- nxi + nz
  nn <- nxiz + nnull
  if (!is.null(wt)) {
    y <- wt*y
    s <- wt*s
    r <- wt*r
    if (!is.null(random)) random$z <- wt*random$z
  }
  ## cv function
  cv <- function(lambda) {
    if (is.null(random)) q.wk <- 10^(lambda+theta)*q
    else {
      q.wk <- matrix(0,nxiz,nxiz)
      q.wk[1:nxi,1:nxi] <- 10^(lambda[1]+theta)*q
      q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
        10^(2*ran.scal)*random$sigma$fun(lambda[-1],random$sigma$env)
    }
    if (qr.trace) {
      suppressWarnings(qq.wk <- chol(q.wk,pivot=TRUE))
      sr <- cbind(s,10^theta*r[,attr(qq.wk,"pivot")])
      sr <- rbind(sr,cbind(matrix(0,nxiz,nnull),qq.wk))
      sr <- qr(sr,tol=0)
      rss <- mean(qr.resid(sr,c(y,rep(0,nxiz)))[1:nobs]^2)
      trc <- sum(qr.Q(sr)[1:nobs,]^2)/nobs
      if (method=="u") score <- rss + alpha*2*varht*trc
      if (method=="v") score <- rss/(1-alpha*trc)^2
      alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
      alpha.wk <- min(alpha.wk,3)
      if (alpha.wk>alpha) {
        if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*trc
        if (method=="v") score <- rss/(1-alpha.wk*trc)^2
      }
      if (return.fit) {
        z <- .Fortran("reg",
                      as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                      as.double(q.wk), as.integer(nxiz), as.double(y),
                      as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                      as.double(alpha), varht=as.double(varht),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                      wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
        z$score <- score
        assign("fit",z[c(1:5,7)],inherits=TRUE)
      }
    }
    else {
      z <- .Fortran("reg",
                    as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                    as.double(q.wk), as.integer(nxiz), as.double(y),
                    as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                    as.double(alpha), varht=as.double(varht),
                    score=double(1), dc=double(nn),
                    as.double(.Machine$double.eps),
                    chol=double(nn*nn), double(nn),
                    jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                    wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                    PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
      if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
      assign("fit",z[c(1:5,7)],inherits=TRUE)
      score <- z$score
      alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
      alpha.wk <- min(alpha.wk,3)
      if (alpha.wk>alpha) {
        if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*z$wk[2]
        if (method=="v") score <- z$wk[1]/(1-alpha.wk*z$wk[2])^2
      }
    }
    score
  }
  cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
  ## initialization
  tmp <- sum(r^2)
  if (is.null(s)) theta <- 0
  else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
  log.la0 <- log10(tmp/sum(diag(q))) + theta
  if (!is.null(random)) {
    ran.scal <- theta - log10(sum(random$z^2)/nz/tmp*nxi) / 2
    r <- cbind(r,10^(ran.scal-theta)*random$z)
  }
  else ran.scal <- NULL
  ## lambda search
  return.fit <- FALSE
  fit <- NULL
  if (is.null(random)) la <- log.la0
  else la <- c(log.la0,random$init)
  if (length(la)-1) {
    counter <- 0
    ## scale and shift cv
    tmp <- abs(cv(la))
    cv.scale <- 1
    cv.shift <- 0
    if (tmp<1&tmp>10^(-4)) {
      cv.scale <- 10/tmp
      cv.shift <- 0
    }
    if (tmp<10^(-4)) {
      cv.scale <- 10^2
      cv.shift <- 10
    }
    repeat {
      zz <- nlm(cv.wk,la,stepmax=1,ndigit=7)
      if (zz$code<=3) break
      la <- zz$est
      counter <- counter + 1
      if (counter>=5) {
        warning("gss warning in ssanova: iteration for model selection fails to converge")
        break
      }
    }
  }
  else {
    mn0 <- log.la0-6
    mx0 <- log.la0+6
    repeat {
      mn <- max(la-1,mn0)
      mx <- min(la+1,mx0)
      zz <- nlm0(cv,c(mn,mx))
      if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
          (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
      else la <- zz$est
    }
  }
  ## return
  return.fit <- TRUE
  jk1 <- cv(zz$est)
  if (is.null(random)) q.wk <- 10^theta*q
  else {
    q.wk <- matrix(0,nxiz,nxiz)
    q.wk[1:nxi,1:nxi] <- 10^theta*q
    q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
      10^(2*ran.scal-zz$est[1])*random$sigma$fun(zz$est[-1],random$sigma$env)
  }
  se.aux <- regaux(s,10^theta*r,q.wk,zz$est[1],fit)
  c <- fit$dc[nnull+(1:nxi)]
  if (nnull) d <- fit$dc[1:nnull]
  else d <- NULL
  if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
  else b <- NULL
  c(list(method=method,theta=theta,ran.scal=ran.scal,c=c,d=d,b=b,
         nlambda=zz$est[1],zeta=zz$est[-1]),fit[-3],list(se.aux=se.aux))
}
## Auxiliary Quantities for Standard Error Calculation
regaux <- function(s,r,q,nlambda,fit)
{
  nnull <- dim(s)[2]
  nn <- nnull +  dim(q)[1]
  zzz <- eigen(q,symmetric=TRUE)
  rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
  val <- zzz$val[1:rkq]
  vec <- zzz$vec[,1:rkq,drop=FALSE]
  if (nnull) {
    wk1 <- qr(s)
    wk1 <- (qr.qty(wk1,r%*%vec))[-(1:nnull),]
  }
  else wk1 <- r%*%vec
  wk2 <- t(t(wk1)/sqrt(val))
  wk2 <- t(wk2)%*%wk2
  wk2 <- solve(wk2+diag(10^nlambda,dim(wk2)[1]),wk2)
  wk2 <- (wk2+t(wk2))/2
  wk2 <- t(wk2/sqrt(val))/sqrt(val)
  wk2 <- diag(1/val,dim(wk2)[1])-wk2
  z <- .Fortran("regaux",
                as.double(fit$chol), as.integer(nn),
                as.integer(fit$jpvt), as.integer(fit$rkv),
                drcr=as.double(t(cbind(s,r))%*%r%*%vec), as.integer(rkq),
                sms=double(nnull^2), as.integer(nnull), double(nn*nnull),
                PACKAGE="gss")[c("drcr","sms")]
  drcr <- matrix(z$drcr,nn,rkq)
  dr <- drcr[1:nnull,,drop=FALSE]
  sms <- 10^nlambda*matrix(z$sms,nnull,nnull)
  wk1 <- matrix(0,nnull+rkq,nnull+rkq)
  wk1[1:nnull,1:nnull] <- sms
  wk1[1:nnull,nnull+(1:rkq)] <- -t(t(dr)/val)
  wk1[nnull+(1:rkq),nnull+(1:rkq)] <- wk2
  suppressWarnings(z <- chol(wk1,pivot=TRUE))
  wk1 <- z
  rkw <- attr(z,"rank")
  while (wk1[rkw,rkw]<wk1[1,1]*sqrt(.Machine$double.eps)) rkw <- rkw-1
  wk1[row(wk1)>col(wk1)] <- 0
  if (rkw<nnull+rkq)
    wk1[(rkw+1):(nnull+rkq),(rkw+1):(nnull+rkq)] <- diag(0,nnull+rkq-rkw)
  hfac <- wk1
  hfac[,attr(z,"pivot")] <- wk1
  list(vec=vec,hfac=hfac)
}
# function for evaluating coefficient function estimate from 2D
# functional linear regression at given grid points.
# Output: an nt*ns matrix
estimate.flr2d <- function(object,tgrid,se.fit = FALSE)
{
  nnull <- length(object$d)
  nt <- length(tgrid)
  quad <- object$quad
  xdomain <- object$xdomain
  qpt <- (quad$pt-xdomain[1])/(xdomain[2]-xdomain[1])
  twk <- (tgrid-xdomain[1])/(xdomain[2]-xdomain[1])

  Ktmat <- outer(twk,qpt[object$id.basis],partK)
  
  tgrid.stand <- (tgrid-xdomain[1])/(xdomain[2]-xdomain[1])
  btest <- object$d[1]+object$d[2]*(tgrid.stand)
  btestrk <- Ktmat%*%(10^object$theta*object$c)
  btest <- btest+btestrk
  
  if (se.fit) {
    b <- object$varht/10^object$nlambda
    pse <- NULL
    
    ## Compute posterior variance
    ss <- rbind(rep(1,nt),tgrid.stand)
    r.wk <- 10^object$theta*Ktmat
    rr <- t(r.wk%*%object$se.aux$vec)
    wk <- object$se.aux$hfac%*%rbind(ss,rr)
    pse <- sqrt(b*apply(wk^2,2,sum))
  }
  
  
  if (se.fit) list(fit = btest, se.fit = pse)
  else btest

}