Soft <-
  function(a,b){
    if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
    return(sign(a)*pmax(0,abs(a)-b))
  }

updatePi <-
  function(B,sqB,A,H,Gamma,nu,rho,Pi,tau){
    C <- Pi + 1/tau*A-nu/tau*B%*%Pi%*%B+nu/tau*sqB%*%(H-Gamma)%*%sqB
    D <- rho/tau
    return(Soft(C,D))
  }

updateH <-
  function(sqB,Gamma,nu,Pi,K){
    
    temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
    temp <- (temp+t(temp))/2
    svdtemp <- eigen(temp)
    d <- svdtemp$values
    p <- length(d)
    if(sum(pmin(1,pmax(d,0)))<=K){
      dfinal <- pmin(1,pmax(d,0))
      return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
    }
    fr <- function(x){
      sum(pmin(1,pmax(d-x,0)))
    }
    # Vincent Vu Fantope Projection
    knots <- unique(c((d-1),d))
    knots <- sort(knots,decreasing=TRUE)
    temp <- which(sapply(knots,fr)<=K)
    lentemp <- tail(temp,1)
    a=knots[lentemp]
    b=knots[lentemp+1]
    fa <- sum(pmin(pmax(d-a,0),1))
    fb <- sum(pmin(pmax(d-b,0),1))
    theta <- a+ (b-a)*(K-fa)/(fb-fa)
    dfinal <- pmin(1,pmax(d-theta,0))
    res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
    return(res)
  }

hard <-
  function(U, k, r){
    if(r>1){
      truncate.value <- sort(apply(U, 1, FUN = function(x) sum(x^2)),decreasing=TRUE)[k]
      U[which(apply(U, 1, FUN = function(x) sum(x^2))<truncate.value)] <- 0
    }else{
      truncate.value <- sort(abs(U),decreasing=TRUE)[k]
      U[which(abs(U)<truncate.value)] <- 0
    }
    return(U)
  }

