# updated 09/25/2013
# contains the following main function:
# initial method: SSVD.initial
# iterative thresholding method
# a routine that combines them both
# rank estimation

error.est <-
  function (x, u, v, n.boot = 100)
    # this is the function to get threshold level by bootstrap
    # n.boot: number of bootstrap
    # u and v: current estimates
    # the value of u does not matter except for the location of nonzero entries
  {
    x <- as.matrix(x)
    u <- as.matrix(u)
    v <- as.matrix(v)
    pu <- nrow(x)
    pv <- ncol(x)
    r <- ncol(u)
    # the rows that are zero everywhere
    sel.u <- apply(u == 0, 1, all)
    sel.v <- apply(v == 0, 1, all)
    # number of rows that are zero everywhere
    n.sel.u <- sum(sel.u)
    n.sel.v <- sum(sel.v)
    # the size of approximate pure noise matrix
    n.sel <- n.sel.u * n.sel.v
    print(c("n.sel", n.sel))
    if (n.sel < (log(pu * (pv - n.sel.v)) * pu * (pv - n.sel.v))) {
      warning("error.est: dense")
      return(sqrt(2 * log(pu)))
    }
    else {
      samp <- sample(n.sel, (pv - n.sel.v) * pu * n.boot, replace = TRUE)
      ans <- abs(matrix(x[sel.u, sel.v][samp], ncol = pv -
                          n.sel.v) %*% v[!sel.v, ])
      return(apply(ans, 2, function(x) {
        median(apply(matrix(x, n.boot, pu), 1, max))
      }))
    }                   
  }

get.subset <-
  function (x, method = c("theory", "method"), alpha.method = 0.05,
            alpha.theory = 1.5, sigma = NA, df = Inf)
    # to be called by SSVD.initial
  {
    if (method == "theory") {
      ans <- which(x/sigma^2/df > (1 + alpha.theory * sqrt(log(df)/df)))
    }
    if (method == "method") {
      ans <- holm.robust(x = x, alpha = alpha.method)
    }
    ans
  }

hard.thresh <-
  function (x, thr)
    # thr might be vector
  {
    if (length(thr) == 1) {
      return(hard.thresh.scalar(x, thr))
    }
    else {
      if (ncol(x) != length(thr)) {
        stop("hard.thresh: dimension do not match")
      }
      ans <- x
      for (i in 1:length(thr)) {
        ans[, i] <- hard.thresh.scalar(x[, i], thr[i])
      }
      return(ans)
    }
  }

hard.thresh.scalar <-
  function (x, thr)
    # thr is scalar
  {
    ans <- x
    ans[x < thr & x > -thr] <- 0
    ans
  }

holm.robust <-
  function (x, alpha = 0.05)
    # to be called by get.subset
    # x is a vector
    # alpha is the level of the hypothesis test
    # value: location of the outliers in the vector x
  {
    n <- length(x)
    ord <- order(x, decreasing = TRUE)
    pvalue <- 1 - pnorm(x, mean = median(x), sd = mad(x))
    alpha.adjust <- alpha/n:1
    TF <- (pvalue[ord] < alpha.adjust)
    if (TF[1] == TRUE) {
      ans <- ord[1:(which(TF == FALSE)[1] - 1)]
    }
    else {
      ans <- numeric(0)
    }
    ans
  }

huberize <-
  function (x, huber.beta = 0.95)
    # huber rho function
  {
    x.huber <- x^2
    delta <- quantile(x.huber, huber.beta)
    sel <- (x.huber > delta)
    x.huber[sel] <- 2 * sqrt(delta) * abs(x[sel]) - delta
    x.huber
  }

soft.thresh <-
  function (x, thr)
    # thr might be vector
  {
    if (length(thr) == 1) {
      return(soft.thresh.scaler(x, thr))
    }
    else {
      if (ncol(x) != length(thr)) {
        stop("soft.thresh: dimension do not match")
      }
      ans <- x
      for (i in 1:length(thr)) {
        ans[, i] <- soft.thresh.scaler(x[, i], thr[i])
      }
      return(ans)
    }
  }

soft.thresh.scaler <-
  function (x, thr)
    # thr is scaler
  {
    sign(x) * pmax(abs(x) - thr, 0)
  }

my.ssvd <-
  function (x, Sigma_u, Sigma_v, x_full = NULL,
            method = c("theory", "method"), alpha.method = 0.05,
            alpha.theory = 1.5, huber.beta = 0.95, sigma = NA, r = 1,
            gamma.u = sqrt(2), gamma.v = sqrt(2), dothres = "hard", tol = 1e-08,
            n.iter = 100, n.boot = 100, non.orth = FALSE, reps = 1, 
            init_norm = FALSE, init_type = "row_selection")
    # the main function
    # init_norm: whether to normalize the initial estimator.
  {
    if (init_type == "row_selection"){
      ans.initial <- ssvd.initial(x, Sigma_u = Sigma_u, Sigma_v = Sigma_v, 
                                  method = method, alpha.method = alpha.method,
                                  alpha.theory = alpha.theory, huber.beta = huber.beta,
                                  sigma = sigma, r = r, init_norm = init_norm,
                                  init_type = "row_selection")
    }else{
      if(is.null(x_full)){
        ans.initial <- ssvd.initial(x, Sigma_u = Sigma_u, Sigma_v = Sigma_v, 
                                    method = method, alpha.method = alpha.method,
                                    alpha.theory = alpha.theory, huber.beta = huber.beta,
                                    sigma = sigma, r = r, init_norm = init_norm,
                                    init_type = "full")
      }else{
        ans.initial <- ssvd.initial(x_full, Sigma_u = Sigma_u, Sigma_v = Sigma_v, 
                                    method = method, alpha.method = alpha.method,
                                    alpha.theory = alpha.theory, huber.beta = huber.beta,
                                    sigma = sigma, r = r, init_norm = init_norm,
                                    init_type = "full")
      }

    }

    my.ssvd.iter.thresh(x, Sigma_u, Sigma_v,
                        method = method, u.old = ans.initial$u,
                       v.old = ans.initial$v, gamma.u = gamma.u, gamma.v = gamma.v,
                     dothres = dothres, r = r, tol = tol, n.iter = n.iter,
                     n.boot = n.boot, sigma = ans.initial$sigma.hat, 
                     non.orth = non.orth, reps = reps)

  }

ssvd.initial <-
  function (x, Sigma_u = NA, Sigma_v = NA, x_full = NULL,
            method = c("theory", "method"), alpha.method = 0.05,
            alpha.theory = 1.5, huber.beta = 0.95, sigma = NA, r = 1, 
            init_norm = FALSE, init_type = "row_selection")
    # implement SSVD initial selection in both the methodology and theoretical paper
    # to use theoretical one, set method to "theory", and set alpha.theory
    # to use methodology, set method to the "method", and set alpha.method, huber.beta
    # sigma can be provided if known, otherwise will be estimated
    # r is the desired rank
  {
    x <- as.matrix(x)
    pu <- nrow(x)
    pv <- ncol(x)
    # estimate sigma
    if (is.na(sigma)) {
      sigma.hat <- mad(as.vector(x))
    }
    else {
      sigma.hat <- sigma
    }
    # huberize x^2
    if (method == "theory") {
      x.huber <- x^2
    }
    else {
      x.huber <- huberize(x, huber.beta = huber.beta)
    }
    # row and column selection
    rownorm2 <- apply(x.huber, 1, sum)
    colnorm2 <- apply(x.huber, 2, sum)
    if (init_type == "row_selection"){
      I.row <- get.subset(rownorm2, method = method, alpha.method = alpha.method,
                          alpha.theory = alpha.theory, sigma = sigma.hat, df = pv)
      print(I.row)
      I.col <- get.subset(colnorm2, method = method, alpha.method = alpha.method,
                          alpha.theory = alpha.theory, sigma = sigma.hat, df = pu)
    }else{
      if (is.null(x_full)){
        pu <- nrow(Sigma_u)
        pv <- ncol(Sigma_v)
        rownorm2 <- apply(x.huber, 1, sum)[1:pu]
        colnorm2 <- apply(x.huber, 2, sum)[(pu+1):(pu + pv)]
      }
      I.row <- (order(rownorm2, decreasing = TRUE))[1:min(0.15 * n, pu)]
      I.col <- (order(colnorm2, decreasing = TRUE))[1:min(0.15 * n, pv)]

    }

    # sanitary
    if (length(I.row) < r) {
      warning("SSVD.initial: Number of selected rows less than rank!")
      I.row <- (order(rownorm2, decreasing = TRUE))[1:min(r +
                                                            10, pu)]
    }
    if (length(I.col) < r) {
      warning("SSVD.initial: Number of selected cols less than rank!")
      I.col <- (order(colnorm2, decreasing = TRUE))[1:min(r +
                                                            10, pv)]
    }
    # CCA on selected submatrix
    Sigma_u_inv = sqrtm(Sigma_u[I.row, I.row])$Binv
    Sigma_v_inv = sqrtm(Sigma_v[I.col, I.col])$Binv
    x.sub.svd <- svd(Sigma_u_inv %*% x[I.row, I.col, drop = FALSE] %*% Sigma_v_inv, 
                     nu = r, nv = r)
    # expanding
    u.hat <- matrix(0, pu, r)
    v.hat <- matrix(0, pv, r)
    u.hat[I.row, ] <- Sigma_u_inv %*% x.sub.svd$u
    v.hat[I.col, ] <- Sigma_v_inv %*% x.sub.svd$v
    d.hat <- x.sub.svd$d[1:r]
    
    if(init_norm){
      norm_u = t(u.hat) %*% Sigma_u %*% u.hat
      index_zero = which(diag(norm_u) < 1e-8)
      for (ind in index_zero){
        norm_u[ind, ind] = 1
      }
      u.hat <- u.hat %*% (sqrtm(norm_u)$Binv)
      
      norm_v = t(v.hat) %*% Sigma_v %*% v.hat
      index_zero = which(diag(norm_v) < 1e-8)
      for (ind in index_zero){
        norm_v[ind, ind] = 1
      }
      v.hat <- v.hat %*% (sqrtm(norm_v)$Binv)
      
    }
    
    list(u = u.hat, v = v.hat, d = d.hat, sigma.hat = sigma.hat)
  }

my.ssvd.iter.thresh <-
  function (x, Sigma_u, Sigma_v,
            method = c("theory", "method"), u.old, v.old, gamma.u = sqrt(2),
            gamma.v = sqrt(2), dothres = "hard", r = ncol(u.old), tol = 1e-08,
            reps = 1,
            n.iter = 100, n.boot = 100, sigma = NA, non.orth = FALSE)
    # x is the observed matrix
    # u.old and v.old are the starting points
    # gamma.u and gamma.v are for theoretical computation:
    # threshold level is set to be gamma*sqrt(log(p))
    # n.iter: max number of iteration
    # n.boot: number of bootstrap
    # sigma: noise level
    # r: rank
    # tol: tolerance level for the convergence
    # non.orth: if set to be TRUE, then the last step does not involve orthoganalization
  {
    x <- as.matrix(x)
    pu <- nrow(x)
    pv <- ncol(x)
    puv <- max(pu, pv)
    # estimate sigma
    if (is.na(sigma)) {
      sigma.hat <- mad(as.vector(x))
    }
    else {
      sigma.hat <- sigma
    }
    # rescaling
    x.scaled <- x/sigma.hat
    # thresholding rule
    if (dothres == "hard") {
      thresh <- hard.thresh
    }
    else if (dothres == "soft") {
      thresh <- soft.thresh
    }
    else {
      warning("SSVD.iter.thresh: argument dothres not recognized! Use hard-thresholding as default")
      thresh <- hard.thresh
    }
    # initialization
    dist.u <- 1
    dist.v <- 1
    i.iter = 1
    u.cur <- u.old
    v.cur <- v.old
    u_init <- u.old
    v_init <- v.old
    
    
    
    while (i.iter <= n.iter & max(dist.u, dist.v) > tol) {
      u.old <- u.cur
      # multiplication
      sel.v <- apply(v.cur == 0, 1, all)
      
      u.cur <- x.scaled[,!sel.v, drop = FALSE] %*% v.cur[!sel.v,, drop = FALSE]
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Orthogonalize Y using economic QR decomposition: Y=QR
      #If q > 0 perfrom q subspace iterations
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if( reps > 1 ) {
        for( i in 1:reps) {
          u.cur <- qr.Q( qr(u.cur, complete = FALSE) , complete = FALSE )
          Z <- crossprod_help(x.scaled , u.cur )
          Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
          u.cur <- u.cur %*% Z
        }#End for
        remove(Z)
      }#End if
      
      #x.mul = x.scaled[,!sel.v, drop = FALSE]
      #if (reps > 1){
      #  ### increase signal
      #  svd_result <- svd(x.mul)
      #  x.mul <- svd_result$u %*% diag(svd_result$d^reps) %*% t(svd_result$v)
      #}
      
      #u.cur <- x.mul %*% v.cur[!sel.v,, drop = FALSE]
      
      # thresholding
      if (method == "theory") {
        u.cur <- thresh(u.cur, gamma.u * sqrt(log(pu)))
      }
      else {
        threshold.u <- error.est(x.scaled, u.old, v.cur,
                                 n.boot = n.boot)
        u.cur <- thresh(u.cur, threshold.u)
      }
      if (all(u.cur == 0)) {
        # every entry is zero
        warning("SSVD.iter.thresh: Overthresh: all zero!")
        u.cur <- u.old
        dist.u <- 0
        dist.v <- 0
        break
      }
      else {
        #### Normalize
        #norm_u <- apply(u.cur, 2, norm)
        #inv_norm <- sapply(norm_u, function(x){ ifelse(x > 1e-4, 1/x, 1)})
        #u.cur <- u.cur %*% diag(inv_norm) ### Sigma_X U Lambda= Sigma_XY V
        # QR decomposition
        #u.cur <- qr.Q(qr(u.cur))
        #selected_col = which(apply(u.cur^2, 2, sum) >0)
        #norm_u = t(u.cur) %*% Sigma_u %*% u.cur
        #index_zero = which(diag(norm_u) < 1e-8)
        index_non_zero = which(apply(u.cur^2, 1, sum) > 1e-6)
        # #print(index_zero)
        # #for (ind in index_zero){
        # #  norm_u[ind, ind] = 1
        # #}
        #u.cur[index_non_zero, ] <- solve(Sigma_u[index_non_zero, index_non_zero]) %*% u.cur #%*% (sqrtm(norm_u)$Binv)
        u.cur[index_non_zero,] <- solve(Sigma_u[index_non_zero, index_non_zero]) %*% u.cur[index_non_zero, ]
        norm_u = t(u.cur[index_non_zero, ]) %*% Sigma_u[index_non_zero, index_non_zero] %*% u.cur[index_non_zero, ]
        u.cur[index_non_zero, ] <- u.cur[index_non_zero, ] %*% ( sqrtm(norm_u + 1e-6 * diag(nrow(norm_u)))$Binv)
        u.cur[-index_non_zero, ] <- 0
        # #norm_u <- apply(u.cur, 2, norm)
        # #inv_norm <- sapply(norm_u, function(x){ ifelse(x > 1e-6, 1/x, 1)})
        # #u.cur <- u.cur %*% diag(inv_norm)
        # print("check norm")
        # print(t(u.cur) %*% Sigma_u %*% u.cur)
        #norm_u <- apply(u.cur, 2, norm)
        #inv_norm <- sapply(norm_u, function(x){ ifelse(x > 1e-4, 1/x, 1)})
        #u.cur <- u.cur %*% diag(inv_norm)
        # QR decomposition
        #u.cur <- qr.Q(qr(u.cur))
        #selected_col = which(apply(u.cur^2, 2, sum) >0)
        #norm_u = t(u.cur) %*% Sigma_u %*% u.cur
        #index_zero = which(diag(norm_u) < 1e-8)
        #print(index_zero)
        #for (ind in index_zero){
        #  norm_u[ind, ind] = 1
        #}
        #u.cur <- u.cur %*% (sqrtm(norm_u)$Binv)
        #dist.u <- subsp.dist.orth(u.cur, u.old)
        dist.u <- subsp.dist.orth(u.cur, u.old)
        #print(c("Dim u", dim(u.cur)))
      }
      v.old <- v.cur
      # multiplication
      sel.u <- apply(u.cur == 0, 1, all)
      
      t_x.mul = t(x.scaled[!sel.u,, drop = FALSE]) 
      if (reps > 1){
        svd_result <- svd(t_x.mul)
        x.mul <- svd_result$u %*% diag(svd_result$d^reps) %*% t(svd_result$v)
      }
      v.cur <- t_x.mul %*% u.cur[!sel.u,, drop = FALSE]
      # thresholding
      if (method == "theory") {
        v.cur <- thresh(v.cur, gamma.v * sqrt(log(pv)))
      }
      else {
        threshold.v <- error.est(t(x.scaled), v.old, u.cur,
                                 n.boot = n.boot)
        v.cur <- thresh(v.cur, threshold.v)
      }
      if (all(v.cur == 0)) {
        # every entry is zero
        warning("SSVD.iter.thresh: Overthresh: all zero!")
        v.cur <- v.old
        dist.u <- 0
        dist.v <- 0
        break
      }
      else {
        norm_v = t(v.cur) %*% Sigma_v %*% v.cur
        # #index_zero = which(diag(norm_u) < 1e-8)
        index_non_zero = which(apply(v.cur^2, 1, sum) > 1e-6)
        # #v.cur[index_non_zero, ] <- solve(Sigma_v[index_non_zero, index_non_zero]) %*% v.cur #%*% (sqrtm(norm_u)$Binv)
        v.cur[index_non_zero, ] <- solve(Sigma_v[index_non_zero, index_non_zero]) %*% v.cur[index_non_zero, ]
        norm_v = t(v.cur[index_non_zero, ]) %*% Sigma_v[index_non_zero, index_non_zero] %*% v.cur[index_non_zero, ]
        v.cur[index_non_zero, ] <- v.cur[index_non_zero, ] %*% sqrtm(norm_v)$Binv
        v.cur[-index_non_zero, ] <- 0
        # print(v.cur)
        #inv_norm <- sapply(norm_v, function(x){ ifelse(x > 1e-6, 1/x, 1)})
        #v.cur <- v.cur %*% diag(inv_norm)
        #norm_v <- apply(v.cur, 2, norm)
        #inv_norm <- sapply(norm_v, function(x){ ifelse(x > 1e-4, 1/x, 1)})
        #v.cur <- v.cur %*% diag(inv_norm)
        # QR decomposition
        #v.cur <- qr.Q(qr(v.cur))
        ### check if some rows are 0
        #norm_v = t(v.cur) %*% Sigma_v %*% v.cur
        #index_zero = which(diag(norm_v) < 1e-8)
        #for (ind in index_zero){
        #  norm_v[ind, ind] = 1
        #}
        #v.cur <- v.cur %*% (sqrtm(norm_v)$Binv)
        print("check norm v")
        print(t(v.cur) %*% Sigma_v %*% v.cur)
        
        #norm_v <- apply(v.cur, 2, norm)
        #inv_norm <- sapply(norm_v, function(x){ ifelse(x > 1e-4, 1/x, 1)})
        #v.cur <- v.cur %*% diag(inv_norm)
        # QR decomposition
        #v.cur <- qr.Q(qr(v.cur))
        ### check if some rows are 0
        #norm_v = t(v.cur) %*% Sigma_v %*% v.cur
        #index_zero = which(diag(norm_v) < 1e-8)
        #for (ind in index_zero){
        #  norm_v[ind, ind] = 1
        #}
        #v.cur <- v.cur %*% (sqrtm(norm_v)$Binv)
        dist.v <- subsp.dist.orth(v.cur, v.old)
        #print(c("Dim v", dim(v.cur)))
        #print(c("Dim u after v", dim(u.cur)))
      }
      i.iter <- i.iter + 1
      #print("iteration is")
      #print(i.iter)
      #print("dim U")
      #print(dim(u.cur))
    }
    
    if(non.orth == TRUE){
      u.old <- u.cur
      u.cur <- x.scaled %*% v.cur
      # thresholding
      if (method == "theory") {
        u.cur <- thresh(u.cur, gamma.u * sqrt(log(pu)))
      }
      else {
        threshold.u <- error.est(x.scaled, u.old, v.cur,
                                 n.boot = n.boot)
        u.cur <- thresh(u.cur, threshold.u)
      }
      u.cur <- apply(u.cur, 2, function(x){x/sqrt(sum(x^2))})
      v.old <- v.cur
      # multiplication
      v.cur <- t(x.scaled) %*% u.old
      # thresholding
      if (method == "theory") {
        v.cur <- thresh(v.cur, gamma.v * sqrt(log(pv)))
      }
      else {
        threshold.v <- error.est(t(x.scaled), v.old, u.cur,
                                 n.boot = n.boot)
        v.cur <- thresh(v.cur, threshold.v)
      }
      v.cur <- apply(v.cur, 2, function(x){x/sqrt(sum(x^2))})
    }
    d.cur <- diag(t(u.cur) %*% x %*% v.cur)
    # To make singular values positive
    u.cur <- u.cur %*% diag(sign(d.cur), r, r)
    if (i.iter == n.iter) {
      warning("increase n.iter")
    }
    if (non.orth == TRUE){
      return(list(u = u.cur, v = v.cur,
                  u.orth = u.old %*% diag(sign(apply(u.cur*u.old, 2, sum)), r, r),
                  v.orth = v.old %*% diag(sign(apply(v.cur*v.old, 2, sum)), r, r),
                  d = abs(d.cur), niter = i.iter - 1,
                  sigma.hat = sigma.hat, dist.u = dist.u, dist.v = dist.v,
                  u_init = u_init,
                  v_init = v_init))
    } else{
      return(list(u = u.cur, v = v.cur, d = abs(d.cur), niter = i.iter - 1,
                  sigma.hat = sigma.hat, dist.u = dist.u, dist.v = dist.v,
                  u_init = u_init,
                  v_init = v_init))
    }
  }

subsp.dist.orth <-
  function (A.q, B.q)
    # the squared distance between two subspaces
    # the columns of A and B are orthonormal already
  {
    overlap <- svd(t(A.q) %*% B.q, nu = 0, nv = 0)$d
    overlap <- overlap[length(overlap)]
    1 - overlap^2
}

