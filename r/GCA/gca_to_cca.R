gca_to_cca <-
  function(a_estimate, S, pp){
    p1 = pp[1];
    p2 = pp[2];
    p = p1 + p2;
    nnz_indices = apply(a_estimate, 1, norm) 
    nnz_indices_x = nnz_indices[which(nnz_indices<(p1+1))]
    nnz_indices_y = nnz_indices[which(nnz_indices>(p1))]
    u_estimate = a_estimate[1:p1,]
    v_estimate = a_estimate[(p1+1):(p1+p2),]
    if (length(which(nnz_indices<(p1+1)))==0){
      sigmaxhat = S[nnz_indices_x,nnz_indices_x];
      u_estimate[nnz_indices_x,] = a_estimate[nnz_indices_x,] %*% pracma::sqrtm(t(a_estimate[nnz_indices_x,]) %*% sigmaxhat %*% a_estimate[nnz_indices_x,])$Binv;
    }
    if (length(which(nnz_indices>(p1)))==0){
      sigmayhat = S[nnz_indices_y,nnz_indices_y];
      v_estimate[nnz_indices_y,] = a_estimate[nnz_indices_y,] %*% pracma::sqrtm(t(a_estimate[nnz_indices_y,]) %*% sigmayhat %*% a_estimate[nnz_indices_y,])$Binv;
    }
    l = list("u" = u_estimate, "v" = v_estimate)
    # sigmaxhat = S[1:p1,1:p1];
    # sigmayhat = S[(p1+1):p,(p1+1):p];
    # u_estimate = a_estimate[1:p1,] %*% pracma::sqrtm(t(a_estimate[1:p1,]) %*% sigmaxhat %*% a_estimate[1:p1,])$Binv;
    # v_estimate = a_estimate[(p1+1):p,] %*% pracma::sqrtm(t(a_estimate[(p1+1):p,]) %*% sigmayhat %*% a_estimate[(p1+1):p,])$Binv;
    # l = list("u" = u_estimate, "v" = v_estimate)
    return(l)
  }