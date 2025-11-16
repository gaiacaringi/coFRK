SRE_mv.predict_multiresGenerale <- function(obj, C_P) {

  p <- length(obj$response_vars)
  n <- obj$num_observations
  T_mat <- obj$T_mat  
  T_mat_BAUs <- obj$T_mat_BAUs
  S <- obj$S_BAU # B x R
  beta_hat  <- obj$beta_hat  
  mu_c  <- obj$mu_c      
  Sigma_c <- obj$Sigma_c 
  V_xi_blocks <- lapply(seq_along(obj$V_xi_list), function(j) {
    (obj$V_xi_list_BAUs)[[j]] * (obj$sigma2_xi)[[j]]
  })
  V_xi <- Matrix::bdiag(V_xi_blocks)
  V_xi <- V_xi + Diagonal(nrow(V_xi)) * 1e-6 #ho visto che quale V_xi ha degli autovalori nulli e questo dovrebbe toglierli, perche dopo mi serve che sia SDP
  V_eps <- obj$V_eps
  Cz_list <- obj$C_Z_list
  Z <- obj$Z
  V_eps_inv <- Matrix::solve(V_eps)
  
  
  #Matrice di mapping
  C_P_big <- Matrix::kronecker(Matrix::Diagonal(p), C_P) 
  Cz <- Matrix::bdiag(Cz_list)
  
  #Matrice Π = [S | I]
  n_total <- nrow(S)  # numero totale di BAUs
  I_mat <- Matrix::Diagonal(n_total)
  Pi <- cbind(S, I_mat)
  
  # Q precision matrix 
  BBt_list <- build_Qlist(obj$tau, sapply(obj$S0_list, ncol),obj$dims_basis_lattice)
  Q <- build_Q_total(BBt_list = BBt_list, obj$sigma2_s, obj$r0, obj$r1, obj$alpha, length(BBt_list))
  
  #Q <- Matrix::forceSymmetric(Q, uplo = "U")
  
  #V_xi inverse
  Vxi_inv_blocks <- lapply(seq_along(obj$V_xi_list_BAUs), function(j) {   ## <<
    Bj <- (obj$V_xi_list_BAUs)[[j]]
    s2 <- max(obj$sigma2_xi[[j]], 1e-8)
    if (inherits(Bj, "diagonalMatrix")) {
      Matrix::Diagonal(x = 1 / (s2 * Matrix::diag(Bj)))
    } else {
      # generic (symmetric) inverse with tiny jitter for safety
      Bj_ <- Matrix::forceSymmetric(Bj, uplo = "U")
      inv_block <- try(Matrix::solve(Bj_), silent = TRUE)
      if (inherits(inv_block, "try-error")) {
        inv_block <- Matrix::solve(Bj_ + 1e-8 * Matrix::Diagonal(nrow(Bj_)))
      }
      (1 / s2) * inv_block
    }
  })
  
  
  
  Vxi_inv <- Matrix::bdiag(Vxi_inv_blocks)
  #Matrice W = (c, ξ)
  # Λ = blockdiag(K, V_ξ)
  Lambda_inv <- Matrix::bdiag(Q, Vxi_inv)  
  
  #Sigma_W
  print(dim(Pi))
  print(dim(Cz))
  print(dim(V_eps))
  tPi_SigZinv <- Matrix::t(Pi) %*% Matrix::t(Cz) %*% V_eps_inv %*% Cz
  Sigma_W <- Matrix::solve(tPi_SigZinv %*% Pi + Lambda_inv)
  
  #mu_W
  resid <- Z - T_mat %*% beta_hat
  mu_W <- Sigma_W %*%  Matrix::t(Pi) %*% Matrix::t(Cz) %*% V_eps_inv %*% resid
  
  #E[Y|Z]
  Y_hat <- T_mat_BAUs %*% beta_hat + Pi %*% mu_W
  
  #Var(Y|Z)
  Var_Y <- Pi %*% Sigma_W %*% Matrix::t(Pi)
  
  #Aggrego la predzione sul dominio predizionale indicato da C_P
  print(dim(C_P_big))
  print(dim(Y_hat))
  Y_P_mean <- C_P_big %*% Y_hat
  Y_P_var  <- C_P_big %*% Var_Y %*% Matrix::t(C_P_big)
  
  return(list(mean = Y_P_mean, var = Y_P_var))
  
}

