### Update for sigma2_s -----

## Computes the inverse of the equicorrelation matrix
inv_equicorr <- function(rho, p) {
  a <- 1 / (1 - rho)
  b <- - rho / ((1 - rho) * (1 - rho + p * rho))
  A <- diag(a, p)
  A[] <- A + b  
  return(A)
}

## Builds G_ell^P from S_eta (process-major) and BBt_list 
# G_l^P[i,j] = tr( (B_l B_l^T) * S_{eta,l}^{(i,j)} )
# S_{eta,l}^{(i,j)} is the r_l x r_l block for (process i, process j) at level l
build_G_list_P <- function(S_eta, BBt_list, p) {
  L   <- length(BBt_list)
  Rl  <- vapply(BBt_list, nrow, integer(1))
  Rtot <- sum(Rl)
  proc_start   <- (0:(p - 1)) * Rtot
  level_offset <- c(0L, cumsum(Rl))[1:L]
  
  G_list <- vector("list", L)
  for (l in seq_len(L)) {
    BBt_l <- BBt_list[[l]]
    r_l <- nrow(BBt_l)
    G_l <- matrix(0, nrow = p, ncol = p)
    
    for (i in seq_len(p)) {
      rows_l_i <- (proc_start[i] + level_offset[l] + 1L):(proc_start[i] + level_offset[l] + r_l)
      for (j in seq_len(p)) {
        cols_l_j <- (proc_start[j] + level_offset[l] + 1L):(proc_start[j] + level_offset[l] + r_l)
        Sij <- as.matrix(S_eta[rows_l_i, cols_l_j, drop = FALSE])
        G_l[i, j] <- sum(BBt_l * Sij)
      }
    }
    G_list[[l]] <- G_l
  }
  G_list
}

## Computes sigma2_s objective to maximize 
# alpha_mat: p x L with alpha[i,l] >= 0
# r0, r1 define rho_l = r0 * exp(-r1 * (l-1))
# BBt_list: list of B_l B_l^T
# G_list:   list of G_l^P (p x p)
Q_sigma_obj <- function(sigma2_s, alpha_mat, r0, r1, BBt_list, G_list, clip_rho = 1e-8, lambda_ridge = 0) {

  p <- length(sigma2_s)
  L <- length(BBt_list)
  Rl <- vapply(BBt_list, nrow, integer(1))
  

  eps <- 1e-12
  sigma2_s <- pmax(as.numeric(sigma2_s), eps)
  
  total <- 0.0
  for (l in seq_len(L)) {
    r_l <- Rl[l]
    d_l <- sqrt(alpha_mat[, l] * sigma2_s)      
    d_l <- pmax(d_l, sqrt(eps))
    Dinv <- diag(1 / d_l, p)
    
    rho_l <- r0 * exp(-r1 * (l - 1))
    rho_min <- -1/(p - 1) + clip_rho
    rho_max <-  1 - clip_rho
    rho_l <- max(min(rho_l, rho_max), rho_min)
    C_l <- inv_equicorr(rho_l, p)  
    
    H_l <- Dinv %*% G_list[[l]] %*% Dinv
    
    total <- total + (-2 * r_l * sum(log(d_l)) - sum(diag(C_l %*% H_l)))
  }
  
  if (lambda_ridge > 0) {
    total <- total - (lambda_ridge / 2) * sum(sigma2_s^2)
  }
  
  return(total)  
}

## Updates sigma2_s 
update_Q_parameters1 <- function(obj, control = list(maxit = 60)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))

  r_vec    <- sapply(obj$S0_list, ncol)
  BBt_list <- build_Qlist(obj$tau, r_vec, obj$dims_basis_lattice)
  L        <- length(BBt_list)
  

  p <- length(obj$sigma2_s)
  G_list <- build_G_list_P(S_eta, BBt_list, p)
  
  alpha_mat <- obj$alpha
  if (!is.matrix(alpha_mat) || any(dim(alpha_mat) != c(p, L))) {
    stop("obj$alpha must be a p x L matrix of level weights (alpha[i,l]).")
  }
  

  start <- pmax(as.numeric(obj$sigma2_s), 1e-8)
  
  fn <- function(x) {
    -Q_sigma_obj(x, alpha_mat, obj$r0, obj$r1, BBt_list, G_list, obj$ridge_lambda)
  }

  lower <- rep(1e-8, p)
  upper <- rep(1e+6, p)
  
  opt <- optim(
    par     = start,
    fn      = fn,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 60), control)
  )
  
  if (!is.null(opt$par) && all(is.finite(opt$par))) {
    obj$sigma2_s <- opt$par
  }
  obj$opt_value       <- if (is.finite(opt$value)) -opt$value else NA_real_
  obj$opt_convergence <- opt$convergence
  obj$opt_message     <- if (!is.null(opt$message)) opt$message else ""
  
  return(obj)
}

#### Update for tau ------

tau_objective_modello <- function(tau, obj) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))
  r_vec <- obj$dims_basis_lattice$n_x * obj$dims_basis_lattice$n_y
  L       <- length(r_vec)
  SigInvL <- build_Sigma_list(obj$sigma2_s, obj$r0, obj$r1, obj$alpha, L)  # p×p per level

  p            <- length(obj$sigma2_s)
  Rl           <- r_vec
  Rtot         <- sum(Rl)
  proc_start   <- (0:(p - 1)) * Rtot
  level_offset <- c(0L, cumsum(Rl))[1:L]
  
  val <- 0.0
  
  for (ell in seq_len(L)) {
    tau_ell <- exp(tau * ell)
    B_ell   <- build_B(obj$dims_basis_lattice$n_x[ell],
                       obj$dims_basis_lattice$n_y[ell],
                       tau_ell)
    Rchol_B <- Matrix::Cholesky(B_ell, LDL = FALSE, super = TRUE)
    logdetB <- sum(log(Matrix::diag(Rchol_B)))  
    val <- val + 2 * p * logdetB               

    BBt <- B_ell %*% Matrix::t(B_ell)
    
    for (i in seq_len(p)) {
      rows_i_start <- proc_start[i] + level_offset[ell]
      rows_i       <- (rows_i_start + 1L):(rows_i_start + Rl[ell])
      for (j in seq_len(p)) {
        cols_j_start <- proc_start[j] + level_offset[ell]
        cols_j       <- (cols_j_start + 1L):(cols_j_start + Rl[ell])
        
        S_ij_ell <- S_eta[rows_i, cols_j, drop = FALSE]
        tr_BBt_S <- sum(BBt * S_ij_ell)
        val <- val - SigInvL[[ell]][i, j] * tr_BBt_S
      }
    }
  }
  
  as.numeric(val)
}

tau_objective_modello <- function(tau, obj) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))
  r_vec <- obj$dims_basis_lattice$n_x * obj$dims_basis_lattice$n_y
  L       <- length(r_vec)
  SigInvL <- build_Sigma_list(obj$sigma2_s, obj$r0, obj$r1, obj$alpha, L)  # p×p per level
  
  p            <- length(obj$sigma2_s)
  Rl           <- r_vec
  Rtot         <- sum(Rl)
  proc_start   <- (0:(p - 1)) * Rtot
  level_offset <- c(0L, cumsum(Rl))[1:L]
  
  val <- 0.0
  
  for (ell in seq_len(L)) {
    tau_ell <- exp(tau * ell)
    B_ell   <- build_B(obj$dims_basis_lattice$n_x[ell],
                       obj$dims_basis_lattice$n_y[ell],
                       tau_ell)
    Rchol_B <- Matrix::Cholesky(B_ell, LDL = FALSE, super = TRUE)
    logdetB <- sum(log(Matrix::diag(Rchol_B))) 
    val <- val + 2 * p * logdetB                
    

    BBt <- B_ell %*% Matrix::t(B_ell)
    
    acc <- 0.0
    for (i in seq_len(p)) {
      rows_i_start <- proc_start[i] + level_offset[ell]
      rows_i       <- (rows_i_start + 1L):(rows_i_start + Rl[ell])
      for (j in seq_len(p)) {
        cols_j_start <- proc_start[j] + level_offset[ell]
        cols_j       <- (cols_j_start + 1L):(cols_j_start + Rl[ell])
        
        S_ij_ell <- S_eta[rows_i, cols_j, drop = FALSE]
        tr_BBt_S <- sum(BBt * S_ij_ell)
        acc <- acc + SigInvL[[ell]][i, j] * tr_BBt_S
      }
    }
    val <- val -acc
  }
  
  as.numeric(val) 
}


tau_objective_pen <- function(tau, obj, tau0 = 0.2, s = 0.5) {
  f_raw <- tau_objective_modello(tau, obj)  
  pen   <- (tau - tau0)^2 / (2 * s^2)        
  f_raw - pen
}


update_Q_parameters2 <- function(obj, lower = 0.05, upper = 3) {
  opt <- optimize(function(t) -tau_objective_pen(t, obj),
                  interval = c(lower, upper))
  obj$tau <- opt$minimum
  obj
}


### Update for r0 and r1 --------
.precompute_ABT <- function(obj) {
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))
  
  r_vec    <- obj$dims_basis_lattice$n_x * obj$dims_basis_lattice$n_y
  BBt_list <- build_Qlist(obj$tau, r_vec, obj$dims_basis_lattice)
  L        <- length(r_vec)
  p        <- length(obj$sigma2_s)
  
  Rl           <- r_vec
  Rtot         <- sum(Rl)
  proc_start   <- (0:(p - 1)) * Rtot
  level_offset <- c(0L, cumsum(Rl))[1:L]
  
  A <- numeric(L)
  B <- numeric(L)

  sigma2 <- obj$sigma2_s
  
  for (ell in seq_len(L)) {
    E  <- BBt_list[[ell]]
    T_ell <- matrix(0.0, p, p)
    for (i in seq_len(p)) {
      ri0 <- proc_start[i] + level_offset[ell]
      rows_i <- (ri0 + 1L):(ri0 + Rl[ell])
      for (j in seq_len(p)) {
        cj0 <- proc_start[j] + level_offset[ell]
        cols_j <- (cj0 + 1L):(cj0 + Rl[ell])
        S_ij_ell <- S_eta[rows_i, cols_j, drop = FALSE]
        T_ell[i, j] <- sum(E * S_ij_ell)
      }
    }
    
    alpha_l <- obj$alpha[, ell]
    M_diag  <- 1 / (sigma2 * alpha_l)                         
    u_l     <- 1 / (sqrt(sigma2) * sqrt(alpha_l))
    

    A[ell] <- sum(M_diag * diag(T_ell))
    B[ell] <- as.numeric(crossprod(u_l, T_ell %*% u_l))
  }
  
  list(A = A, B = B, Rl = Rl, L = L, p = p)
}

update_Q_parameters <- function(obj, control = list(maxit = 80)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  cache <- .precompute_ABT(obj)
  A <- cache$A; B <- cache$B; Rl <- cache$Rl; L <- cache$L; p <- cache$p
  

  logdetR <- function(rho) (p - 1) * log1p(-rho) + log1p((p - 1) * rho)   # = log|R(rho)|
  a_fun   <- function(rho) 1 / (1 - rho)
  b_fun   <- function(rho) - rho / ((1 - rho) * (1 + (p - 1) * rho))
  

  f_obj <- function(par) {
    r0 <- par[1]
    r1 <- par[2]                
    lidx <- 0:(L - 1)
    rho  <- r0 * exp(-r1 * lidx)
    
    if (any(rho <= -1/(p - 1) + 1e-10) || any(rho >= 1 - 1e-10)) return(-Inf)

    val_det  <- - Rl * ( (p - 1) * log1p(-rho) + log1p((p - 1) * rho) )
    val_tr   <- - a_fun(rho) * A - b_fun(rho) * B
    
    sum(val_det + val_tr)
  }
  

  rho_min <- -1 / (p - 1) + 1e-6
  lower   <- c(rho_min, 0.0)
  upper   <- c(0.95,     10.0)
  
  opt <- optim(
    par     = c(obj$r0, obj$r1),
    fn      = function(par) -f_obj(par),    # minimizer -> negate
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 80), control)
  )
  
  obj$r0 <- opt$par[1]
  obj$r1 <- opt$par[2]
  obj
}

