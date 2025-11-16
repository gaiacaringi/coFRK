#funzione che crea oggetto SRE_mv_multiresProva 
#Versione generale (non per forza tutti i processi hanno stesso numero di osservazioni) 
SRE_mv_multiresProva <- function(
    f, data, basis, BAUs, response_vars, n_levels, dims_basis_lattice, S0_list,
    est_error = TRUE,
    average_in_BAU = TRUE,
    sum_variables = NULL,
    normalise_wts = TRUE,
    fs_model = "ind", vgm_model = NULL,
    K_type = c("block-exponential", "precision", "unstructured"),
    normalise_basis = TRUE,
    response = c("gaussian", "poisson", "gamma",
                 "inverse-gaussian", "negative-binomial", "binomial"),
    link = c("identity", "log", "sqrt", "logit", "probit", "cloglog", "inverse", "inverse-squared"),
    include_fs = TRUE, fs_by_spatial_BAU = FALSE,
    ...
) {
  # Controllo input
  if (length(f) != length(response_vars)) {
    stop("Il numero di formule deve coincidere con il numero di variabili risposta")
  }
  
  p <- length(response_vars)
  all_vars <- names(data[[1]]@data)
  covariate_names <- setdiff(all_vars, response_vars)
  covariate_names <- c("fs","altitude","total_population")
  #covariate_names <- c("fs","altitude")
  F_block <- model.matrix(delete.response(terms(f[[1]])), data = BAUs@data)
  #F_block <- as.matrix(BAUs@data[, covariate_names, drop = FALSE])
  F_mat <- Matrix::kronecker(Matrix::Diagonal(p), F_block)
  
  # Liste degli oggetti che voglio salvarmi
  SRE_list        <- vector("list", p)
  Z_list <- list()
  C_Z_list <- list()
  V_eps_list <- list()
  V_xi_list <- list()
  V_xi_list_BAUs <- list()
  S_list <- list()
  S_BAU_list <- list()
  T_list <- list()
  sigma2_xi_hat   <- list()
  num_observations <- numeric(p)
  beta_list <- list()
  eta_list <- list()

  # Ciclo su ogni variabile risposta
  for (i in seq_along(response_vars)) {
    vname <- response_vars[i]
    covariate_names <- setdiff(names(data), response_vars)
    
    #all_BAU_vars   <- names(BAUs@data)
    #covariate_names <- setdiff(all_BAU_vars, c( "w", "offset", response_vars))
    
    temp <- data[[i]]  # Estrai l'oggetto Spatial della i-esima variabile
    # Chiama SRE per questa variabile
    model <- SRE(
      f = f[[i]],
      data = list(temp),
      basis = basis,
      BAUs = BAUs,
      average_in_BAU = FALSE
    )
    model <- SRE.fit(model)
    SRE_list[[i]] <- model
    Z_list[[i]]   <- model@Z 
    V_eps_list[[i]]   <- model@Ve
    V_xi_list[[i]] <- model@Vfs  
    V_xi_list_BAUs[[i]] <- model@Vfs_BAUs
    C_Z_list[[i]] <- model@Cmat 
    sigma2_xi_hat[[i]]  <- as.numeric(model@sigma2fshat)
    num_observations[i] <- length(model@Z)
    S_list[[i]] <- model@Cmat %*% model@S0
    S_BAU_list[[i]] <- model@S0
    T_list[[i]] <- model@Cmat %*% F_block
    beta_list[[i]] <- model@alphahat # dim: (1 + num covariates) x 1
    eta_list[[i]] <- model@mu_eta
  }
  
  ## Aggrego le matrici in blocchi
  #Z <- do.call(cbind, lapply(Z_list, as.numeric))
  #Z <- Matrix::Matrix(Z, sparse = FALSE)
  Z <- Matrix::Matrix(do.call(c, lapply(Z_list, as.vector)), sparse = FALSE)
  V_eps <- Matrix::bdiag(V_eps_list)
  V_xi <- Matrix::bdiag(V_xi_list)
  V_xi_BAUs <- Matrix::bdiag(V_xi_list_BAUs)
  S_mat <- Matrix::bdiag(S_list)
  S_BAU <- Matrix::bdiag(S_BAU_list)
  T_mat <- Matrix::bdiag(T_list)
  beta_hat <- do.call(rbind, beta_list)
  

  ## Inizializzazione dei parametri
  response_matrix <- do.call(cbind, lapply(Z_list, as.vector))
  obs_var_estimates <- apply(response_matrix, 2, var)
  
  sigma2_s_init <-  do.call(rbind, sigma2_xi_hat)
  sigma2_xi_init <- obs_var_estimates * 0.2
  
  nu <- rep(0.2,p)
  alpha <- alpha_init(n_levels, nu)
  
  r0 <- 0.9
  r1 <- 0.5
  
  
  
  SRE_mv_obj <- list(
    SRE_list = SRE_list,
    Z = Z,
    C_Z_list = C_Z_list,
    V_eps = V_eps,
    V_xi = V_xi,
    V_xi_BAUs = V_xi_BAUs,
    V_xi_list = V_xi_list,
    V_xi_list_BAUs = V_xi_list_BAUs,
    S0_list = S0_list,
    S = S_mat,
    S_BAU = S_BAU,
    response_vars = response_vars,
    T_mat_BAUs = F_mat,
    T_mat = T_mat,
    loglik = NULL,
    beta_hat = beta_hat, # dim: (F+1) x 1
    sigma2_xi = sigma2_xi_init,
    #tau = rep(0.2, n_levels),
    tau = 0.2,
    sigma2_s = sigma2_s_init,
    nu = nu,
    alpha = alpha,
    r0 = r0,
    r1= r1,
    mu_c = NULL,
    Sigma_c = NULL,
    dims_basis_lattice = dims_basis_lattice,
    num_observations = num_observations,
    time = 0,
    ridge_lambda =0
  )
  
  return(SRE_mv_obj)
}

alpha_init <- function(L, nu) {
  # nu: length-p vector of smoothness (one per process)
  # returns: p x L matrix with rows summing to 1
  p <- length(nu)
  lev <- seq_len(L)
  base <- outer(nu, lev, function(n, l) 2^(-2 * n * l))  # p x L
  base / rowSums(base)
}


#Function that initialized sigma2_s
init_sigma2_s_from_obj <- function(obj) {
  p <- length(obj$response_vars)
  m <- nrow(obj$C_Z)                                   # obs per process
  
  # build Q_l = B_l B_l^T and r_vec from S0_list
  r_vec  <- sapply(obj$S0_list, ncol)
  Q_list <- build_Qlist(obj$tau, r_vec)
  
  sigma2_s <- numeric(p)
  
  for (j in seq_len(p)) {
    # indices for process j
    idx  <- ((j-1)*m + 1):(j*m)
    
    # empirical variance of Z_j (one number)
    Zj   <- as.numeric(obj$Z[idx])
    varZ <- var(Zj)
    
    # average measurement-error variance for process j
    sigma2_eps_bar <- mean(diag(obj$V_eps[idx, idx]))
    
    # K0_j = blockdiag_l( alpha_{j,l} * Q_l^{-1} )
    K0_j <- Matrix::bdiag(lapply(seq_along(Q_list), function(l) {
      obj$alpha[j, l] * Matrix::solve(Q_list[[l]])
    }))
    
    # S block for process j
    Rj  <- ncol(K0_j)
    S_j <- obj$S[idx, ((j-1)*Rj + 1):(j*Rj)]
    
    # denominator: (1/m) * tr(S K0 S^T)  -- computed without forming dense m×m
    den <- mean(rowSums((S_j %*% K0_j) * S_j))
    
    # numerator: Var(Z_j) - sigma2_xi_j - avg sigma2_eps_j
    num <- varZ - obj$sigma2_xi[j] - sigma2_eps_bar
    
    sigma2_s[j] <- max(num / den, 0)
  }
  
  obj$sigma2_s <- sigma2_s
  return(obj)
}

