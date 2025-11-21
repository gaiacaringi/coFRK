### helper functions ----
is_spd <- function(M) {
  # Forza simmetria se non già
  if (!Matrix::isSymmetric(M)) return(FALSE)
  
  # Prova a calcolare autovalori
  eigenvals <- tryCatch(
    eigen(as.matrix(M), symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) return(NA)
  )
  
  if (any(is.na(eigenvals))) return(FALSE)
  
  return(all(eigenvals > 0))
}

### -----
SRE_mv.fit <- function(SRE_mv_obj, n_EM = 100L, tol = 0.01, print_lik = FALSE) {
  start_time <- Sys.time()
  
  llk <- rep(NA, n_EM)  #vettore che tiene traccia del valore della log-verosimiglianza ad ogni iterazione
  Z <- SRE_mv_obj$Z                      #osservazioni (tutti i processi incolonnati nello stesso vettore)
  p <- length(SRE_mv_obj$response_vars)  # numero di processi
  n_obs <- SRE_mv_obj$num_observations                     # numero di osservazioni (assumo ogni osservazione "completa",i.e. con tutti i processi)
  
  for (iter in 1:n_EM) {
    
    # E-step
    SRE_mv_obj <- E_step_mv(SRE_mv_obj)
    assign("SRE_snapshot", SRE_mv_obj, envir = .GlobalEnv)
    
    # M-step
    SRE_mv_obj <- M_step_mv(SRE_mv_obj)
    assign("SRE_snapshot", SRE_mv_obj, envir = .GlobalEnv)
    
    # Log-likelihood
    llk[iter] <- logLik_mv(SRE_mv_obj)
    
    if (print_lik) cat("Iter", iter, ": logLik =", llk[iter], "\n")
    
    if (iter > 1 && abs(llk[iter] - llk[iter - 1]) < tol) {
      message("Convergenza raggiunta dopo ", iter, " iterazioni.")
      break
    }
  }
  
  end_time <- Sys.time()     
  elapsed <- end_time - start_time
  SRE_mv_obj$time <- as.numeric(elapsed)
  
  return(SRE_mv_obj)
}

### E step ------
E_step_mv <- function(obj) {
  
  Z <- obj$Z
  T_mat <- obj$T_mat   # blocco diagonale con T1,...,Tp
  beta <- obj$beta_hat
  S <- obj$S       # (I_p ⊗ S0)
  V_xi_blocks <- lapply(seq_along(obj$V_xi_list), function(j) {
    obj$sigma2_xi[j] * obj$V_xi_list[[j]]
  })
  V_xi <- Matrix::bdiag(V_xi_blocks)
  D <- obj$V_eps + V_xi       # block-diagonal (errore intra-BAU xi  + measurement eps)
  L <- length(obj$S0_list)
  
  BBt_list <- build_Qlist(obj$tau, sapply(obj$S0_list, ncol), obj$dims_basis_lattice)
  Q <- build_Q_total(BBt_list = BBt_list, obj$sigma2_s, obj$r0, obj$r1, obj$alpha, L)
  
  # Inversione D diretta
  #D_inv <- Matrix::solve(D)
  # Inversione D versione ottimizzata 
  chol_D <- Matrix::Cholesky(D)
  D_inv <- Matrix::solve(chol_D)  
  dim(Q)
  
  # E-step: calcolare Sigma_c e mu_c
  Sigma_c_inv <- Matrix::t(S) %*% D_inv %*% S + Q
  Sigma_c <- Matrix::solve(Sigma_c_inv) #varianza condizionata 
  resid <- Z - T_mat %*% beta
  mu_c <- Sigma_c %*% Matrix::t(S) %*% D_inv %*% resid #media condizionata
  obj$mu_c <- mu_c
  obj$Sigma_c <- Sigma_c
  cat("E step fatto")
  
  return(obj)
}



# Function to build B_l SAR matrix(q is the number of basis in level l: q x q grid  ->  B is (q^2 x q^2) )
build_B <- function(n_x,n_y, tau) {
  n <- n_x * n_y
  a  <- 4 + tau^2
  ids <- matrix(seq_len(n), nrow = n_y, ncol = n_x, byrow = FALSE)
  
  if (n_x > 1) {
    eH <- cbind(as.vector(ids[, 1:(n_x - 1), drop = FALSE]),
                as.vector(ids[, 2:n_x,       drop = FALSE]))
  } else eH <- cbind(integer(0), integer(0))
  if (n_y > 1) {
    eV <- cbind(as.vector(ids[1:(n_y - 1), , drop = FALSE]),
                as.vector(ids[2:n_y,       , drop = FALSE]))
  } else eV <- cbind(integer(0), integer(0))
  edges <- rbind(eH, eV)
  if (nrow(edges) == 0) {
    return(a * Diagonal(n))  # degenerate case: single node
  }
  O <- sparseMatrix(
    i = c(edges[,1], edges[,2]),
    j = c(edges[,2], edges[,1]),
    x = -1,
    dims = c(n, n)
  )
  
  B <- a * Diagonal(n) + O
  return(B)
}

#Function to build Q0 = (B*Bt)
#Multilevel version: list of Q0 of length L 
#Input: tau (scalar), r (vector of length L -> number of basis per level )
library(Matrix)
build_Qlist <- function(tau, r_vec, dims_basis_lattice) {
  L <- length(r_vec)
  Q0_list <-  vector("list", L)
  for (l in seq_len(L)) {
    tau_l <- exp(tau*l)
    B <- build_B(dims_basis_lattice$n_x[l],dims_basis_lattice$n_y[l],tau_l)
    BBt <- B %*% Matrix::t(B)
    
    Q0_list[[l]] <- BBt
  }
  
  return(Q0_list)
}


#Function to build Sigma_l inverse, l=1...L
#Extension multilevel: list of length L of Sigma_l
# NB: function return inverse of Sigma_l !
#Input: 
# -sigma2_s (vector of length p )
# -r0 (scalar), r1 (scalar)
# -L (scalar): numbers of levels
build_Sigma_list <- function(sigma2_s, r0, r1, alpha, L,
                                jitter0 = 1e-8, max_attempt = 5, clip_eps = 1e-10) {
  p <- length(sigma2_s)
  rho_lo <- -1/(p - 1) + clip_eps   
  rho_hi <-  1 - clip_eps
  if (p < 1L) stop("sigma2_s deve avere almeno un elemento")
  if (L < 1L) stop("L deve essere >= 1")
  
  s <- sqrt(sigma2_s)
  SigmaInv_list <- vector("list", L)
  

  for (l in 1:L) {
    rho_l <- r0 * exp(-r1 * (l-1))
    rho_l <- max(min(rho_l, rho_hi), rho_lo)   
    # |rho_l| < 1
    if (abs(rho_l) >= 1) {
      rho_l <- sign(rho_l) * (1 - clip_eps)
      warning(sprintf("rho_l clippato a %.12f al livello l=%d per garantire SPD.", rho_l, l))
    }
   
    alpha_l <- alpha[,l]
    
    Sigma <- matrix(0, nrow = p, ncol = p)
    diag(Sigma) <- sigma2_s * alpha_l
    
    if (p > 1) {
      for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
          sij <- s[i] * s[j]
          alphaij <- sqrt(alpha[i,l]*alpha[j,l])
          Sigma[i, j] <- rho_l * sij *alphaij
          Sigma[j, i] <- Sigma[i, j]
        }
      }
      Sigma <- 0.5 * (Sigma + t(Sigma))
    }
    
    
    attempt <- 0
    while (!is_spd(Sigma) && attempt < max_attempt) {
      Sigma <- Sigma + diag(jitter0 * 10^attempt, p)
      attempt <- attempt + 1
    }
    if (!is_spd(Sigma)) stop(sprintf("Sigma_l (l=%d) non è SPD anche dopo la stabilizzazione.", l))
    
    
    Rchol <- chol(Sigma)                
    SigmaInv <- chol2inv(Rchol)         
    
  
    SigmaInv <- 0.5 * (SigmaInv + t(SigmaInv))
    
    SigmaInv_list[[l]] <- SigmaInv
  }

  return(SigmaInv_list)
}


## Function to build Q total 
#Input: BBt_list, r0,r1,L
build_Q_total <- function(BBt_list, sigma2_s, r0, r1, alpha, L) {
  SigmaInv_list <- build_Sigma_list(sigma2_s, r0, r1, alpha, L)
  stopifnot(length(BBt_list) == length(SigmaInv_list))
  L <- length(BBt_list)
  if (L < 1L) stop("Liste vuote.")
  
  Rl <- vapply(BBt_list, nrow, integer(1))
  if (any(vapply(BBt_list, ncol, integer(1)) != Rl))
    stop("Ogni BBt_l deve essere quadrata.")
  Rtot <- sum(Rl)
  
  p <- nrow(SigmaInv_list[[1]])
  if (any(vapply(SigmaInv_list, ncol, integer(1)) != p) ||
      any(vapply(SigmaInv_list, nrow, integer(1)) != p))
    stop("Ogni SigmaInv_l deve essere p x p, stesso p per tutti i livelli.")
  
  # offsets: ordinamento per processi, dentro ogni processo concateno i livelli
  # c = (c1, c2, ..., cp), con c_j = (c_{j1}, ..., c_{jL})
  proc_start   <- (0:(p - 1)) * Rtot                   # inizio blocco di ogni processo
  level_offset <- c(0L, cumsum(Rl))[1:L]               # inizio di ogni livello dentro al blocco di processo
  
  # matrice finale sparsa inizializzata a zero
  Q <- Matrix::Matrix(0, nrow = p * Rtot, ncol = p * Rtot, sparse = TRUE)
  
  # funzione di comodo per aggiungere un blocco su (rows, cols)
  add_block <- function(Q, rows, cols, Block) {
    Block <- methods::as(Block, "dgCMatrix")
    Q[rows, cols] <- Q[rows, cols] + Block
    Q
  }
  
  for (l in seq_len(L)) {
    BBt_l <- BBt_list[[l]]
    if (!inherits(BBt_l, "sparseMatrix")) BBt_l <- Matrix::Matrix(BBt_l, sparse = TRUE)
    
    n_l <- nrow(BBt_l) 
    if (!is.finite(n_l) || n_l <= 0L) stop(sprintf("BBt_list[[%d]] ha dimensione non valida.", l))
    
    for (i in seq_len(p)) {
      # righe del sotto-vettore (processo i, livello l)
      row_start <- proc_start[i] + level_offset[l]
      rows_l    <- (row_start + 1L):(row_start + Rl[l])
      
      # blocco diagonale (i,i): aggiungo SigmaInv[ii] * BBt_l
      coef_ii <- SigmaInv_list[[l]][i, i]
      if (coef_ii != 0) {
        Q <- add_block(Q, rows_l, rows_l, coef_ii * BBt_l)
      }
      
      # blocchi tra processi (i,j) con j > i: metto blocco e anche il simmetrico
      if (p > 1L) {
        for (dj in seq_len(max(0L, p - i))) {
          j <- i +dj
          col_start <- proc_start[j] + level_offset[l]
          cols_l    <- (col_start + 1L):(col_start + n_l)
          
          coef_ij <- SigmaInv_list[[l]][i, j]
          if (coef_ij != 0) {
            Block <- coef_ij * BBt_l
            # parte superiore
            Q <- add_block(Q, rows_l, cols_l, Block)
            # parte inferiore (simmetrica): stesso valore perché BBt_l è simmetrica
            Q <- add_block(Q, cols_l, rows_l, Block)
          }
        }
      }
    }
  }
  
  # forza simmetria numerica (non strettamente necessario perché abbiamo riempito entrambi i lati)
  Q <- Matrix::forceSymmetric(Q, uplo = "U")
  return(Q)
}


### M step ------
M_step_mv <- function(obj) {
  Z <- obj$Z
  T_mat <- obj$T_mat
  S <- obj$S
  D <- obj$V_eps + obj$V_xi
  mu_c <- obj$mu_c
  Sigma_c <- obj$Sigma_c
  
  #Aggiorno beta
  obj$beta_hat <- update_beta(Z, T_mat, S, mu_c, D)
  
  #Inizializzazione parametri per ottimizzazione
  init_params <- list(
    log_sigma2 = log(obj$sigma2_s),
    rho = obj$rho,
    tau = obj$tau
  )
  p <- length(obj$response_vars)
  r <- ncol(obj$S0)
  
  #Aggiorno Q
  obj <- update_Q_parameters1(obj)
  obj <- update_Q_parameters2(obj)
  #obj <- update_Q_parameters3(obj)
  obj <- update_Q_parameters(obj)
  
  #Aggioro sigma_xi 
  #sigma2_xi_opt<- update_sigma2_xi(obj)
  obj <- update_sigma2_xi(obj)
  
  #Salvo stime
  
  #obj$optim_info <- list(value = K_opt$value, convergence = K_opt$convergence)
  #obj$sigma2_xi <- sigma2_xi_opt 
  
  
  return(obj)
}

#Funzione per aggiornare beta (effetti fissi)
update_beta <- function(Z, T_mat, S, mu_c, D) {
  resid_eta <- Z - S %*% mu_c 
  D_inv <- Matrix::solve(D)
  XtDX <- Matrix::t(T_mat) %*% D_inv %*% T_mat
  XtDy <- Matrix::t(T_mat) %*% D_inv %*% resid_eta
  beta_new <- Matrix::solve(XtDX, XtDy)
  return(beta_new)
}

#Funzione per agigonrare i parametri di Q
#Function to update parameters of Q: r0, r1
#It updaes the parameters in obj
update_Q_parameters <- function(obj, control = list(maxit = 50)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  # S_eta = E[c c^T | Z] = Sigma_c + mu mu^T
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))
  
  BBt_list <- build_Qlist(obj$tau, sapply(obj$S0_list, ncol), obj$dims_basis_lattice)
  #Q <- build_Q_total(BBt_list = BBt_list, obj$sigma2_s, r0, r1, obj$alpha, L)
  L <- length(BBt_list)
  
  r0 <- obj$r0
  r1 <- obj$r1
  tau <- obj$tau
  
  # objective: f(params) = log|Q| - tr(Q S_eta), to MAXIMIZE
  obj_fun <- function(par) {
    r0 <- par[1]; r1 <- par[2]; 
    #BBt_list <- build_Qlist(tau, sapply(obj$S0_list, ncol), obj$dims_basis_lattice)
    Q <- build_Q_total(BBt_list = BBt_list, obj$sigma2_s, r0, r1, obj$alpha, L)
    
    R  <- Matrix::chol(Q)
    ld <- 2 * sum(log(Matrix::diag(R)))
    
    # tr(Q S_eta) = somma Hadamard
    trQS <- sum(Q * S_eta)
    
    ld - trQS
    
  }
  
  # constaints: r0 in (-0.95, 0.95), r1 >= 0
  p <- length(obj$sigma2_s)
  rho_min <-  -1/(p - 1) + 1e-6
  lower <- c(rho_min, 0.0)   
  upper <- c( 0.95, 10.0)
  
  opt <- optim(
    par    = c(obj$r0, obj$r1),
    function(par) -obj_fun(par),  # minimiser -> negate
    method = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 60), control)
  )
  print(opt$par[1])
  # write back the parameters, in the same order we packed them
  obj$r0 <- opt$par[1]
  obj$r1 <- opt$par[2]
  
  return(obj)
}

#Funzione per agigonrare i parametri di Q
#Function to update parameters of Q: r0, r1, nu
#It updaes the parameters in obj
update_Q_parameters <- function(obj, control = list(maxit = 50)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  # S_eta = E[c c^T | Z] = Sigma_c + mu mu^T
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))
  
  BBt_list <- build_Qlist(obj$tau, sapply(obj$S0_list, ncol), obj$dims_basis_lattice)
  #Q <- build_Q_total(BBt_list = BBt_list, obj$sigma2_s, r0, r1, obj$alpha, L)
  L <- length(BBt_list)
  
  r0 <- obj$r0
  r1 <- obj$r1
  tau <- obj$tau
  
  # objective: f(params) = log|Q| - tr(Q S_eta), to MAXIMIZE
  obj_fun <- function(par) {
    r0 <- par[1]; r1 <- par[2]; 
    nu <- exp(par[-c(1,2)])
    alpha_now <- alpha_init(L, nu)
    #BBt_list <- build_Qlist(tau, sapply(obj$S0_list, ncol), obj$dims_basis_lattice)
    Q <- build_Q_total(BBt_list = BBt_list, obj$sigma2_s, r0, r1, alpha_now, L)
    
    R  <- Matrix::chol(Q)
    ld <- 2 * sum(log(Matrix::diag(R)))
    
    # tr(Q S_eta) = somma Hadamard
    trQS <- sum(Q * S_eta)
    
    ld - trQS
    
  }
  
  # constaints: r0 in (-0.95, 0.95), r1 >= 0
  p <- length(obj$sigma2_s)
  rho_min <-  -1/(p - 1) + 1e-6
  lower <- c(rho_min, 0.0)   
  upper <- c( 0.95, 10.0)
  
  opt <- optim(
    par    = c(obj$r0, obj$r1, log(pmax(obj$nu, 1e-3))),
    function(par) -obj_fun(par),  # minimiser -> negate
    method = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 60), control)
  )
  print(opt$par[1])
  # write back the parameters, in the same order we packed them
  obj$nu <- exp(opt$par[-c(1,2)])
  obj$r0 <- opt$par[1]
  obj$r1 <- opt$par[2]
  obj$alpha <- alpha_init(L, obj$nu)
  
  return(obj)
}

# Funzione che aggiorna i parametri: tau
update_Q_parameters2 <- function(obj, control = list(maxit = 60)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c)) # S_eta = E[cc^T | Z]
  
  r_vec <- sapply(obj$S0_list, ncol)
  L     <- length(r_vec)
  
  # start: tau per livello (positivo) -> ottimizziamo su log(tau)
  #if (length(obj$tau) != L) stop("obj$tau deve avere lunghezza L.")
  theta0 <- as.numeric(log(obj$tau))
  
  # obiettivo: f(tau) = log|Q(tau)| - tr(Q(tau) * S_eta)  (massimizzare)
  # -> minimizziamo il negativo
  obj_fun <- function(theta) {
    tau <- exp(theta) # garantisce positività
    
    # ricostruisci B_l B_l^T a questi tau
    BBt_list <- build_Qlist(tau, r_vec,obj$dims_basis_lattice)
    
    # costruisci Q(tau) con sigma2_s, r0, r1 fissi
    Q <- build_Q_total(
      BBt_list = BBt_list,
      sigma2_s = obj$sigma2_s,
      r0 = obj$r0, r1 = obj$r1, obj$alpha, L = L
    )
    
    # log|Q| via chol con penalizzazione se non SPD
    R <- try(Matrix::chol(Q), silent = TRUE)
    if (inherits(R, "try-error")) return(1e15)  # penalizza
    
    ld <- 2 * sum(log(Matrix::diag(R)))
    trQS <- sum(Q * S_eta)  # tr(Q S_eta) via Hadamard
    
    # minimizzatore -> negativo dell'obbiettivo da massimizzare
    -(ld - trQS)
  }
  
  # bounds su log(tau) SE VETTORE 
  #lower <- rep(log(0.05), L)
  #upper <- rep(log(10.0), L)
  
  #bounds su log(tau) SE PARAMETRO UNICO 
  lower <- log(0.05)
  upper <- log(10.0)
  
  
  opt <- optim(
    par     = theta0,
    fn      = obj_fun,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 60), control)
  )
  
  # aggiorna tau solo se ha senso
  if (!is.null(opt$par) && all(is.finite(opt$par))) {
    obj$tau <- exp(opt$par)
  }
  
  # info utili
  obj$opt_value        <- if (is.finite(opt$value)) -(opt$value) else NA_real_
  obj$opt_convergence  <- opt$convergence  # 0 = ok
  obj$opt_message      <- if (!is.null(opt$message)) opt$message else ""
  
  return(obj)
}

pc_logprior_theta <- function(theta, u = 1.0, alpha = 0.1) {
  # P(sigma_s > u) = alpha  -> lambda = -log(alpha)/u
  lambda <- -log(alpha) / u
  sigma_s <- exp(theta/2)          # since theta = log(sigma^2_s)
  # log prior up to additive constant:
  return( - lambda * sigma_s )
}

ridge_logprior_theta <- function(theta, mu = log(1.0), tau = 1.0) {
  
  return( - (theta - mu)^2 / (2 * tau^2) )
}

update_Q_parameters1 <- function(obj, control = list(maxit = 60)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))  # E[cc^T | Z]
  r_vec    <- sapply(obj$S0_list, ncol)
  BBt_list <- build_Qlist(obj$tau, r_vec, obj$dims_basis_lattice)
  L        <- length(BBt_list)
  
  p <- length(obj$sigma2_s)
  if (p < 1L) stop("obj$sigma2_s deve avere lunghezza >= 1")
  theta0 <- log(as.numeric(obj$sigma2_s))  # theta = log sigma2_s
  
  # objective: maximize  f = log|Q| - tr(Q S_eta) + log prior
  # -> our fn returns NEGATIVE of that (minimizer)
  obj_fun <- function(theta) {
    sigma2_s <- exp(theta)
    Q <- build_Q_total(
      BBt_list = BBt_list,
      sigma2_s = sigma2_s,
      r0 = obj$r0, r1 = obj$r1, obj$alpha, L = L
    )
    R <- try(Matrix::chol(Q), silent = TRUE)
    if (inherits(R, "try-error")) return(1e15)
    ld   <- 2 * sum(log(Matrix::diag(R)))
    trQS <- sum(Q * S_eta)
    
    # --- pick ONE prior: PC (A) OR Ridge (B) ---
    # (A) PC prior on sigma_s (SD): set u, alpha for gentle shrinkage
    #logprior <- pc_logprior_theta(theta, u = 1.0, alpha = 0.1)
    # (B) Ridge on log sigma2_s (comment the previous line if you prefer ridge)
     logprior <- ridge_logprior_theta(theta, mu = log(0.7^2), tau = 0.4) 
    
    # return negative penalized objective for minimizer
    return( - (ld - trQS + logprior) )
  }
  
  lower <- rep(log(1e-6), p)
  upper <- rep(log(1e+2), p)
  
  opt <- optim(
    par     = theta0,
    fn      = obj_fun,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 60), control)
  )
  
  if (!is.null(opt$par) && all(is.finite(opt$par))) {
    obj$sigma2_s <- exp(opt$par)
  }
  obj$opt_value       <- if (is.finite(opt$value)) -(opt$value) else NA_real_
  obj$opt_convergence <- opt$convergence
  obj$opt_message     <- if (!is.null(opt$message)) opt$message else ""
  return(obj)
}


# Funzione che aggiorna i parametri sigma2_s
update_Q_parameters1 <- function(obj, control = list(maxit = 60)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  # S_eta = E[c c^T | Z]
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))
  
  # B B^T per livello (dipende da tau ma NON da sigma2_s)
  r_vec    <- sapply(obj$S0_list, ncol)
  BBt_list <- build_Qlist(obj$tau, r_vec, obj$dims_basis_lattice)
  L        <- length(BBt_list)
  
  # start sui log (per imporre sigma2_s > 0)
  p <- length(obj$sigma2_s)
  if (p < 1L) stop("obj$sigma2_s deve avere lunghezza >= 1")
  theta0 <- log(as.numeric(obj$sigma2_s))
  
  # objective: f(sigma2_s) = log|Q| - tr(Q S_eta)  (MAXIMIZE)
  # -> minimizziamo il negativo
  obj_fun <- function(theta) {
    sigma2_s <- exp(theta)
    
    # Q dipende da sigma2_s, r0, r1, (tau è già in BBt_list)
    Q <- build_Q_total(
      BBt_list = BBt_list,
      sigma2_s = sigma2_s,
      r0 = obj$r0, r1 = obj$r1, obj$alpha, L = L
    )
    
    # log|Q| via Cholesky; penalizza se non SPD
    R <- try(Matrix::chol(Q), silent = TRUE)
    if (inherits(R, "try-error")) return(1e15)
    
    ld   <- 2 * sum(log(Matrix::diag(R)))
    trQS <- sum(Q * S_eta)   # tr(Q S_eta) con prodotto di Hadamard
    
    -(ld - trQS)  # minimizer
  }
  
  # bounds sui log(sigma2_s)
  lower <- rep(log(1e-6), p)
  upper <- rep(log(1e+2), p)
  
  opt <- optim(
    par     = theta0,
    fn      = obj_fun,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 60), control)
  )
  
  # aggiorna se tutto ok
  if (!is.null(opt$par) && all(is.finite(opt$par))) {
    obj$sigma2_s <- exp(opt$par)
  }
  obj$opt_value       <- if (is.finite(opt$value)) -(opt$value) else NA_real_
  obj$opt_convergence <- opt$convergence   # 0 = ok
  obj$opt_message     <- if (!is.null(opt$message)) opt$message else ""
  
  return(obj)
}


# Funzione che aggiorna i parametri nu
update_Q_parameters3 <- function(obj, control = list(maxit = 60)) {
  stopifnot(!is.null(obj$mu_c), !is.null(obj$Sigma_c))
  
  # S_eta = E[c c^T | Z]
  S_eta <- obj$Sigma_c + tcrossprod(as.vector(obj$mu_c))
  
  # B B^T per livello (dipende da tau ma NON da sigma2_s)
  r_vec    <- sapply(obj$S0_list, ncol)
  BBt_list <- build_Qlist(obj$tau, r_vec, obj$dims_basis_lattice)
  L        <- length(BBt_list)
  
  # start sui log (per imporre sigma2_s > 0)
  p <- length(obj$nu)
  if (p < 1L) stop("obj$nu deve avere lunghezza >= 1")
  theta0 <- log(pmax(as.numeric(obj$nu), 1e-6)) 
  
  # objective: f(sigma2_s) = log|Q| - tr(Q S_eta)  (MAXIMIZE)
  # -> minimizziamo il negativo
  obj_fun <- function(theta) {
    nu <- pmax(exp(theta), 1e-6)
    alpha_now <- alpha_init(L, nu)
    # Q dipende da sigma2_s, r0, r1, (tau è già in BBt_list)
    Q <- build_Q_total(
      BBt_list = BBt_list,
      sigma2_s = sigma2_s,
      r0 = obj$r0, r1 = obj$r1, alpha_now, L = L
    )
    
    # log|Q| via Cholesky; penalizza se non SPD
    R <- try(Matrix::chol(Q), silent = TRUE)
    if (inherits(R, "try-error")) return(1e15)
    
    ld   <- 2 * sum(log(Matrix::diag(R)))
    trQS <- sum(Q * S_eta)   # tr(Q S_eta) con prodotto di Hadamard
    pen <- 1e-3 * sum((nu - 0.5)^2)
    -(ld - trQS -pen )  # minimizer
  }
  
  # bounds sui log(sigma2_s)
  lower <- rep(log(1e-3), p)
  upper <- rep(log(10.0), p)
  
  opt <- optim(
    par     = theta0,
    fn      = obj_fun,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = modifyList(list(maxit = 60), control)
  )
  
  # aggiorna se tutto ok
  if (!is.null(opt$par) && all(is.finite(opt$par))) {
    obj$nu    <- pmax(exp(opt$par), 1e-6)
    obj$alpha <- alpha_init(L, obj$nu)
  }
  obj$opt_value       <- if (is.finite(opt$value)) -(opt$value) else NA_real_
  obj$opt_convergence <- opt$convergence   # 0 = ok
  obj$opt_message     <- if (!is.null(opt$message)) opt$message else ""
  
  return(obj)
}


# Funzione che aggiorna i parametri sigma2_xi -----
trace_Omega <- function(S, Sigma_c, mu_c, Z, T_mat, beta) {
  mu_c <- as.numeric(mu_c)
  resid <- as.numeric(Z - T_mat %*% beta)
  
  #Omega is the sum of 4 parts -> we compute the trace as the sum of the traces of these 4 components
  A <- sum( (S %*% Sigma_c) * S ) 
  Smu <- as.numeric(S %*% mu_c)       
  B <- sum(Smu * Smu)
  C <- -2 * sum(resid * Smu)
  D <- sum(resid * resid)
  
  return( A + B + C + D)
}

# ----- trace_Omega per variable (returns a length-p vector) -----
trace_Omega_blocks <- function(obj, m_per_var) {
  mu_c <- obj$mu_c
  Sigma_c <- obj$Sigma_c
  beta <- obj$beta_hat
  Z <- obj$Z
  T_mat <- obj$T_mat
  S <- obj$S
  resid <- as.numeric(Z - T_mat %*% beta)
  
  idx_ends   <- cumsum(m_per_var)
  idx_starts <- c(1, head(idx_ends, -1) + 1)
  
  p <- length(m_per_var)
  tr_vec <- numeric(p)

  for (i in seq_len(p)) {
    idx  <- idx_starts[i]:idx_ends[i]
    S_i  <- S[idx, , drop = FALSE]
    r_i  <- resid[idx]
    
    # Components for tr(Ω_ii):
    # A_i = tr(S_i Sigma_c S_i^T)
    A_i <- sum((S_i %*% Sigma_c) * S_i)     # efficient trace trick
    
    # B_i = || S_i mu_c ||^2
    Smu_i <- as.numeric(S_i %*% mu_c)
    B_i <- sum(Smu_i * Smu_i)
    
    # C_i = -2 * r_i^T (S_i mu_c)
    C_i <- -2 * sum(r_i * Smu_i)
    
    # D_i = || r_i ||^2
    D_i <- sum(r_i * r_i)
    
    tr_vec[i] <- A_i + B_i + C_i + D_i
  }
  
  tr_vec
}

# ----- multivariate update for sigma2_xi (length p) -----
update_sigma2_xi <- function(obj) {
  p        <- length(obj$response_vars)
  Z        <- obj$Z
  T_mat    <- obj$T_mat
  beta     <- obj$beta_hat
  S        <- obj$S
  mu_c     <- obj$mu_c
  Sigma_c  <- obj$Sigma_c
  V_eps    <- obj$V_eps
  
  
  m_per_var <- rep(length(Z)/p,p)
  gamma1_vec <- rep(1,p)
  dV <- diag(V_eps)
  idx_ends   <- cumsum(m_per_var)
  idx_starts <- c(1, head(idx_ends, -1) + 1)
  gamma2_vec <- vapply(seq_len(p), function(i) {
    idx <- idx_starts[i]:idx_ends[i]
    mean(dV[idx])
  }, numeric(1))
  
  tr_Omega_vec <- trace_Omega_blocks(obj, m_per_var)
  s2_new <- (tr_Omega_vec / m_per_var - gamma2_vec) / gamma1_vec
  s2_new <- pmax(s2_new, 0)
  
  obj$sigma2_xi <- as.numeric(s2_new)
  obj
}

#funzione di update pe rsigma2_xi che funziona per p=1
update_sigma2_xi <- function(obj){
  p <- length(obj$response_vars)
  Z <- obj$Z
  T_mat <- obj$T_mat
  beta <- obj$beta_hat
  S <- obj$S
  mu_c <- obj$mu_c
  Sigma_c <- obj$Sigma_c
  C_Z <- obj$C_Z
  gamma1 <- 1
  gamma2 <- obj$V_eps[1,1]
  m <- dim(obj$V_eps)[1]
  
  tr_Omega <- trace_Omega(S,Sigma_c,mu_c,Z,T_mat,beta)
  s2_new <- (tr_Omega / m -  gamma2) / (gamma1)
  obj$sigma2_xi <- s2_new
  return(obj)
}





### Likelihood (for convergence) -------
logLik_mv <- function(obj) {
  Z     <- obj$Z
  T_mat <- obj$T_mat
  S     <- obj$S
  V_xi_blocks <- lapply(seq_along(obj$V_xi_list), function(j) {
    obj$sigma2_xi[j] * obj$V_xi_list[[j]]
  })
  V_xi <- Matrix::bdiag(V_xi_blocks)
  D     <- obj$V_eps + V_xi
  beta  <- obj$beta_hat
  BBt_list <- build_Qlist(obj$tau, sapply(obj$S0_list, ncol), obj$dims_basis_lattice)
  L <- length(BBt_list)
  Q <- build_Q_total(BBt_list = BBt_list, obj$sigma2_s, obj$r0, obj$r1, obj$alpha, L)
  
  resid <- Z - T_mat %*% beta
  out   <- compute_Vinv_logdet_quad(S, D, Q, resid)
  
  n <- length(Z)
  llk <- -0.5 * (n * log(2*pi) + out$logdetV + out$quad)
  as.numeric(llk)
}

# Computation of:
#   - V^{-1} = D^{-1} - D^{-1} S (Q + S^T D^{-1} S)^{-1} S^T D^{-1}
#   - log|V|  = log|D| - log|Q| + log|Q + S^T D^{-1} S|
# and the quadratic form (Z - Tβ)^T V^{-1} (Z - Tβ)
#

compute_Vinv_logdet_quad <- function(S, D, Q, resid) {
  cholD <- try(Matrix::Cholesky(D, LDL = FALSE, super = TRUE), silent = TRUE)
  if (inherits(cholD, "try-error")) {
    D <- D + 1e-8 * Matrix::Diagonal(nrow(D))  # tiny jitter if needed
    cholD <- Matrix::Cholesky(D, LDL = FALSE, super = TRUE)
  }
  
  A        <- Matrix::solve(cholD, S)            # A = D^{-1} S  (solve D * X = S)
  StDinvS  <- Matrix::crossprod(S, A)            # S^T D^{-1} S
  
  # G = Q + S^T D^{-1} S   
  G    <- Matrix::forceSymmetric(Q + StDinvS, uplo = "U")
  cholG <- try(Matrix::Cholesky(G, LDL = FALSE, super = TRUE), silent = TRUE)
  if (inherits(cholG, "try-error")) {
    G <- G + 1e-8 * Matrix::Diagonal(nrow(G))
    cholG <- Matrix::Cholesky(G, LDL = FALSE, super = TRUE)
  }
  
  # V^{-1} = D^{-1} - A %*% (G^{-1} S^T D^{-1})   where A = D^{-1} S
  StDinv <- Matrix::crossprod(S, Matrix::solve(cholD, Matrix::Diagonal(nrow(D))))  # S^T D^{-1}
  X      <- Matrix::solve(cholG, StDinv)   # X = G^{-1} S^T D^{-1}
  Vinv   <- Matrix::forceSymmetric(
    Matrix::solve(cholD, Matrix::Diagonal(nrow(D))) - A %*% X,
    uplo = "U"
  )
  
  # log|V| via determinant lemma 
  logdetD <- as.numeric(Matrix::determinant(D, logarithm = TRUE)$modulus)
  logdetQ <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  logdetG <- as.numeric(Matrix::determinant(G, logarithm = TRUE)$modulus)
  logdetV <- logdetD - logdetQ + logdetG
  
  # quadratic form
  quad <- as.numeric(Matrix::crossprod(resid, Vinv %*% resid))
  
  list(Vinv = Vinv, logdetV = logdetV, quad = quad)
}



