# modules



#' co-clustering based on gaussian latent block
#'
#' This function performs a co-clustering by using a gaussian latent block model.
#' A gaussian distribution is used, and the parameters inference is realized
#' with latent block variational expectation-maximisation algorithm.
#'
#' @param X Data matrix
#' @param g number of row clusters
#' @param m number of col clusters
#' @param n_init number of initialization for VEM and kmeans.
#' @param niter number of iterations for the co clustering algorithm
#'
#' @return list with model outputs
#' ICL: ICL value
#' Z_tild: a posteriori probabilities for row clusters
#' W_tild: a posteriori provabilities for col clusters
#' pi_k: proportions of row clusters
#' rho_l: proportions of col clusters
#' mu: co clusters mean matrix
#' sigma: co clusters variance matrix
#'
#' @importFrom stats kmeans
#' @examples
#' #library(blockcluster)
#' #data(gaussiandata)
#' #model <- lbvem(gaussiandata, g=3, m=4, niter=5)
#' @export
lbvem <- function(X, g, m, n_init = 4, niter = 4){
  n <- nrow(X)
  d <- ncol(X)

  z_init_list <- list()
  w_init_list <- list()
  i_init <- 1
  ICL_init <- c()

  while(i_init < n_init){
    ## Init
    #### initialisation kmeans en lignes
    z_kmeans <- stats::kmeans(X, centers = g, nstart = 1)
    z_init <- z_kmeans$cluster
    z_init_list[[i_init]] <- z_init

    #### initialisation kmeans en colonnes
    w_kmeans <- stats::kmeans(t(X), centers = m, nstart = 1)
    w_init <- w_kmeans$cluster
    w_init_list[[i_init]] <- w_init

    # calcul du modele pour chaque init de kmeans
    model <- LBVEM_core(X, z_init, w_init, g, m, n, d, 1)
    ICL_init <- c(ICL_init, model$ICL)
    i_init <- i_init + 1
  }
  ind_best_init <- which.max(ICL_init)

  model <- LBVEM_core(X, z_init_list[[ind_best_init]], w_init_list[[ind_best_init]], g, m, n, d, niter)
  return( model )
}

LBVEM_core <- function(X, z_init, w_init, g, m, n, d, niter){
  # initialisation de Z_tild, W_tild, pi_k, rho_l, mu_kl, sigma2_kl
  Z_tild <- t(sapply(z_init, function(i) one_hot(i, g)))
  W_tild <- t(sapply(w_init, function(j) one_hot(j, m)))

  # somme des probas que ligne i appartienne au cluster k
  z_tild_ <- sapply(1:g, function(k) sum(Z_tild[, k]))

  # somme des probas que colonne j appartienne au cluster l
  w_tild_ <- sapply(1:m, function(l) sum(W_tild[, l]))

  # proportion des cluster en ligne et en colonne
  pi_k  <- z_tild_ / n
  rho_l  <- w_tild_ / d

  # calcul les mu_kl
  x_wz_tild <- t(Z_tild) %*% X %*% W_tild
  wz_tild <- z_tild_ %*% t(w_tild_)
  mu_kl <- x_wz_tild / wz_tild

  # initialisation des mu_kl precedents (pour la convergence dans les sous-boucles)
  mu_kl_z_old <- matrix(Inf, g, m)
  mu_kl_w_old <- matrix(Inf, g, m)

  # calcul des sigma_kl2
  sigma_kl2 <- ( ( t(Z_tild) %*% (X ** 2) %*% W_tild ) / wz_tild ) - mu_kl** 2

  iter <- 0
  eps_ICL <- 0.1
  ICL <- c(Inf, -Inf)
  eps <- 0.01
  niter_wz <- 5

  while( sum( abs( ICL[length(ICL)] - ICL[length(ICL)-1] ) ) > eps_ICL && iter < niter ){

    if(iter == 1) ICL <- ICL[3]

    x_il_w_tild <- t( apply(X %*% W_tild, 1, function(x) x/ w_tild_) )
    u_il_w_tild <- t( apply((X**2) %*% W_tild, 1, function(x) x/ w_tild_) )
    mu_kl_z_old <- matrix(Inf, g, m)
    iter_wz <- 0

    while( sum( abs(mu_kl - mu_kl_z_old) ) > eps && iter_wz < niter_wz ){
      # step1
      print("E-step: Z")
      Z_tild <- compute_Z_tild(pi_k, w_tild_, mu_kl, sigma_kl2, u_il_w_tild, x_il_w_tild)
      Z_tild <- Z_tild / apply(Z_tild, 1, function(x) sum(x) )

      # step2
      print("M-step: Z")
      z_tild_ <- sapply(1:g, function(k) sum(Z_tild[, k]) )

      pi_k  <- z_tild_ / n
      mu_kl_z_old <- mu_kl
      mu_kl <- apply( t(Z_tild) %*% x_il_w_tild, 2, function(x) x / z_tild_ )
      tmp <- apply( t(Z_tild) %*% u_il_w_tild, 2, function(x) x / z_tild_)
      sigma_kl2 <- tmp - mu_kl ** 2

      # increment iteration
      iter_wz <- iter_wz + 1
      print(paste('iteration: ',iter_wz))
    }

    x_kj_z_tild <- apply(t(Z_tild) %*% X, 2, function(x) x/ z_tild_)
    v_kj_z_tild <- apply(t(Z_tild) %*% (X ** 2) , 2, function(x) x/ z_tild_)
    mu_kl_w_old <- matrix(Inf, g, m)
    iter_wz <- 0

    while( sum( abs(mu_kl - mu_kl_w_old) ) > eps && iter_wz < niter_wz ){
      # step3
      print("E-step: W")
      # W_tild <- compute_W_tild(rho_l, z_tild_, mu_kl, sigma_kl2, v_kj_z_tild, x_kj_z_tild)
      W_tild <- compute_W_tild(rho_l, z_tild_, mu_kl, sigma_kl2, v_kj_z_tild, x_kj_z_tild )
      W_tild <- W_tild / apply(W_tild, 1, sum)

      # step4
      print("M-step: W")
      w_tild_ <- sapply(1:m, function(l) sum(W_tild[, l]))

      rho_l  <- w_tild_/ d
      mu_kl_w_old <- mu_kl
      mu_kl <- t( apply(x_kj_z_tild %*% W_tild, 1, function(x) x/w_tild_) )
      tmp <- t( apply(v_kj_z_tild %*% W_tild, 1, function(x) x/w_tild_) )
      sigma_kl2 <- tmp - mu_kl ** 2

      # increment iteration
      iter_wz <- iter_wz + 1
      print(paste('iteration: ',iter_wz))
    }

    ## critere de convergence
    ICL <- c(ICL, icl(X, Z_tild, W_tild, pi_k, rho_l, mu_kl, sigma_kl2, g, m) )
    iter <- iter + 1
    print(paste("iteration globale: ", iter))
  }
  return( list(ICL=ICL, Z_tild=Z_tild, W_tild=W_tild, pi_k=pi_k, rho_l=rho_l, mu=mu_kl, sigma=sigma_kl2) )
}





one_hot <- function(n_one, len_vector){
  v <- matrix(0, len_vector)
  v[n_one] <- 1
  return(v)
}


compute_Z_tild <- function(pi_, w_tild, mu_, sigma_, u, x_w_tild){
  g <- length(pi_)
  n <- dim(u)[1]
  m <- dim(u)[2]
  Z <- matrix(0, n, g)
  # sigma_ <- sigma_ + 1E-3 * matrix( abs( runif(length(sigma_)) ), g, m)
  sigma_ <- sigma_ + 1E-5 * matrix( 1, g, m)
  for(i in 1:n){
    for(k in 1:g){
      Z[i,k] <- log(pi_[k]) - 0.5 * sum( sapply(1:m, function(l) w_tild[l] * ( log(sigma_[k,l]) + ( u[i,l] - 2*mu_[k,l]*x_w_tild[i,l] + mu_[k,l]**2 )/sigma_[k,l] )))
    }
  }
  if(min(Z) < -100 ){
    print("norm")
    Z <- Z / abs(min(Z)) * 100
  }
  return( exp(Z) )
}

compute_W_tild <- function(rho_, z_tild, mu_, sigma_, v, x_z_tild){

  m <- length(rho_)
  g <- dim(v)[1]
  d <- dim(v)[2]
  W <- matrix(0, d, m)
  # sigma_ <- sigma_ + 1E-3 * matrix( abs(runif(length(sigma_)) ), g, m)
  sigma_ <- sigma_ + 1E-5 * matrix( 1, g, m)
  for(j in 1:d){
    for(l in 1:m){
      W[j,l] <- log(rho_[l]) - 0.5 * sum( sapply(1:g, function(k) z_tild[k] * ( log(sigma_[k,l]) + ( v[k,j] - 2*mu_[k,l]*x_z_tild[k,j] + mu_[k,l]**2 )/sigma_[k,l] )))
    }
  }
  if(min(W) < -100 ){
    print("norm")
    W <- W / abs(min(W)) * 100
  }
  return( exp(W)  )
}

# complete data log likelihood
cdll <- function(X, z, w, pi_, rho_, mus, sigmas){
  ## premier membre de la fonction
  row.entropy <- sum( z %*% log(pi_) )
  ## deixueie membre de fonction
  col.entropy <- sum( w %*% log(rho_) )

  z_ <- apply(z, 2, sum)
  w_ <- apply(w, 2, sum)
  membre3 <- 0
  ## double somme sur k,l
  for(k in 1:ncol(z)){
    for(l in 1:ncol(w)){
      tmp <- 0
      ## double somme sur i,j
      for(i in 1:nrow(X)) {
        for(j in 1:ncol(X)){
          tmp <- tmp + z[i,k]*w[j,l]*((X[i,j] - mus[k,l])**2)
        }
      }
      ## compile des deux doubles sommes
      member3 <- z_[k] * w_[l] * log(sigmas[k,l] + tmp/sigmas[k,l])
    }
  }

  ## compile des 3 membres de la fonctions
  res <- row.entropy + col.entropy - 0.5 * member3
}

icl <- function(X, z, w, pi_, rho_, mus, sigmas, g, m){
  n <- nrow(X)
  d <- ncol(X)
  nu <- 2

  res <- cdll(X, z, w, pi_, rho_, mus, sigmas) -
    0.5*(g-1)*log(n) - 0.5*(m-1)*log(d) - 0.5*g*m*nu*log(n*d)
  return(res)
}
