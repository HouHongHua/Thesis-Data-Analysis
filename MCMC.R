## General Parallel Metropolis-Hastings ##
rGPMH <- function(yjc, nu, mu_co, sigma_cco, do.j, p, gpmh_n, gpmh_N)   
{
  pjc <- length(yjc)
  yjc_gen <- matrix(NA, nrow = pjc, ncol = gpmh_n*(gpmh_N + 1))
  yjc_gen[, 1] <- yjc
  yjc_gpmh <- matrix(NA, nrow = pjc, ncol = gpmh_N+1)
  yjc_gpmh[, 1] <- yjc_gen[, 1]
  
  for(j in 2:(gpmh_N+1)) {
    yjc_gpmh[, j] <- c(rtmsl(n = 1, mu_co, sigma_cco, nu, 
                             distr = 'MSL', 
                             lower = rep(-Inf, pjc),
                             upper = yjc)$Y)  
  }  ## for loop
  
  I <- 1 
  K <- numeric(gpmh_N+1)
  A <- numeric(gpmh_N+1) 
  begin = proc.time()[1]
  pb <- txtProgressBar(min = 0, max = gpmh_n, style = 3) 
  
  for (i in 1:gpmh_n) { 
    for(j in 1:(gpmh_N+1)) {
      for (k in 1:(gpmh_N+1)) {
        K[k] <- dtmsl(yjc_gpmh[, j], as.vector(mu_co), sigma_cco, nu, 
                      distr = 'MSL', 
                      a.low = rep(-Inf, pjc),
                      a.upp = yjc)
      } #k
      A[j] <- prod(K[-j]) * dyco.MSL(yjc_gpmh[, j], mu2.1=as.vector(mu_co), 
                                     R22.1=sigma_cco, nu=nu, do.j=as.vector(do.j), p=p)
    }   #j
    
    yjc_gen[, ((i-1)*gpmh_N+i):(i*(gpmh_N+1))] <- yjc_gpmh[, sample(seq(1:(gpmh_N+1)), replace = TRUE, prob = A)]
    
    I <- sample(((i-1)*gpmh_N+i):(i*(gpmh_N+1)), 1)  
    yjc_gpmh[, 1] <- yjc_gen[, I]
    
    for(j in 2:(gpmh_N+1)) {
      yjc_gpmh[, j] <- c(rtmsl(n = 1, mu_co, sigma_cco, nu,
                               distr = 'MSL', 
                               lower = rep(-Inf, pjc),
                               upper = yjc)$Y)
    }
    setTxtProgressBar(pb, i)
  } ## for loop of MH 
  
  burn_in <- matrix(yjc_gen[, -1:-(gpmh_n*(gpmh_N+1)/2)], ncol = gpmh_n*(gpmh_N+1)/2) 
  choose2 <- seq(2, gpmh_n*(gpmh_N+1)/2, 2)
  if(pjc > 1) burn_in <- burn_in[, choose2]
  else burn_in <- matrix(burn_in[, choose2], nrow = 1)
  
  EWW_gen <- array(0, dim = c(pjc, pjc, gpmh_n*(gpmh_N+1)/(2*2)))
  EWW_gen <- array(apply(burn_in, 2, function(x) x %*% t(x)), dim = c(pjc, pjc, gpmh_n*(gpmh_N+1)/(2*2)))
  
  EW <- rowMeans(burn_in)
  EWW <- apply(EWW_gen, c(1, 2), mean)
  
  close(pb)
  end = proc.time()[1]
  GPMH.run.sec = end - begin
  cat("Steal my time", GPMH.run.sec, "seconds.\n")
  return(list(yjc_gen = yjc_gen, burn_in = burn_in, EW = EW, EWW = EWW))
  ## https://github.com/tom-jin/GPMH
}

## Metropolis-Hastings algorithm ##
MH <- function(yjc, nu, mu_co, sigma_cco, do.j, p, mhn)
{
  # begin = proc.time()[1]
  
  pjc <- length(yjc)
  yjc_gen <- matrix(NA, nrow = pjc, ncol = mhn+1) 
  yjc_gen[, 1] <- yjc
  yjc_mh <- c(rtmsl(n = 1, mu_co, sigma_cco, nu, 
                    distr = 'MSL', 
                    lower = rep(-Inf, pjc),
                    upper = yjc)$Y)  ## let it exist
  
  for (i in seq(1, (mhn))) {
    
    yjc_mh <- c(rtmsl(n = 1, mu_co, sigma_cco, nu, 
                      distr = 'MSL', 
                      lower = rep(-Inf, pjc),
                      upper = yjc)$Y)
    
    qyc.old <- dtmsl(yjc_mh, as.vector(mu_co), sigma_cco, nu, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
    pco.old <- dyco.MSL(yjc_mh, mu2.1=as.vector(mu_co), R22.1=sigma_cco, nu=nu, do.j=as.vector(do.j), p=p)
    
    yjc_mh_new <- c(rtmsl(n = 1, mu_co, sigma_cco, nu, 
                          distr = 'MSL', 
                          lower = rep(-Inf, pjc),
                          upper = yjc)$Y)
    
    qyc.new <- dtmsl(yjc_mh_new, as.vector(mu_co), sigma_cco, nu, distr='MSL', a.low=rep(-Inf, pjc), a.upp=yjc)
    pco.new <- dyco.MSL(yjc_mh_new, mu2.1=as.vector(mu_co), R22.1=sigma_cco, nu=nu, do.j=as.vector(do.j), p=p)
    
    a <- min(1, (qyc.old*pco.new) /(qyc.new*pco.old)) 
    if (runif(1) < a) yjc_gen[, i + 1] <- yjc_mh_new else yjc_gen[, i + 1] <- yjc_gen[, i]
  }
  
  yjc_gen <- matrix(yjc_gen[, -1], ncol = mhn, nrow = pjc)
  burn_in <- matrix(yjc_gen[, -1:-(mhn/2)], ncol = mhn/2, nrow = pjc) 
  choose2 <- seq(2, (mhn/2), 2)
  if(pjc > 1) burn_in <- burn_in[, choose2]
  else burn_in <- matrix(burn_in[, choose2], nrow = 1)
  
  EWW_gen <- array(0, dim = c(pjc, pjc, mhn/(2*2)))
  EWW_gen <- array(apply(burn_in, 2, function(x) x %*% t(x)), dim = c(pjc, pjc, mhn/(2*2)))
  
  EW <- rowMeans(burn_in)
  EWW <- apply(EWW_gen, c(1, 2), mean)
  
  # end = proc.time()[1]
  # GPMH.run.sec = end - begin
  # cat("Steal my time", GPMH.run.sec, "seconds.\n")
  
  return(list(yjc_gen = yjc_gen, burn_in = burn_in, EW = EW, EWW = EWW))
}