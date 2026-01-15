# UNDEFINED FOR m > 0 and nu = 0? 
decompos <- function( time, thalf, nu, m, complet)
{
  #_______________________________________________________________________________
  #
  # DEFINITIONS OF CALLING ARGUMENTS: DEFAULT:
  #
  # INPUT
  # time = Time function time
  # thalf = Half-time for integrated cum haz thalf
  # nu = Exponent of time nu
  # m = Exponent of entire function m
  # complet = Output complete dist. function complet formultated:
  # 0=incomplete, 1=complete
  #
  # OUTPUT
  # surviv = Survival surviv
  # cumhaz = Cumulative hazard cumhaz
  # hazard = Hazard function hazard
  # density = Probability density function density
  #
  #__________________________________________________________________________________
  # Compute general factors needed for all models
  #__________________________________________________________________________________
  #
  if (m !=0) {mm1 <- -(1/m)-1;}
  if (nu !=0) {num1 <- -(1/nu)-1;}
  #
  #
  #_________________________________________________________________________________
  # Case 1: M>0, NU>0 model
  #_________________________________________________________________________________
  if (m>0 & nu>0) {
    rho <- nu*thalf*((((2^m)-1)/m)^nu);
    bt <- nu*time/rho;
    btnu <- 1 + m*(bt^(-1/nu));
    capgt <- btnu^(-1/m);
    gt <- (btnu^mm1)*(bt^num1)/rho;
  }
  # ******************************************************************************
  # Limiting Case 1: M=0, NU>0
  else if (m==0 & nu>0) {
    rho <- nu*thalf*(log(2)^nu);
    bt <- nu*time/rho;
    btnu <- bt^(-1/nu);
    capgt <- exp(-btnu);
    gt <- capgt*(bt^num1)/rho;
  }
  # # Limiting Case 2: M>0, NU=0
  # else if (m > 0 & nu == 0) {
  #   FILL IN HERE 
  # }
  #
  #
  #_______________________________________________________________________________
  # Case 2: M<0, NU>0 model
  #_______________________________________________________________________________
  else if (m<0 & nu>0) {
    rho <- nu*thalf/(((1-2^m)^(-nu))-1);
    bt <- 1+nu*time/rho;
    btnu <- 1 - bt^(-1/nu);
    capgt <- btnu^(-1/m);
    gt <- -(btnu^mm1)*(bt^num1)/(m*rho);
  }
  # ******************************************************************************
  # Limiting Case 2: M<0, NU=0
  else if (m<0 & nu==0) {
    rho <- -thalf/(log(1-2^ m));
    bt <- exp(-time/rho);
    btm <- (1-bt);
    capgt <- btm^(-1/m);
    gt <- -(btm^mm1)*bt/(m*rho);
  }
  #
  #
  #_________________________________________________________________________________
  # Case 3: M>0, NU<0 model
  #_________________________________________________________________________________
  else if (m>0 & nu<0) {
    rho <- -nu*thalf*(((2^m)-1)^(nu));
    bt <- -nu*time/rho;
    btnu <- 1 + m*(bt^(-1/nu));
    capgt <- 1-(btnu^(-1/m));
    gt <- (btnu^mm1)*(bt^num1)/rho;
  }
  # ******************************************************************************
  # Limiting Case 3: M=0, NU<0
  #
  # M=0, NU<0 limiting exponential case
  else if (m==0 & nu<0) {
    rho <- -nu*thalf*(log(2)^nu);
    bt <- -nu*time/rho;
    btnu <- bt^(-1/nu);
    capgt <- 1-exp(-btnu);
    gt <- exp(-btnu)*(bt^num1)/rho;
  }
  else {
    stop("Generic function undefined for m/gamma < 0 and nu/eta < 0: m = ", m, ", nu = ", nu)
  }
  
  ht = gt / (1 - capgt)
  # #___________________________________________________________________________________
  # # Complete complete and incomplete functions
  # #___________________________________________________________________________________
  # #
  # # Incomplete distribution function
  # if (complet==0) {
  #   cumhaz <- capgt;
  #   surviv <- exp(-capgt);
  #   hazard <- gt;
  #   density <- gt*capgt;
  # }
  # # ******************************************************************************
  # # Complete distribution function
  # if (complet==1) {
  #   surviv <- 1-capgt;
  #   cumhaz <- -log( surviv);
  #   density <- gt;
  #   hazard <- density/surviv;
  # }
  # Output
  # list(hazard=hazard, density=density, cumhaz=cumhaz, surviv=surviv)
  list(capgt = capgt, gt = gt, ht = ht)
}

get_early_phase <- function(time, thalf, eta = 1, gamma = 0) {
  
  # Note: nu = eta, m = gamma 
  gt <- decompos(time, thalf, nu = eta, m = gamma, complet = 1)$gt
  
  # if(any(is.na(gt) | !is.finite(gt))) {
  #   stop(paste0("gt is not a finite real number: gt =", gt, "nu/eta = ", eta, "m/gamma = ", gamma))
  # }
  # print(paste0("gt = ", gt))
  return(gt)
}

get_late_phase <- function(time, thalf, eta, gamma = 0) {
  
  # Note: nu = eta, m = gamma
  ht <- decompos(time, thalf, nu = eta, m = gamma, complet = 1)$ht

  # if(any(is.na(ht) | !is.finite(ht))) {
  #   stop(paste0("ht is not a finite real number: ht =", ht, "nu/eta = ", eta, "m/gamma = ", gamma))
  # }
  # print(paste0("ht = ", ht))
  return(ht)

  
  # # Note: nu = eta, m = gamma 
  # gt <- decompos(time, thalf, nu = eta, m = gamma, complet = 1)$gt
  # 
  # if(any(is.na(gt) | !is.finite(gt))) {
  #   stop(paste0("gt is not a finite real number: gt =", gt, "nu/eta = ", eta, "m/gamma = ", gamma))
  # }
  # # print(paste0("gt = ", gt))
  # return(gt)
}



