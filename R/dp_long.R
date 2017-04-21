#'Longitudinal mixed model with a Dirichlet process prior
#'
#'  @description Fill in later.
#' 
#'  @param y vector of outcomes. may have multiple (longitudinal outcomes) for each subject.
#'  @param trt vector of treatments (may be NULL). only one treatment is allowed per subject.
#'  @param newtrt vector of treatments for hold out subjects whose outcome is predicted.
#'  @param x matrix of covariates. accepts binary and continuous covariates. covariates must be ordered with binary variables first. each row must correspond to one subject so nrow(x) = length(unique(id))
#'  @param newx matrix of covariates for hold out subjects whose outcome is predicted
#'  @param id vector of ids corresponding to y. must have same length as y.
#'  @param timepoints vector of timepoints on which outcome y is drawn.
#'  @param prior list of prior parameters
#'  @param mcmc list of mcmc parameters
#'  @param spline list of spline rules for timepoints. may be NULL if none are desired.
#'  @param verbose logical indicating if the user wants to see printed output in the console. default is FALSE.
#'  @param printevery numeric indicating how often to print updates during MCMC algorithm. used only if verbose = TRUE. default is 100.
#'  
#'  @return Returns a list containing posterior draws and predictions for parameters. Size of output depends on value of mcmc$nsave.
#'  
#'  @export
#'  
#'  @useDynLib EDPlong
#'  
#'  @importFrom splines bs
#'  @importFrom Rcpp evalCpp
#'  @import dplyr
#'  @import MCMCpack
#'  @import mvnfast
dp.long <- function(y, trt, newtrt, x, newx, id, timepoints, prior, mcmc, spline, verbose = FALSE, printevery = 100) {

  N         <- length(id)
  unique.id <- unique(id)
  n         <- length(unique.id)
  n2        <- nrow(newx)   # yields NULL if newx = NULL
	nospline  <- is.null(spline)
	nopred    <- as.logical(mcmc$npred == 0)

	## check validity of arguments
	if( length(y) != length(timepoints)           ) stop( "length of y must be of same length as timepoints" )
	if( length(y) != length(id)                   ) stop( "length of y must be of same length as id" )
	if( length(id) != length(timepoints)          ) stop( "length of id must be of same length as timepoints" )
	if( nrow(x) != n                              ) stop( "number of rows in x must be same as number of unique values in id" )
	if( length(trt) != n & !is.null(trt)          ) stop( "length of trt must equal number of unique values in id" )
	if( any( trt != 0 & trt != 1) & !is.null(trt) ) stop( "treatment must be binary variable with values of 0 or 1" )
	if( !is.list(prior) & !is.null(prior)         ) stop( "prior must be a list or NULL")
	if( !is.list(mcmc) & !is.null(mcmc)           ) stop( "mcmc must be a list")
	if( !is.list(spline) & !is.null(spline)       ) stop( "spline must be a list")
  if( !is.logical(verbose)                      ) stop( "verbose must be a logical variable" )
  if( !is.numeric(printevery) | printevery <= 0 ) stop( "printevery must be an integer > 0")

	printevery <- as.integer(printevery)	

  ## unlist values or set defaults if NULL
	
	# priors
	beta0   <- prior$beta0
	prec0   <- prior$prec0
	beta.a0 <- prior$beta.a0
	beta.b0 <- prior$beta.b0
	a0      <- prior$a0
	b0      <- prior$b0
	a0.u    <- prior$a0.u
	b0.u    <- prior$b0.u
	a0.b    <- prior$a0.b
	b0.b    <- prior$b0.b
	alpa0   <- prior$alpa0
  alpb0   <- prior$alpb0
  mu0     <- prior$mu0
  c0      <- prior$c0
  tau0    <- prior$tau0
  nu0     <- prior$nu0
  prop.a  <- prior$prop.a
  prop.b  <- prior$prop.b

	# splines -- only perform these if splines have been specified
	if (!nospline) {
		nknots <- as.integer(spline$nknots)
  	knots  <- spline$knots
  	degree <- as.integer(spline$degree)
	}
	
	# mcmc
	burnin    		<- as.integer(mcmc$burnin)
	ngibbs    		<- as.integer(mcmc$ngibbs)
	m         		<- as.integer(mcmc$m)
  npred     		<- as.integer(mcmc$npred)
  pred.time 		<- mcmc$pred.time
	pred.time.new <- mcmc$pred.time.new
  nsave     		<- as.integer(mcmc$nsave)
	startY    		<- as.integer(mcmc$startY)
	# startX    		<- as.integer(mcmc$startX)

	pbi       <- ngibbs - burnin  ## post burn-in draws
  pred.rate <- floor( pbi / npred )  ## predict every pred.rate post burn in draws
  save.rate <- floor( pbi / nsave )  ## save every save.rate post burn in draws

	if (verbose) {
		cat("Total MCMC iterations: ", ngibbs, "\n")
	  cat("  Burn-in draws: ", burnin, "\n")
	  cat("  Post burn-in draws: ", pbi, "\n\n")
		cat("Saving parameter values every ", save.rate, " draws post burn-in\n") 
		cat("Predicting outcome every ", pred.rate, " draws post burn-in\n\n")
		cat("-------------------------------------------------------------------\n")
		cat(n, " subjects are contributing ", N, "observations\n")
		if (is.null(n2)) cat("There are no additional subjects to predict (from newx argument)\n")
		else cat("There are an additional ", n2, "subjects to predict (from newx argument)\n")
	}

	## check validity of mcmc list values
  if( !is.numeric(burnin) | burnin < 0 | burnin > ngibbs       ) stop( "burnin must be an integer >=0 and less than ngibbs" )
  if( !is.numeric(ngibbs) | ngibbs <= 0                        ) stop( "ngibbs must be an integer > 0" )
  if( !is.numeric(m) | m <= 0                                  ) stop( "m must be an integer > 0" )
  if( !is.numeric(startY) | startY <= 0                        ) stop( "startY must be an integer > 0" )
  # if( !is.numeric(startX) | startX <= 0                        ) stop( "startX must be an integer > 0" )
  if( !is.numeric(npred) | npred < 0 | npred > ngibbs - burnin ) stop( "npred must be an integer >= 0 and be less than ngibbs - burnin" )
  if( !is.numeric(nsave) | nsave < 0 | nsave > ngibbs - burnin ) stop( "nsave must be an integer >= 0 and be less than ngibbs - burnin" )
	if( length(pred.time) != 1 & length(pred.time) != n          ) stop( "pred.time must be of length 1 or length(unique(id))" )
	if( !is.null(newx) ) {
		if( length(pred.time.new) != 1 & length(pred.time.new) != n2 ) stop( "pred.time2 must be of length 1 or nrow(newx)" )
  }

	## check validity of prior list values
  #if( length(beta0) != ncol(x) + 1                ) stop( "beta0 must be of length ncol(x) + 1" )
  if( ncol(prec0) != length(beta0)                ) stop( "ncol(prec0) must be length of beta0" )
  if( nrow(prec0) != length(beta0)                ) stop( "nrow(prec0) must be length of beta0" )
  if( beta.a0 <= 0 | beta.b0 <= 0                 ) stop( "beta.a0 and beta.b0 must be positive numbers" ) 
  if( !is.numeric(beta.a0) | !is.numeric(beta.b0) ) stop( "beta.a0 and beta.b0 must be numeric" ) 
  if( a0 <= 0 | b0 <= 0                           ) stop( "a0 and b0 must be positive numbers" ) 
  if( !is.numeric(a0) | !is.numeric(b0)           ) stop( "a0 and b0 must be numeric" ) 
  if( a0.u <= 0 | b0.u <= 0                       ) stop( "a0.u and b0.u must be positive numbers" ) 
  if( !is.numeric(a0.u) | !is.numeric(b0.u)       ) stop( "a0.u and b0.u must be numeric" ) 
  if( a0.b <= 0 | b0.b <= 0                       ) stop( "a0.b and b0.b must be positive numbers" ) 
  if( !is.numeric(a0.b) | !is.numeric(b0.b)       ) stop( "a0.b and b0.b must be numeric" ) 
  if( alpa0 <= 0 | alpb0 <= 0                     ) stop( "alpa0 and alpb0 must be positive numbers" ) 
  if( !is.numeric(alpa0) | !is.numeric(alpb0)     ) stop( "alpa0 and alpb0 must be numeric" ) 
	if( !is.numeric(mu0)                            ) stop( "mu0 must be numeric" )
  if( !is.numeric(c0) | c0 <=0                    ) stop( "c0 must be numeric and greater than 0" )
  if( !is.numeric(tau0) | tau0 <=0                ) stop( "tau0 must be numeric and greater than 0" )
  if( !is.numeric(nu0) | nu0 <=0                  ) stop( "nu0 must be numeric and greater than 0" )
  if( prop.a <= 0 | prop.b <= 0                   ) stop( "prop.a and prop.b must be positive numbers" ) 
  if( !is.numeric(prop.a) | !is.numeric(prop.b)   ) stop( "prop.a and prop.b must be numeric" ) 

	## check validity of spline list values
	if (!nospline) {
  	if( !is.numeric(nknots) | nknots < 0                   ) stop( "nknots must be an integer >=0" )
  	if( !is.numeric(knots) | length(knots) != nknots       ) stop( "knots must be a numeric vector of length nknots" )
  	if( nknots > 0 & ( !is.numeric(degree) | degree <= 0 ) ) stop( "degree must be integer > 0" )
	}

	if( !is.null(trt) ) {
		x    <- cbind(trt, x)
		ptrt <- 1
    if ( !is.null(n2) ) {
      newx <- cbind(newtrt, newx)
    }
	} else {
		ptrt <- 0
	}

  
	## calculate variable types in x matrix	
	num.unique <- apply(x, 2, countUnique)
	first.cont <- firstContinuous(num.unique)
	p1         <- first.cont - 1 - ptrt       ## number of discrete non-treatment variables
	p2         <- ncol(x) - first.cont + 1    ## number of continuous variables

	## calucate number of timepoints per person
	max.per.person <- data.frame(id = id) %>% 
  				group_by(id) %>% 
  				summarise(n = n())
	mpp <- max.per.person$n


	## make matrix with X for regression
	long.rows <- rep( 1:nrow(x), times = mpp )
	mat.all   <- x[ long.rows , ]
	mat.all   <- cbind( 1, mat.all )  ## add intercept
  mat.all   <- cbind( mat.all, timepoints )  ## add main effect for time
	## make matrix for spline
	if (!nospline) {
		Z  <- splines::bs(timepoints, knots = knots, degree = degree)
  	nZ <- ncol(Z)
	} else {
		Z <- NULL
		nZ <- 0
	}

	## make matrix for predictions
	if (!nopred) {
		if( length(pred.time) == 1) pt <- rep(pred.time, n)
		else                        pt <- pred.time
	} else {
		pt <- NULL
	}

	if (!is.null(n2)) {
		if( length(pred.time.new) == 1) pt.new <- rep(pred.time.new, n2)
		else                        pt.new <- pred.time.new
	} else {
		pt.new <- NULL
	}

	if (!nospline & !nopred) pZ <- predict(Z, pt)
	else pZ <- NULL

	if (!nospline & !is.null(n2)) pZ.new <- predict(Z, pt.new)
	else pZ.new <- NULL

	
	#### initialize cluster
	s <- numeric(n)

	## min and max y values for initial clustering
	yvals <- data.frame(y = y, id = id) %>% group_by(id) %>% summarise( max.y = max(y) , min.y = min(y) )
	max.y <- as.vector(yvals[ ,2])
  min.y <- as.vector(yvals[ ,3])

	## y cluster
	distyx    <- dist(cbind(x, max.y, min.y))
	hi        <- hclust(distyx, method = "ward.D")
	s         <- cutree(hi, k = startY )


	#### initialize containers for parameters
	#columns: # cont, # binary, # trt, intercept, # time coefficients
	beta.reg <- matrix(NA , nrow = length( unique(s) ) , ncol = length(beta0) )
	sig.reg <- numeric( length( unique(s) ) )

	if (!nospline) {
		b.reg  <- matrix(0 , nrow = length( unique(s) ) , ncol = nZ )  # regression parameters on knots
		sig2.b <- numeric( length( unique(s) ) )
	} else {
		b.reg  <- NULL
		sig2.b <- NULL
	}
  
	x.pi.pars  <- NULL
	x.mu.pars  <- NULL
	x.sig.pars <- NULL
	
	if (p1 + ptrt > 0) {
	  x.pi.pars  <- matrix(0 , nrow = length(unique(s)), ncol = p1 + ptrt)
	}
	
	if (p2 > 0) {
	  x.mu.pars  <- matrix(0 , nrow = length(unique(s)), ncol = p2)
	  x.sig.pars <- matrix(0 , nrow = length(unique(s)), ncol = p2)
	}
	
  # initialize random intercept
	u.int <- rep(0, n)

	#### initialize parameter values
	## regression parameters
	get.reg <- samp_reg(y,
											as.matrix(mat.all),
											id,
											Z,
											rep(u.int, mpp),
											beta.reg, sig.reg, b.reg,
											s ,
											beta0, prec0, beta.a0, beta.b0,
											sig2.b, a0.b, b0.b)

	beta.reg <- get.reg$betaY
	sig.reg  <- as.vector(get.reg$sig2)
	if (!nospline) {
		b.reg    <- get.reg$b
		sig2.b   <- as.vector(get.reg$sig2_b)
	}

	## todo:all thhis not needed. pare down
	uniques             <- unique(s)
	sortedunique        <- uniques[ order(uniques) ]
	
	covs <- update_covs(x, s, sortedunique, x.pi.pars,
	                    x.mu.pars, x.sig.pars,
	                    p1, p2, ptrt,
	                    a0, b0, 
	                    mu0, nu0, tau0, c0)
	
	if (p1 + ptrt > 0 & p2 > 0) {
	  x.pi.pars <- covs$xpipars
	  x.mu.pars <- covs$xmupars
	  x.sig.pars <- covs$xsigpars
	} else if (p2 == 0) {
	  x.pi.pars <- covs$xpipars
	} else if (p1 + ptrt == 0) {
	  x.mu.pars <- covs$xmupars
	  x.sig.pars <- covs$xsigpars
	}
	


  # initialize
  sig2.u <- 1.0
  
	## random intercept
	get.u <- samp_u(y, as.matrix(mat.all), Z,
       						beta.reg, b.reg,
       						id, s,
       						sig2.u, sig.reg,
       						n)
	u.int <- get.u$u
	sig2.u <-  1/rgamma(1,a0.u+(n/2),b0.u + (0.5)*crossprod(u.int,u.int))


	#### calculations for predictions
	if ( !is.null(n2) ) {
		h0x <- matrix( NA , nrow = n2 , ncol = ptrt + p1 + p2 )
		for(i in 1:n2) {
      if (ptrt + p1 > 0) {
        for(j in 1:( ptrt + p1 ) ) {
          h0x[i,j] <- beta(a0 + newx[ i , j ] , b0 - newx[ i , j ] + 1 )
        }
      }
 
      if (p2 > 0) {
        for(j in 1:p2) {
          denom <- ( sqrt( 2 * pi ) / sqrt( c0 ) ) * gamma( nu0 / 2 )*( 2 / ( nu0 * tau0 ) ) ^ ( nu0 / 2 )
          cn <- c0 + 1
          nun <- nu0 + 1
          taun <- ( 1 / nun ) * ( nu0 * tau0 + ( c0 / ( c0 + 1 ) ) * ( mu0 - newx[ i , p1 + ptrt + j ] ) ^ 2 )
          num <- ( sqrt( 2 * pi ) / sqrt( cn ) ) * gamma( nun / 2 ) * ( 2 / ( nun * taun ) ) ^ ( nun / 2 )
          h0x[i, ptrt + p1 + j] <- num / ( denom * sqrt( 2 * pi ) )
        }
      }
		}
	} else h0x <- NULL

  ## initialization
  alpha   <- 1
  
	#### containers for saving stuff
	s.save           <- vector("list", nsave)
  beta.reg.save    <- vector("list", nsave)
  b.reg.save       <- vector("list", nsave)
  sig.reg.save     <- vector("list", nsave)
  u.int.save       <- vector("list", nsave)
  x.pi.save        <- vector("list", nsave)
  x.mu.save        <- vector("list", nsave)
  x.sig.save       <- vector("list", nsave)
  sig.u.save       <- numeric(ngibbs)
	alpha.rep        <- numeric(ngibbs)
	total.clusters   <- numeric(ngibbs)
	total.Y.clusters <- numeric(ngibbs)
	
	pred.w.data      <- matrix(NA, nrow = npred, ncol = n)
	
	if(!is.null(n2)) {
    pred.wo.data         <- matrix(NA, nrow = npred, ncol = n2)
    pred.wo.data.cluster <- matrix(NA, nrow = npred, ncol = n2)
  } else {
    pred.wo.data         <- NULL
    pred.wo.data.cluster <- NULL
  }
	

	#### iterate through MCMC draws
  if (verbose) cat("Beginning MCMC algorithm\n")
	beg.time     <- Sys.time()
	int.beg.time <- Sys.time()
	for(i in 1:ngibbs) {
		
  	uniques <- unique(s)
  	unique0 <- uniques[order(uniques)]

		## choose cluster

  	cluster_res <- cluster(y,  #y values
													 as.matrix(x),     #matrix with n rows
													 as.matrix(mat.all),      #matrix with row for each time point
													 mpp,  # time points per person
													 id,    #ids 
													 Z,
													 u.int,
													 beta.reg, sig.reg, b.reg,  #Y params
													 x.pi.pars, x.mu.pars, x.sig.pars, #X params
													 ptrt, p1, p2,
													 alpha,
													 s,
													 unique0, #cluster
													 beta0, prec0, beta.a0, beta.b0, #priors on Y params
													 sig2.b, a0.b, b0.b,
													 a0, b0, #priors on discrete X
													 mu0, nu0, tau0, c0, #priors on cont X
													 alpa0, alpb0, #priors on alpha
													 m # no. of aux params
  												 )
  	
  	s          <- cluster_res$Sy 
  	beta.reg   <- cluster_res$betaY
  	sig.reg    <- as.vector(cluster_res$sig2)
  	x.pi.pars  <- cluster_res$xpipars
  	x.mu.pars  <- cluster_res$xmupars
  	x.sig.pars <- cluster_res$xsigpars
		if (!nospline) {
  		b.reg      <- cluster_res$b
  		sig2.b     <- as.vector(cluster_res$sig2_b)
		}

  	uniques             <- unique(s)
  	sortedunique        <- uniques[order(uniques)]
  	k                   <- length(sortedunique) #number of y clusters
  	total.clusters[i]   <- k


		####### update regression parameters #######
  	get.reg <- samp_reg(y,
												as.matrix(mat.all),
												id,
												Z,
												rep(u.int, mpp),
    	                  beta.reg, sig.reg, b.reg,
      	                s,
        	              beta0, prec0, beta.a0, beta.b0,
          	            sig2.b, a0.b, b0.b)
 
  
  	beta.reg <- get.reg$betaY
  	sig.reg  <- as.vector(get.reg$sig2)
  	if (!nospline) {
			b.reg    <- get.reg$b
  		sig2.b   <- as.vector(get.reg$sig2_b)
  	}

		####### update random intercept #######
  	get.u <- samp_u( y , mat.all , Z ,
    	              beta.reg , b.reg ,
      	            id , s ,
        	          sig2.u , sig.reg ,
          	        n )
  	u.int <- get.u$u




		####### update covariate parameters #######
  	covs <- update_covs(x, s, sortedunique, x.pi.pars,
  	                    x.mu.pars, x.sig.pars,
  	                    p1, p2, ptrt,
  	                    a0, b0, 
  	                    mu0, nu0, tau0, c0)
  	
  	if (p1 + ptrt > 0 & p2 > 0) {
  	  x.pi.pars <- covs$xpipars
  	  x.mu.pars <- covs$xmupars
  	  x.sig.pars <- covs$xsigpars
  	} else if (p2 == 0) {
  	  x.pi.pars <- covs$xpipars
  	} else if (p1 + ptrt == 0) {
  	  x.mu.pars <- covs$xmupars
  	  x.sig.pars <- covs$xsigpars
  	}
  
	

		####### update variance for random intercept ########
  	sig2.u <- sig.u.save[i] <- 1/rgamma(1,a0.u+(n/2),b0.u + (0.5)*crossprod(u.int,u.int))

		####### update alpha parameters #######
  	#  alpha
  	alpha        <- newalpha(k , alpha , alpa0 , alpb0 , n )
  	alpha.rep[i] <- alpha


		##########################################
		##########################################
		######### predictions ####################
		##########################################
		##########################################

		## todo: update function to handle additional parameters
  	if(i > burnin & (i - burnin) %% pred.rate == 0) {
    	preds <- pred(Xonly       = as.matrix(x),
      	            Xonly2b     = newx,
        	          h0xb        = h0x,
          	        betaY       = beta.reg,
            	      sig2        = sig.reg,
              	    bregb       = b.reg,
                	  xpipars     = x.pi.pars,
                 		xmupars     = x.mu.pars,
                  	xsigpars    = x.sig.pars,
                  	ptrt        = ptrt,
                  	p1          = p1,
                  	p2          = p2,
                  	alpha       = alpha,
                  	Sy          = s,
                  	uniqueS     = sortedunique,
                  	beta0       = beta0,
                  	prec0       = prec0,
                  	beta_a0     = beta.a0,
                  	beta_b0     = beta.b0,
                  	a0_b        = a0.b,
                  	b0_b        = b0.b,
                  	timepoint   = pt,
										timepoint2b = pt.new,
                  	tZb         = pZ,
										tZ2b        = pZ.new )

	    pred.w.data[(i - burnin) / pred.rate, ]          <- preds$pred1
  	  if (!is.null(n2)) {
        pred.wo.data[(i - burnin) / pred.rate, ]         <- preds$pred2
  	    pred.wo.data.cluster[(i - burnin) / pred.rate, ] <- preds$pred2_clust
      }
  	}



		##########################################
		##########################################
		########### save values ##################
		##########################################
		##########################################

  	if(i > burnin & (i - burnin) %% save.rate == 0) {
			s.save[[(i - burnin) / save.rate]]           <- s
  		beta.reg.save[[(i - burnin) / save.rate]]    <- beta.reg
  		b.reg.save[[(i - burnin) / save.rate]]       <- b.reg
  		sig.reg.save[[(i - burnin) / save.rate]]     <- sig.reg
  		u.int.save[[(i - burnin) / save.rate]]       <- u.int
  		x.pi.save[[(i - burnin) / save.rate]]        <- x.pi.pars
  		x.mu.save[[(i - burnin) / save.rate]]        <- x.mu.pars
  		x.sig.save[[(i - burnin) / save.rate]]       <- x.sig.pars
		}
		##########################################
		##########################################
		########## print status  #################
		##########################################
		##########################################

  	if ( (verbose) & (i %% printevery == 0) ) {
  	  cat("Iteration: ",i," of ",ngibbs,"\n")
  	  cat("Time elapsed: ", Sys.time() - int.beg.time, "\n\n")
  	  int.beg.time <- Sys.time()
    }

	}  ## end of gibbs loop

  if (verbose) {
    cat("end of MCMC algorithm\n")
    cat("Total time elapsed: ", Sys.time() - beg.time, "\n")
  }
  return( list( s                = s.save, 
								beta.reg         = beta.reg.save, 
								b.reg            = b.reg.save,
						  	sig.reg          = sig.reg.save, 
								u.int            = u.int.save,
								x.pi.pars        = x.pi.save, 
								x.mu.pars        = x.mu.save, 
								x.sig.pars       = x.sig.save,
								sig.u            = sig.u.save, 
								alpha      = alpha.rep, 
								pred.w.data      = pred.w.data,
								pred.wo.data     = pred.wo.data,
                pred.wo.datc     = pred.wo.data.cluster,
                total.clusters   = total.clusters) )

}
