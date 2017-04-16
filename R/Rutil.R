
dinvchisq <- function (
  ### density function and for the (scaled) inverse-chi-squared distribution.
  x    			##<< vector of quantiles
  , df				##<< degrees of freedom parameter, usually represented as nu
  , scale = 1/df		##<< scale parameter, usually represented as lambda.
  , log = FALSE		##<< Logical. If log=TRUE, then the logarithm of the density is returned.
){
  # adopted from copyLeft pakcage LaplaceDemon, as not available for R 2.10 which runs on the cluster
  ##seealso<< \code{\link{dinvchisq}}
  x <- as.vector(x)
  df <- as.vector(df)
  scale <- as.vector(scale)
  if (any(x <= 0)) 
    stop("x must be positive.")
  if (any(df <= 0)) 
    stop("The df parameter must be positive.")
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  NN <- max(length(x), length(df), length(scale))
  x <- rep(x, len = NN)
  df <- rep(df, len = NN)
  scale <- rep(scale, len = NN)
  nu <- df/2
  dens <- nu * log(nu) - log(gamma(nu)) + nu * log(scale) - 
    (nu + 1) * log(x) - (nu * scale/x)
  if (log == FALSE) 
    dens <- exp(dens)
  dens
}
attr(dinvchisq,"ex") <- function(){
  x <- dinvchisq(1,1,1)
  x <- rinvchisq(10,1)
  
  #Plot Probability Functions
  x <- seq(from=0.1, to=5, by=0.01)
  plot(x, dinvchisq(x,0.5,1), ylim=c(0,1), type="l", main="Probability Function",
       ylab="density", col="red")
  lines(x, dinvchisq(x,1,1), type="l", col="green")
  lines(x, dinvchisq(x,5,1), type="l", col="blue")
  legend(3, 0.9, expression(paste(nu==0.5, ", ", lambda==1),
                            paste(nu==1, ", ", lambda==1), paste(nu==5, ", ", lambda==1)),
         lty=c(1,1,1), col=c("red","green","blue"))	
}


rinvchisq <- function (
  ### random number function and for the (scaled) inverse-chi-squared distribution.
  n				## the number of observations. If length(n) > 1, then the length is taken to be the number required.
  , df
  , scale = 1/df
){
  ##seealso<< \code{\link{rinvchisq}}
  df <- rep(df, len = n)
  scale <- rep(scale, len = n)
  if (any(df <= 0)) 
    stop("The df parameter must be positive.")
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  x <- (df * scale)/rchisq(n, df = df)
  return(x)
}


## see https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf (section 5) for more details

############CONTINUOUS COVARIATES UPDATES##################################
updatevar<-function(tempx,nu0,tau0,c0,mu0){
  newdf<-nu0+length(tempx)
  if(length(tempx)==1){varx<-0}
  if(length(tempx)>1){varx<-var(tempx)}
  numer<-nu0*tau0+(length(tempx)-1)*varx+(c0*length(tempx)/(c0+length(tempx)))*(mean(tempx)-mu0)^2
  newval<-rinvchisq(1,newdf,numer/newdf)
  return(newval)
}

updatemean<-function(tempx,curtau,c0,mu0){
  newvar<-1/(c0/curtau+length(tempx)/curtau)
  newmean<-(mu0*c0/curtau+mean(tempx)*length(tempx)/curtau)*newvar
  newval<-rnorm(1,newmean,sqrt(newvar))
  return(newval)
}
###########################################################################



##############MASS PARAMETER UPDATES#######################################
# mass parameter for y clusters
newalpthet <- function(numclus, prioralp, alpa0, alpb0, nobs){
  eta      <- rbeta(1, prioralp + 1, nobs)
  pieta    <- ( numclus / ( nobs * ( 1 - log(eta) ) ) ) / ( 1 + numclus / ( nobs * ( 1 - log(eta) ) ) )
  whichmix <- rbinom(1, 1, pieta)
  newalp   <- whichmix * rgamma(1, (alpa0 + numclus), (alpb0 - log(eta) ) ) +
              ( 1 - whichmix ) * rgamma(1, (alpa0 + numclus - 1), (alpb0 - log(eta)))
  return(newalp)
}

# mass parameter for x clusters
#update alpha_psi using MH
newalppsi <- function(alpha, cluster, alpa0, alpb0, prop.a, prop.b) {
  sortuniquey <- sort(unique(cluster[ , 1])) ## vector of unique y clusters
  l <- length(sortuniquey)
  k <- nrow(unique(cluster)) ## total number of clusters
  #ss<-numeric(length(sortuniquey))
  #for(j in 1:length(sortuniquey)){
  #  ss[j]<-sum(cluster[,1]==sortuniquey[j])
  #}
  ss       <- as.numeric(table(cluster[,1])) ## number of subjects in each y cluster
  likecur  <- dgamma(alpha, alpa0, alpb0) * ( alpha^(k - l) ) * 
              prod( (alpha + ss) * beta( alpha + 1, ss))
  proposed <- rgamma(1, prop.a, prop.b)
  likeprop <- dgamma(proposed, alpa0, alpb0) * (proposed ^ (k - l)) * 
              prod((proposed + ss) * beta(proposed + 1, ss))
  rat      <- likeprop / likecur
  if (runif(1, 0, 1) < rat) alpha <- proposed
  return(alpha)
}
############################################################################


## count unique number of elements in vector
countUnique <- function(x) return( length( unique(x) ) )

## find first element in vector that is > 2
firstContinuous <- function(x) {
  return( Position( function(p) p > 2,
                    x,
                    right = FALSE,
                    nomatch = length(x) + 1 )  
  )
}
