##################################################################
# Beta Rank Function (BRF)
#
# defined by this quantile function:
#  x = p^b/(1-p)^a 
# where p is the cumulative probability, x is the quantile value at p
# when b=0, equivalent to (one-sided) power-law
#
# Example:
#
# from cumulative probability (p) to quantile value (x)
# qbrf(0.9, a=0.5, b=1.5) -> 2.7
# similar to Davies' function  qdavies(0.9, c(1,1.5, 0.5))
#
# from quantile value (x) to cumulative probability (p)
# pbrf(2.7, a=0.5, b=1.5) -> 0.9
# similar to Davies' function pdavies(2.7, c(1,1.5, 0.5))
#
# from quantile value (x) to pdf
# dbrf(2.7, a=0.5, b=1.5) -> 0.055555..
# similar to Davies' function ddavies(2.7, c(1,1.5, 0.5))
#
#  if a,b are not specified, a=1, b=1
##################################################################



####################################
# x = qbrf(p, a=1, b=1)
# require: 0<p<1, a>=0, b>=0
####################################
qbrf <- function(p,a=1, b=1){                                                                                                      
 if(any(p < 0) | any(p > 1) ){
  cat("p=", p, "has to be within (0,1)\n")
  return(NA)
 }
 if(a < 0 | b < 0){
  warning("a or b should not be negative\n")
 }
 logtmp <- b*log(p)-a*log(1-p)
 return( exp(logtmp) )                                                                                                             
}    

####################################
# pdf = dbrf(x, a=1, b=1, h=1E-9, tol=1E-12)
####################################


dbrf <- function(x, a=1, b=1, h=1e-9, tol=1e-12){
 dbrf.point <- function(x0, a=1, b=1, h=1e-9, tol=1e-12){

  root.fun <- function(x){
   r <- uniroot(function(y) (A*((1-y)^b)/y^a)-x,lower=0,upper=1,tol=tol)$root
   return(r)
  }

  x <- seq(x0-2*h,x0+2*h,length.out=5)
  r <- sapply(x, FUN = function(x) root.fun(x))
  f <- -(-r[5]+8*r[4]-8*r[2]+r[1])/(12*h)
  return(f)
 }
 
 y <- sapply(x, FUN=function(x) dbrf.point(x,a=a,b=b,h=h,tol=tol) )
 return(y)
}
 



# ####################################
# # pdf = dbrf(x, a=1, b=1, h=1E-9, tol=1E-12)
# ####################################
# dbrf<-function(x,a=1,b=1,h=1e-9,tol=1e-12){
#  A<-1
#  x<-c(x-2*h,x-h,x,x+h,x+2*h)
#  r<-rep(0,length(x))
#  for(i in 1:length(r)){
#   r[i]<-uniroot(function(y) (A*((1-y)^b)/y^a)-x[i],lower=0,upper=1,tol=tol)$root
#  }
#  FF<-1-r
#  l<-length(FF)
#  F2<-FF[-c(1,2,3,4)]
#  F1<-FF[-c(1,2,3,l)]
#  F0<-FF[-c(1,2,l-1,l)]
#  Fm1<-FF[-c(1,l-2,l-1,l)]
#  Fm2<-FF[-c(l-3,l-2,l-1,l)]
#  f<-(-F2+8*F1-8*Fm1+Fm2)/(12*(x[2]-x[1]))
#  return(f)
# }

####################################
# random number from Beta Rank Function (BRF)
# rbrf(n, a=1, b=1)
# Example: x <- rbrf(1000, a=0.5,b=1.5)
#          plot( density( log(x)), log="y")
#	 or(see below) loglog.hist(x, partition="coarse")
#
# require n > 0, n is interger, a >0  b> 0
####################################

rbrf <- function(n, a=1, b=1){
 if(n < 0){
  cat("n=", n," has to be a positive integer\n")
  return(NA)
 }

 if(a < 0 | b < 0){
  warning("a or b should not be negative\n")
 }

 n <- floor(n)

 x <- qbrf(runif(n),a,b)
 return(x)
}


# rbrf2 <- function(n, a=1, b=1){
#  qbrf <- function(F,b,a){ (F^b)/(1-F)^a }
#  x <- qbrf(runif(n),b,a)
#  return(x)
# }


####################################
# dbrf(x, a=1, b=1, log = FALSE)
# pbrf(x, a=1, b=1, lower.tail = TRUE, log.p = FALSE)
# qbrf(x, a=1, b=1, lower.tail = TRUE, log.p = FALSE)
####################################


##################################################################
# Histogram plot function, where z is log-transformed, and histogram is in log scale.
#
# loglog.hist(z, partition="fine")
#     partition = "coarse","intermediate","fine"
#
# Example:
#  loglog.hist(rbrf(10000,a=0.5,b=0.5),partition="fine")
#  loglog.hist(rbrf(10000,a=0.0,b=0.5),partition="intermediate")
#  loglog.hist(rbrf(10000,a=0.5,b=0.0),partition="coarse")
#
# require z > 0 (otherwise, only position values are used)
##################################################################

loglog.hist <- function(z,partition="fine", main=""){

 if( any(z < 0)){
  warning("x < 0 values are omitted")
  z <- z[z>0]
 }

 logz <- log(z)

 if(partition=="coarse"){
  bins <- 5
 } else if (partition=="intermediate"){
  bins <- 30
 } else if (partition=="fine"){
  bins <- 50
 }

 n <- length(logz)
 k <- round((log(n)/log(2))+1)
 step <- (max(logz)-min(logz))/(bins*k)
 breaks.seq <- seq(min(logz)-step,max(logz)+step,by=step)

 h <- hist(logz,breaks=breaks.seq,plot=FALSE)

 xx <- h$mid
 yy <- h$density *step*100	# y normalized to 100
 xx <- xx[yy>0]
 yy <- yy[yy>0]
 plot(xx,yy,type="b",log="y", xlab="log(x)", ylab="hist(in log)", main=main)
}






