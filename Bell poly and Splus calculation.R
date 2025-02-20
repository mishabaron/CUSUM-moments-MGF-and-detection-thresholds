### Calculations related to S_n^+ ###

# Normal distribution case. F = N(theta0, sigma), G = Normal(theta1, sigma).
# Let delta = (theta1-theta0)/sigma.
# Then Y = LLR = Normal( -delta^2/2, delta ) under F 
#  and Y = LLR = Normal( +delta^2/2, delta ) under G 

setwd("C:\\Users\\baron\\Documents\\Research\\Malov\\Spitzer formula paper\\R")
delta = 0.1;   # Choose some delta
N = 1000;      # We'll calculate the first N probabilities
n = 1:N;

Prob = pnorm( -delta*sqrt(n)/2, 0, 1 );

plot(n, Prob)

ESplus = 2*(1-Prob);     # ES_n^+ computed exactly as P_F(Sn<0) + P_G(Sn>0)
plot(n, ESplus)


# Calculate ordinary Bell polynomials

# https://cran.r-project.org/web/packages/kStatistics/kStatistics.pdf

# install.packages("kStatistics")
# library(kStatistics)

# eBellPol(d) = complete exponential Bell polynomials
# eBellPol(d,k) = incomplete exponential Bell polynomials
# e_eBellPol(d,k,c( dim = d )) = evaluated
# oBellPol(d,k) = complete ordinary Bell polynomials
# Not used:
# GCBellPol, e_GCBellPol = generalized complete Bell Polynomial

# Find the forms of incomplete (partial) ordinary Bell polynomials:

# oBellPol(1,1)
# oBellPol(2,1)
# oBellPol(2,2)
# oBellPol(3,1)
# oBellPol(3,2)
# oBellPol(3,3)


# Calculate exponential CUSUM moments with partial ordinary Bell polynomials:

# EU1 = ESplus[1];                                 # EU1 = oBellPol(1,1)
# EU2 = ESplus[2]/2 + 1/2*(ESplus[1])^2;           # EU2 = oBellPol(2,1) + 1/2 * oBellPol(2,2) 
# EU3 = 1/3*ESplus[3] + 1/2*ESplus[1]*ESplus[2] + 1/6*(ESplus[1])^3
# c(EU1,EU2,EU3)



# Calculate exponential CUSUM moments with partial exponential Bell polynomials:

# N = 5; n = 1:N;
# Bell.arg = factorial(n-1)*ESplus[n];   # Arguments for Bell polynomials

# EU = rep(0,N);
# for (d in 1:N){
#   Result = 0;
#   for (k in 1:d){ 
#      Result = Result + e_eBellPol( d, k, c(Bell.arg[1:(d-k+1)],rep(1,k-1)) )
#        # This function requires d arguments, although the last (k-1) are not used
#       }
#   EU[d] = 1/factorial(d) * Result;
#   }

# Compare with the previous calculation: c(EU1,EU2,EU3)

# data.frame(n,EU)
# c(EU1,EU2,EU3)



# Define a function EexpW = EexpW(delta, d) to calculate Eexp(W_d) = EU_d 
# with delta = |theta1-theta0|/sigma = standardized difference between bas and disturbed means
# We'll redefine this function below, using recursive computation of Bell polynomials.
# It will be much faster!

# EexpW = function(d,delta){
#   n = 1:d;
#   ESplus = 2*pnorm( delta*sqrt(n)/2, 0, 1 );     # Expected value of S_n^+
#   Bell.arg = factorial(n-1)*ESplus[n];           # Arguments for Bell polynomials
#   Result = 0;                                    # Calculate the sum adding term by term
#   for (k in 1:d){ 
#      Result = Result + e_eBellPol( d, k, c(Bell.arg[1:(d-k+1)],rep(1,k-1)) )
#        # This function requires d arguments, although the last (k-1) are not used
#      }
#   return( 1/factorial(d) * Result );
#   }

# Check:
# for (k in 1:10){print(EexpW(k,0.1))}
  

# Compare exact CUSUM exponential moments with our uniform upper bound of (d+1)

# Example: d = 1:20, delta = 0.1, 0.5, 1.0, 3.0

# N = 20; d = 1:N; 
# DELTA = c(0.1, 0.5, 1.0, 3.0); 
# COLOR = c("red", "blue", "purple", "orange"); 
# plot(d,d+1,col="black",type="b",xlim=c(1,N),ylim=c(0,N),main="CUSUM exponential moment and our upper bound")

# for (j in 1:length(DELTA)){
#     delta = DELTA[j];
#     EU = rep(0,N); d = 1:N; for (k in 1:N){ EU[k] = EexpW(k,delta); }
#     lines(d,EU,col=COLOR[j],type="b")
#   }
# legend(1, N, legend=c("Exact, delta=0.1", "Exact, delta=0.5", "Exact, delta=1.0", "Exact, delta=3.0", "Upper bound"), col=c(COLOR, "black"), lty=1, cex=0.9)


# For the same scenarios, compare Eexp(S_n^+) with our upper bound = 2
#
# N = 20; d = 1:N; 
# plot(d,rep(2,N),col="black",type="b",lwd=3,xlim=c(1,N),ylim=c(0,2),main="Exponential moments of rectified sums S^+ and our upper bound")
#
# for (j in 1:length(DELTA)){
#     delta = DELTA[j];
#     EexpSplus = 2*pnorm( delta*sqrt(d)/2, 0, 1 );     # Expected value of exp(S_d^+)
#     lines(d,EexpSplus,col=COLOR[j],type="b",lwd=3)
#   }
#
# legend(1, 0.75, legend=c("Exact, delta=0.1", "Exact, delta=0.5", "Exact, delta=1.0", "Exact, delta=3.0", "Upper bound"), col=c(COLOR, "black"), lty=1, cex=0.9)



# Compare exact CUSUM exponential moments with our uniform upper bound of (d+1) 
# and delta-dependent upper bounds 1 + d*P(Y < 0)

# N = 100; d = 1:N; 
# DELTA = c(0.1, 0.5, 1.0, 3.0); 
# COLOR = c("red", "blue", "purple", "orange"); 
# plot(d,d+1,col="black",type="b",xlim=c(1,N),ylim=c(0,N+1), lwd = 3, main="CUSUM exponential moment and our upper bound")
# 
# for (j in 1:length(DELTA)){
#     delta = DELTA[j];
#     EU = rep(0,N); d = 1:N; 
#     for (k in 1:N){ 
#        EU[k] = EexpW(k,delta); 
#        	  filename = paste("update",j,k,".txt")
# 	         write.table(EU,file=filename)
#        }
#     lines(d,EU,col=COLOR[j],type="b", lwd = 3)
#     BetterUpperBound = 1 + d*pnorm( delta/2, 0, 1 );
#     lines(d, BetterUpperBound,col=COLOR[j], lwd = 2, lty=2)
#   }
# 
# legend(1, N+1, legend=c("delta=0.1", "delta=0.5", "delta=1.0", "delta=3.0", "upper bound"), col=c(COLOR, "black"), 
# lty=1, cex=0.9, lwd=2)


######### Direct computation, based on the recursive formula ###########

# Define a function to compute the complete exponential Bell polynomial

Bell_polynomial = function(d, x){
  if (d == 0) { return(1) }
  B = numeric(d+1); B[1] = 1; 
  for (k in 1:d){ B[k+1] = 0;
    for (j in 1:k){
      B[k+1] = B[k+1] + choose(k-1, j-1) * B[j] * x[k-j+1]
  }}  # close both for-loops
  return(B[d+1])
} # end of function Bell_polynomial



# Calculate with the new function Bell_polynomial

EexpW = function(d,delta){
   n = 1:d;
   ESplus = 2*pnorm( delta*sqrt(n)/2, 0, 1 );     # Expected value of S_n^+
   Bell.arg = factorial(n-1)*ESplus[n];           # Arguments for Bell polynomials
   return( 1/factorial(d) * Bell_polynomial(d, Bell.arg) );
   }


# Compare exact CUSUM exponential moments with our uniform upper bound of (d+1) 
# and delta-dependent upper bounds 1 + d*P(Y < 0)

N = 170; d = 1:N; 
DELTA = c(0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0); 
Ndelta = length(DELTA); COLOR = 2:(Ndelta+1);
plot(d,d+1, col="black", type="l", lwd = 2, lty=5,    # 6 = longdash
   xlim=c(1,N),ylim=c(0,N+1), 
   main="CUSUM exponential moments and our upper bounds",
   xlab="delta, the standardized difference")
 
for (j in 1:Ndelta){
    delta = DELTA[j];
    EU = rep(0,N); d = 1:N; 
    for (k in 1:N){ 
        EU[k] = EexpW(k,delta); 
        }
    lines(d,EU,col=COLOR[j], type="l", lwd = 3)
    UpperBound = 1 + d*pnorm( delta/2, 0, 1 );
    lines(d, UpperBound, col=COLOR[j], lwd=2, lty=2)
  }

color.legend = paste("delta =",DELTA[1]);
for (j in 2:Ndelta){ color.legend = c(color.legend, paste("delta =",DELTA[j])) }

legend(1, N+1, legend=c(color.legend,"upper bounds","uniform bound"), 
col=c(COLOR,"red","black"), lty=c(rep(1,Ndelta),2,5), cex=0.9, lwd=2)


# Function psi for Normal distribution

ff = function(w,delta){
 return( pnorm((-w+delta^2/2)/delta)-exp(w)*pnorm((-w-delta^2/2)/delta))
 }

deriv = function(w,delta){
 return( -1/delta*dnorm((-w+delta^2/2)/delta)
  -exp(w)*( pnorm((-w-delta^2/2)/delta))
            -1/delta*dnorm((-w-delta^2/2)/delta)) }

par(mfrow=c(2,2))

w = seq(0,3,0.01); delta = 0.5;
plot(w,ff(w,delta),type="l",col="red",ylim=c(-1,1),lwd=3,xlab="x",ylab="psi(x), psi'(x)", main="delta=0.5")
points(w,deriv(w,delta),type="l",col="blue",lwd=3)
legend(0, 1, legend=c("psi","psi'"), 
col=c("red","blue"), cex=0.9, lwd=2)

w = seq(0,3,0.01); delta = 1;
plot(w,ff(w,delta),type="l",col="red",ylim=c(-1,1),lwd=3,xlab="x",ylab="psi(x), psi'(x)", main="delta=1.0")
points(w,deriv(w,delta),type="l",col="blue",lwd=3)
legend(0, 1, legend=c("psi","psi'"), 
col=c("red","blue"), cex=0.9, lwd=2)

w = seq(0,10,0.01); delta = 2;
plot(w,ff(w,delta),type="l",col="red",ylim=c(-1,1),lwd=3,xlab="x",ylab="psi(x), psi'(x)", main="delta=2.0")
points(w,deriv(w,delta),type="l",col="blue",lwd=3)
legend(8, 1, legend=c("psi","psi'"), 
col=c("red","blue"), cex=0.9, lwd=2)

w = seq(0,10,0.01); delta = 3;
plot(w,ff(w,delta),type="l",col="red",ylim=c(-1,1),lwd=3,xlab="x",ylab="psi(x), psi'(x)", main="delta=3.0")
points(w,deriv(w,delta),type="l",col="blue",lwd=3)
legend(8, 1, legend=c("psi","psi'"), 
col=c("red","blue"), cex=0.9, lwd=2)

par(mfrow=c(1,1))