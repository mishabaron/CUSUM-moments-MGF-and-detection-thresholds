### Calculations related to S_n^+ ###

# Normal distribution case. F = N(theta0, sigma), G = Normal(theta1, sigma).
# Let delta = (theta1-theta0)/sigma.
# Then Y = LLR = Normal( -delta^2/2, delta ) under F 
#  and Y = LLR = Normal( +delta^2/2, delta ) under G 

delta = 0.1;   # Choose some delta
N = 1000;      # We'll calculate the first N probabilities
n = 1:N;

Prob = pnorm( -delta*sqrt(n)/2, 0, 1 );

plot(n, Prob)

ESplus = 2*(1-Prob);     # ES_n^+ computed exactly as P_F(Sn<0 + P_G(Sn>0)
plot(n, ESplus)

# Calculate ordinary Bell polynomials

install.packages("kStatistics")
library(kStatistics)

# Find the forms of incomplete (partial) ordinary Bell polynomials:

oBellPol(1,1)
oBellPol(2,1)
oBellPol(2,2)
oBellPol(3,1)
oBellPol(3,2)
oBellPol(3,3)

# Calculate exponential CUSUM moments with partial ordinary Bell polynomials:

EU1 = ESplus[1];                                 # EU1 = oBellPol(1,1)

EU2 = ESplus[2]/2 + 1/2*(ESplus[1])^2;           # EU2 = oBellPol(2,1) + 1/2 * oBellPol(2,2) 

EU3 = 1/3*ESplus[3] + 1/2*ESplus[1]*ESplus[2] + 1/6*(ESplus[1])^3

c(EU1,EU2,EU3)

# Calculate exponential CUSUM moments with partial exponential Bell polynomials:

N = 10; n = 1:N;
Bell.arg = factorial(n-1)*ESplus[n];   # Arguments for Bell polynomials

EU = rep(0,N);
for (d in 1:N){
  Result = 0;
  for (k in 1:d){ 
     Result = Result + e_eBellPol( d, k, c(Bell.arg[1:(d-k+1)],rep(1,k-1)) )
       # This function requires d arguments, although the last (k-1) are not used
      }
  EU[d] = 1/factorial(d) * Result;
  }

# Compare with the previous calculation: c(EU1,EU2,EU3)

data.frame(n,EU)
c(EU1,EU2,EU3)


# We can define a function EexpW = EexpW(delta, d) to calculate Eexp(W_d) = EU_d 
# with delta = |theta1-theta0|/sigma = standardized difference between bas and disturbed means

EexpW = function(d,delta){
   n = 1:d;
   ESplus = 2*pnorm( delta*sqrt(n)/2, 0, 1 );     # Expected value of S_n^+
   Bell.arg = factorial(n-1)*ESplus[n];           # Arguments for Bell polynomials
   Result = 0;                                    # Calculate the sum adding term by term
   for (k in 1:d){ 
      Result = Result + e_eBellPol( d, k, c(Bell.arg[1:(d-k+1)],rep(1,k-1)) )
        # This function requires d arguments, although the last (k-1) are not used
      }
   return( 1/factorial(d) * Result );
   }

# Check:
for (k in 1:10){print(EexpW(k,0.1))}
  
# Compare exact CUSUM exponential moments with our upper bound of (d+1)

# Example: d = 1:20, delta = 0.1, 0.5, 1.0, 3.0

N = 20; d = 1:N; 
DELTA = c(0.1, 0.5, 1.0, 3.0); 
COLOR = c("red", "blue", "purple", "orange"); 
plot(d,d+1,col="black",type="b",xlim=c(1,N),ylim=c(0,N),main="CUSUM exponential moment and our upper bound")

for (j in 1:length(DELTA)){
    delta = DELTA[j];
    EU = rep(0,N); d = 1:N; for (k in 1:N){ EU[k] = EexpW(k,delta); }
    lines(d,EU,col=COLOR[j],type="b")
  }

legend(1, N, legend=c("Exact, delta=0.1", "Exact, delta=0.5", "Exact, delta=1.0", "Exact, delta=3.0", "Upper bound"), col=c(COLOR, "black"), lty=1, cex=0.9)


# For the same scenarios, compare Eexp(S_n^+) with our upper bound = 2

plot(d,rep(2,N),col="black",type="b",xlim=c(1,N),ylim=c(0,2),main="Exponential moments of rectified sums S^+ and our upper bound")

for (j in 1:length(DELTA)){
    delta = DELTA[j];
    EexpSplus = 2*pnorm( delta*sqrt(d)/2, 0, 1 );     # Expected value of exp(S_d^+)
    lines(d,EexpSplus,col=COLOR[j],type="b")
  }

legend(1, 0.75, legend=c("Exact, delta=0.1", "Exact, delta=0.5", "Exact, delta=1.0", "Exact, delta=3.0", "Upper bound"), col=c(COLOR, "black"), lty=1, cex=0.9)







