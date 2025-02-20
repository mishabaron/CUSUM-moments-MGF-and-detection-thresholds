### Recursive formula for the CUSUM mean and variance ###

# Normal distribution case. F = N(theta0, sigma), G = Normal(theta1, sigma).
# Let delta = (theta1-theta0)/sigma.
# Then Y = LLR = Normal( -delta^2/2, delta ) under F 
#  and Y = LLR = Normal( +delta^2/2, delta ) under G 

setwd("C:\\Users\\baron\\Documents\\Research\\Malov\\Spitzer formula paper\\R")


N = 1000; d = 0:N; 
DELTA = c(0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0); 
Ndelta = length(DELTA); COLOR = 2:(Ndelta+1);
plot(d,d+1, col="black", type="l", lwd = 2, lty=5,    # 6 = longdash
   xlim=c(1,N),ylim=c(0,N+1), main="CUSUM exponential moments and our upper bounds",
   xlab="d", ylab="EexpW[d]")


# Recursive computation of EU = EexpW

for (j in 1:Ndelta){
    delta = DELTA[j];
    EU = numeric(N+1);                       # First element is EU[0], then EU[1], etc.
    EU[1] = 1;                               # EU[1] = EexpW[0] = Eexp0 = 1
    EexpSplus = 2*pnorm(delta*sqrt(1:N)/2);  # First element is ESplus[1], ...
    for (n in 1:N){ 
        EU[n+1] = sum(EU[n:1] * EexpSplus[1:n])/n;
        }
    lines(d,EU,col=COLOR[j], type="l", lwd = 3)
    UpperBound = 1 + d*pnorm( delta/2, 0, 1 );
    lines(d, UpperBound, col=COLOR[j], lwd=2, lty=2)
    BetterUpperBound = 1 + d*(pnorm( delta/2, 0, 1 ) - pnorm( -delta/2, 0, 1 ))
    lines(d, BetterUpperBound, col=COLOR[j], lwd=2, lty=3)
  }

color.legend = paste("delta =",DELTA[1]);
for (j in 2:Ndelta){ color.legend = c(color.legend, paste("delta =",DELTA[j])) }

legend(1, N+1, legend=c(color.legend,"old upper bounds","new upper bounds", "uniform bound"), 
col=c(COLOR,"red","red","black"), lty=c(rep(1,Ndelta),2,3,5), cex=0.9, lwd=2)


# Black and white for the paper:
DELTA = c(0.1, 0.5, 1.0, 2.0, 3.0, 5.0); 
Ndelta = length(DELTA); 
plot(2*N, 2*N,          # Outside the plot. This is for limits and titles.
   xlim=c(1,N),ylim=c(0,N+1), main="CUSUM's exponential moments",
   xlab=expression(n), ylab=expression(Ee^(W[n])))

for (j in 1:Ndelta){
    delta = DELTA[j];
    E <- numeric(N+1);                        # First element is E[0], then E[1], etc.
    E[1] <- 1;                                # EU[1] = EexpW[0] = Eexp0 = 1
    x <- 2*pnorm(delta*sqrt(1:N)/2);          # Exponential moments of rectified sums
    for (n in 1:N){                           # Iterative computation
        E[n+1] <- sum(E[n:1] * x[1:n])/n;
        }
    lines(d, E, lty = j, type="l", lwd = 3)
  }

color.legend = substitute(delta == val, list(val = DELTA[1]));
for (j in 2:Ndelta){ color.legend = c(color.legend, substitute(delta == val, list(val = DELTA[j]))) };

legend(1, N+1, legend=color.legend,  lty=1:6, cex=1.5, lwd=2)





# Searching for the asymptote of EU = EexpW

NN = 1000; 
KK = length(NN);
for (kk in 1:KK){
 N = NN[kk]; d = 0:N; 
 DELTA = seq(0.01,4,0.02); 
 Ndelta = length(DELTA); 
 EexpW.N = numeric(Ndelta);

UpperBound = 2*pnorm(DELTA/2)-1;    # Our best upper bound for the increments = P_F(Z>0)+P_G(Z<0)
plot(DELTA,UpperBound,col="red",type="l",lwd=3, 
     xlab="Delta, standardized difference between F and G",
     ylab="Slope of CUSUM's exp. moment and our upper bound",
     main="Slope of EexpW = EexpW(n)/n")

for (j in 1:Ndelta){
    delta = DELTA[j];
    EU = numeric(N+1);                       # First element is EU[0], then EU[1], etc.
    EU[1] = 1;                               # EU[1] = EexpW[0] = Eexp0 = 1
    EexpSplus = 2*pnorm(delta*sqrt(1:N)/2);  # First element is ESplus[1], ...
    for (n in 1:N){ 
        EU[n+1] = sum(EU[n:1] * EexpSplus[1:n])/n;
        }
    EexpW.N[j] = EU[N+1];
  }

 lines(DELTA, EexpW.N/N, col=kk, lwd=3)
 }

legend(0, 1, legend=c("min(N)", "max(N)"), col=c(1,KK), lwd=2)
