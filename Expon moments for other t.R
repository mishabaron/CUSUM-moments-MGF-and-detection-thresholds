### Recursive formula for the CUSUM exponential moment ###

###### Try MGF at other values of t ###########
t = 0.999;
DELTA = c(0.1, 0.5, 1.0, 2.0, 3.0, 5.0); 
Ndelta = length(DELTA); 
N=1000; d = 0:N;
plot(2*N, 2*N,          # Outside the plot. This is for limits and titles.
   xlim=c(1,N),ylim=c(0,270),  
   xlab=expression(n), ylab=expression(M[n](lambda)))
title(expression("CUSUM's moment generating function for " * lambda * "=0.999"))

for (j in 1:Ndelta){
    delta = DELTA[j];
    mu = -delta^2/2; sigma = delta;   # E and Var of Y = log-likelihood ratios
    E <- numeric(N+1);                        # First element is E[0], then E[1], etc.
    E[1] <- 1;                                # EU[1] = EexpW[0] = Eexp0 = 1
    x <- exp((1:N)*(t*mu+t^2*sigma^2/2))*pnorm(sqrt(1:N)*(mu/sigma+t*sigma)) + pnorm(-sqrt(1:N)*mu/sigma);
      # x = Exponential moments of rectified sums (at t)
    for (n in 1:N){                           # Iterative computation
        E[n+1] <- sum(E[n:1] * x[1:n])/n;
        }
    lines(d, E, lty = j, type="l", lwd = 3)
  }

color.legend = substitute(delta == val, list(val = DELTA[1]));
for (j in 2:Ndelta){ color.legend = c(color.legend, substitute(delta == val, list(val = DELTA[j]))) };

legend("topleft", legend=color.legend,  lty=1:6, cex=1, lwd=2, bg = "white")
