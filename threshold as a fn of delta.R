### Compare probabilities and thresholds using
### (a) Simulations of CUSUM, we get approximate LHS. 
### (b) Recursive formula for EexpW, we get exact RHS.
### Normal(0,1) changes to Normal(mu,1), known parameters 
### [equivalent to any Normal to any Normal]

N = c(5000,1000,500,100); alpha = 0.05;    # Sample size and level of significance
NrunsH = 50;           # MC runs for simulations 
Delta = c(seq(0.01,0.09,0.01),seq(0.1,10,0.1));
NDelta = length(Delta); Wmax = numeric(NrunsH); NN=length(N);

### Function estimating the threshold h as the quantile from the null dist of Wmax
threshold = function(n,mu,alpha,nruns){
  # n = max sample size; mu = disturbed distribution mean. 
  # We are actually calculating the null distribution of the test statistic Wn
  for (run in 1:nruns){
   	 Z = rnorm(n, -mu^2/2, mu);      # Generate log-likelihood ratios
   	 S = cumsum(Z);   	             # Partial sums of z
   	 W = S-cummin(S);	           	   # Cusum process  
    Wmax[run] = max(W);
  } # end for-loop by run
  h = quantile(Wmax,1-alpha);         
  return(h);
  }  # End of function "threshold"

Threshold = numeric(NDelta);

for (k in 1:NN){
  n = N[k];
  for (j in 1:NDelta){ mu = Delta[j];     
     Threshold[j] = threshold(n,mu,alpha,NrunsH);
           }  # end of j-loop
  if (k == 1){  legend.text=paste("n =",n);
    plot(Delta, Threshold, type="l", lwd=2, col=k+2, xlab=expression(delta ~ ", the magnitude of change"),
    ylab=expression("Threshold" ~ h), main="Threshold h as a function of the size of change")
             } else { legend.text = c(legend.text,paste("n =",n));
  lines(MU, Threshold, lwd=2, col=k+2) }; # end of if
                 } # end of k-loop
legend(7.5, 10, legend.text, col=3:(NN+2), lwd=2, cex=0.6, seg.len=5);