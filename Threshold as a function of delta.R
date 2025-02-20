### BW version. Compare probabilities and thresholds using
### (a) Simulations of CUSUM, we get approximate LHS. 
### (b) Recursive formula for EexpW, we get exact RHS.
### Normal(0,1) changes to Normal(mu,1), known parameters 
### [equivalent to any Normal to any Normal]

N = c(5000,1000,500,100,50); alpha = 0.05;    # Sample size and level of significance
NrunsH = 1000000;           # 1 mln MC runs for simulations 
Delta = c(seq(0.01,0.09,0.01),seq(0.1,9,0.1));
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
    plot(Delta, Threshold, type="l", lwd=2, lty=k, xlab=expression("Magnitude of change " ~ delta),
    ylab=expression("Threshold " ~ " " ~ h), main="Threshold as a function of magnitude", cex=0.75)
             } else { legend.text = c(legend.text,paste("n =",n));
  lines(Delta, Threshold, lwd=2, lty=k) }; # end of if
                 } # end of k-loop
legend("topright", legend.text, lty=1:NN, lwd=2, cex=0.9, seg.len=6, title="Sample size:");
