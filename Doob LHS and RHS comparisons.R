### Similar to Doob LHS and RHS comparisons.R
### Creates BLACK and WHITE plots in separate ps files.

### Compare probabilities and thresholds using
### (a) Simulations of CUSUM, we get approximate threshold. 
### (b) Upper and lower bounds using recursive comp and other tools from our 2025 paper
### Normal(0,1) changes to Normal(delta,1), known parameters 
### [equivalent to any Normal to any Normal]

N = 1000; alpha = 0.05;   # Sample size and level of significance
NrunsH = 5000;           # MC runs for simulations 
Delta = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 4)
NDelta = length(Delta);


### Function estimating the threshold h as the quantile from the null dist of Wmax
thresholds = function(n,delta,alpha,nruns){
  # n = max sample size; delta = disturbed distribution mean (in general, the standardized difference). 
  # We are actually calculating the null distribution of the test statistic Lambda = max(Wn)
  Wmaxmat = matrix(0,nruns,n);             # running max(W) for 1..n 
  for (run in 1:nruns){
   	 Z = rnorm(n, -delta^2/2, delta);    # Generate log-likelihood ratios
   	 S = cumsum(Z);   	             # Partial sums of z
   	 W = S-cummin(S);	           	       # Cusum process 
     Wmaxmat[run,] = cummax(W);            # running max(W) for this MC run 
  } # end for-loop by run
  h = numeric(n);                          # thresholds h = h(n)
  for (k in 1:n){ h[k] = quantile(Wmaxmat[,k],1-alpha); }        
  return(h);
  }  # End of function "thresholds"

ThresholdSIM = numeric(N);  # Threshold estimated by simulations as an alpha-upper quantile of max(W)(n)
ThresholdUB  = numeric(N);  # Upper bound for the threshold calculated as log(EexpW/alpha), as in Sec. 6.1
ThresholdUBT = numeric(N);  # "Tighter" threshold as in Sec. 6.3
ThresholdLB  = numeric(N);  # Lower bound from 2025 paper, with an approx optimal k      

EU = numeric(N); 

par(mfrow=c(2,4), oma = c(4, 4, 4, 2));       

for (j in 1:NDelta){ delta = Delta[j];        # Recursive computation of EU = EexpW
    EexpSplus = 2*pnorm(delta*sqrt(1:N)/2);   # P_F(Sn<0)+P_G(Sn>0), in paper, formula (3.5)
    EU[1] = EexpSplus[1];
    for (n in 2:N){ 
        EU[n] = (EexpSplus[n] + sum(EU[1:(n-1)] * EexpSplus[(n-1):1]))/n;
          # Now we derive the lower bound, using the method of splitting n into [n/k] segments
        IX = function(x){ return( log(x) + (x+delta^2/2)^2/2/delta^2 - log(n) )}   # Function for finding k, take log of both sides to avoid huge numbers
        if (IX(1) >= 0){k=1} else {k = (uniroot(IX, interval=c(1,n))$root)};       # Actually, need to round this root to an integer
        ThresholdLB[n] = -k*delta^2/2 + delta*sqrt(k)*qnorm((1-alpha)^(1/floor(n/k)));  # Lower bound from 2025 paper, with an approx optimal k      
        }
    ThresholdSIM = thresholds(N,delta,alpha,NrunsH);    # Threshold estimated by simulations as an alpha-upper quantile of max(W)(n)
    ThresholdUB = log(EU/alpha);                        # Upper bound for the threshold calculated as log(EexpW/alpha), as in Sec. 6.1
    D = pnorm(delta/2)-pnorm(-delta/2);                 # Our new discrepancy measure between F and G
    ThresholdUBT = log((1+D*(1:N))/alpha);              # "Tighter" threshold as in Sec. 6.3
    ThresholdUBU = log(((1:N)+1)/alpha);                # "Universal" threshold as in Sec. 6.2 
      # Using our best upper bound EU(n) <= 1 + Mn, 
      # where $M = \psi(0) = \P_F(S\le 0)-P_G(S\le 0)$

    ThresholdLB1 = -delta^2/2 + delta*qnorm((1-alpha)^(1/(1:N)));    # Lower bound from 2023 paper, sec. 5.2 , equation (22)

    YL = 10; # ThresholdUBT[N];

    ThresholdUBU <- ifelse(ThresholdUBU > YL, NA, ThresholdUBU)
    ThresholdLB <- ifelse(ThresholdLB < 0, NA, ThresholdLB)
    ThresholdLB1 <- ifelse(ThresholdLB1 < 0, NA, ThresholdLB1)

    plot(ThresholdUB,type="l", lwd=2, lty=2, main = substitute(delta == val, list(val = delta)),
         ylim=c(0,YL), xlab="n", ylab="h");  
    lines(ThresholdSIM, lwd=3, lty=1);   
    lines(ThresholdUB,  lwd=2, lty=2); 
    lines(ThresholdUBU, lwd=2, lty=3); 
    lines(ThresholdUBT, lwd=2, lty=4); 
#    lines(ThresholdLB, lwd=2, lty=5);
    lines(ThresholdLB1, lwd=2, lty=6);

  if (delta < 3){    lines(ThresholdLB, lwd=2, lty=5);  # Lower bound, based on partitions
# lines(ThresholdLB1,  lwd=2, lty=6 );   # Lower bound, based on 1 obs.
   }
} # End of the k-loop, over scenarios

# Add an overall text
mtext("Threshold h", side = 2, outer = TRUE, cex = 0.9) 
mtext("    Sample size n", side = 1, outer = TRUE, cex = 0.9) 
mtext("Detection threshold. Exact value and bounds", side = 3, outer = TRUE, cex = 1) 

# Add an overall legend
# par(xpd = NA)  # Allow plotting outside the region

par(mfrow=c(1,1))

legend(x = -36, y = -5.7, legend = c("Exact ", "UB-1", "UB-2", "UB-3", "LB-1", "LB-2"),
       lty = (1:6), lwd = 2, cex = 0.6, seg.len = 6, yjust = 0, horiz = TRUE)

