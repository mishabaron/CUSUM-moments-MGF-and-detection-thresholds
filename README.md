List and description of R codes for the project and preprint
“Practical Properties of the CUSUM process”, by Michael Baron, and Sergey V. Malov

	Bell polynomials and Splus calculation.R

 Calculations related to rectified partial sums $S_n^+$ 
 Normal distribution case. Change point from F = N(theta0, sigma) to G = Normal(theta1, sigma).
 Let delta = (theta1-theta0)/sigma be the normalized detectable difference.
 Then Y = LLR = Normal( -delta^2/2, delta ) under F 
  and Y = LLR = Normal( +delta^2/2, delta ) under G

	CUSUM MGF.R

 Recursive formula for the CUSUM exponential moment 
There is an option of choosing any value of $\lambda$ and obtaining and graphing MGF as a function of n.
 
	CUSUM mean and variance by recursion.R
 
Recursive formula for the CUSUM mean and variance.
Change-point situation for Normal distributions.

	Doob LHS and RHS comparisons.R

Comparison of left and right hand sides in the Doob’s maximal inequality for CUSUM.
Compares probabilities and thresholds using
 (a) Simulations of CUSUM, we get approximate threshold. 
 (b) Upper and lower bounds using recursive comp and other tools from the paper.
 Normal(0,1) changes to Normal(delta,1), with known parameters 
 [equivalent to any Normal to any Normal]

	Expon moment by recursion.R

Recursive formula for the CUSUM exponential moment 
Compares exact CUSUM exponential moments with our uniform upper bound of (d+1) 
and delta-dependent upper bounds 1 + d*P(Y < 0)

	Lower bound for exp moment.R

Search for a lower bound for $E \exp W_n$

	LowerBound one observ.R

Lower bound for the tail probability of max(CUSUM), obtained by considering one extreme observation, as in our 2023 article in the Journal of Applied Statistics.

	Recursive increments.R
 
Recursive formula for the CUSUM exponential moments and their increments 

	Threshold as a function of delta.R
 
Compare probabilities and thresholds using
 (a) Simulations of CUSUM, we get approximate LHS. 
 (b) Recursive formula for EexpW, we get exact RHS.
 Normal(0,1) changes to Normal(mu,1), known parameters 
 [equivalent to any Normal to any Normal]

