# Simple population development scenarios

[Rob Hyndman](https://robjhyndman.com/) developed an R function to estimate the existence and length of the time lag in species populations based on species occurrence time series. The [method published in Biological Invasions](https://link.springer.com/article/10.1007/s10530-015-0962-8) detects a period of stagnant population growth prior to an increase. While the original script distinguishes species inhibiting a bi-phasic population growth, 4 more scenarios are of interest:
1. No increase: slow/stagnant growth
2. Instant increase: no lag
3. Instant decrease: no lag
4. Decrease with lag
5. Lag phase: Increase after lag ORIGINAL

This repo fits piecewise linear splines to distinguish between the five scenarios. 

[Scenarios](https://github.com/PhillRob/lag-scenarios/blob/master/lag-scenarios.png)