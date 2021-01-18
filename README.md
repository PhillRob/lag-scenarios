# Adapted script to assess population development scenarios

[Rob Hyndman](https://robjhyndman.com/) developed an R function to estimate the existence and length of the invasion time lag in species populations based on species occurrence time series. The [method published in Biological Invasions](https://link.springer.com/article/10.1007/s10530-015-0962-8) detects a period of stagnant population growth prior to an increase. While the script distinguishes species inhibiting the bi-phasic population growth X more scenarios are of interest:
1. No increase: slow/stagnant growth
2. Instant increase: no lag
3. Instant decrease: no lag
4. Decrease with lag
5. Lag phase: Increase after lag ORIGINAL SCRIPT

This repo uses the method to distinguish between the 5 scenarios. 
