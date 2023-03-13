# EBPM
Empirical Bayes Poisson Matrix
===
EBPM is the algorithm of an empirical Bayes method for the Poisson matrix denoising and completion problems

Installation
---
Please directly clone EBPM.R on your computer.  
Operating systems: Linux/Unix/Windows, with R installed.  
Please install the dependency package: CVXR

Input
---
p × q matrix $Y$ (0 or NA in the unobserved matrix entries)
p x q matrix $\Omega$ (0-1 matrix, indication of observed matrix entries)
number of cores $m$

Output
---
p x q matrix $\widehat M$ (The estimated underlying matrix)


Usage 
---
### Matrix denoising
``` r
EB_denoise(Y, n=$m$)
```
### Matrix completion
``` r
EB_denoise(Y, obs=$\Omega$, n=$m$)
```
