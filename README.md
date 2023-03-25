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
p × q matrix $Y$ (observed matrix, 0 or NA in the unobserved matrix entries)  
p × q matrix $\textit{obs}$ (0-1 matrix, indication of observed matrix entries, 0 in the unobserved entries)  
number of cores $m$ (default value is 1)

Output
---
p x q matrix $\widehat M$ (estimated underlying matrix)


Usage 
---
### Matrix denoising
``` r
source('EBPM.R')
EB_denoise(Y, n=m)
```
### Matrix completion
``` r
source('EBPM.R')
EB_complete(Y, obs, n=m)
```

Demo
---
There are some demo on the simulated data.
Please read the notebook "EBPM.ipynb" in the folder "code of paper".
