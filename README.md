# Hypothesis test on weighted networks

This `R` code performs hypothesis testing on multiple aligned graphs using the method described at [1]. 

## How to use?

The function `netCompMM` requires a list with the adjacency matrices and a  dataset with `user`, `time`, `PopId` and `net.name` variables. Here is an example: 
```R
    user time PopId net.name
1     u1    1     1    u1.t1
20   u20    1     1   u20.t1
55    u5    2     1    u5.t2
97   u47    2     1   u47.t2
141  u41    3     1   u41.t3
210  u60    1     2   u60.t1
```

For non heterogenous scenarios, use `time=1` for all items or a named vector with the population information.

## Synthetic data

You can generate synthetic data by sourcing `generateSynthetiData.R`
```R
source("generateSyntheticData.R")
```

## Full example

```R
source("generateSyntheticData.R")

mm <- netCompMM(graphs, side, family = "binom", link = "logit", 
					heterogenity = TRUE, threshold = NULL, matrix.plot = TRUE, 
					R =15, H = 12, Tsamples = 1500, burnIn = 500, by.print = 15, nCores = 8)	

```
## Reference
[1] G. Gomes, V. Rao, and J. Neville, “Multi-level hypothesis testing for populations of heterogeneous networks,” in IEEE 18th International Conference on Data Mining (ICDM) . IEEE, 2018.
