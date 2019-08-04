# High dimensional VAR model estimation with strucutural assumption

With extra assumption on transition matrix, we can estimate multivariate high-dimensional VAR model accurately. I decomposed transition matrix into "variable" and "location" matrices respectively. 

## Problem Formulation.
We consider the generic VAR(d) model, for transition matrix M, we decompose this matrix into "variable" effects, which is temporary effect, denoted as A, and "spatial" effect, which is denoted as B. We define the new product, which is a combination of Tensor product and Hadamard product to formulate the objective function. 

## Estimation Procedure.
We adopt ADMM to alternatively update two parameter matrices. 
