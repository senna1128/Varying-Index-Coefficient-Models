# Varying-Index-Coefficient-Models
This repository contains the Matlab and R code for the real application in the paper: High-dimensional Varying Index Coefficient Models via Stein’s Identity, Sen Na, Zhuoran Yang, Zhaoran Wang, Mladen Kolar.

The data folder contains the data downloaded from https://datadryad.org/resource/doi:10.5061/dryad.1139fm7, which has been carefully introduced by the paper https://www.nature.com/articles/s41437-018-0105-y.pdf.

The Figures folder contains all output figures.

The script main.m contains Matlab code for recovering sparse parameter matrix. We require CVX package to be installed in advance. Please see http://cvxr.com/ for instruction.

The script preprodata.R contains R code for preprocessing the data. It outputs GenDat.csv which is the final dataset we analyzed.

Please see the comments in each script for brief introduction of other files.
