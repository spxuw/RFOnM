# RFOnM (A statistical physics approach to integrating multi-omics data for disease module detection)

This is a Pytorch and Sklearn implementation of COP, as described in our paper:

Xu-Wen Wang, Min Hyung Ryu, Michael H. Cho, Peter Castaldi, Craig P. Hersh, Edwin K. Silverman, Yang-Yu Liu. [A statistical physics approach to integrating multi-omics data for disease module detection].

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [Data type for DKI](#Data-type-for-DKI)
- [How the use the DKI framework](#How-the-use-the-DKI-framework)

# Overview


# Repo Contents
(1) R codes to generated figures in paper.

(2) Python code to predict the colonized abudnace from the baseline taxanomic profile (Pytorch and Sklearn).


# Data type for COP
(1) X_train.csv: matrix of taxanomic profile of size N*M1, where N is the number of taxa and M1 is the sample size (without header) used to training the prediction models.

(2) X_test.csv: matrix of taxanomic profile of size N*M2, where N is the number of taxa and M2 is the sample size (without header) used to test the prediction model.

(3) y_train.csv: vector of size M1 representing the colonized abudnace of invasion species used to training the prediction models.

(4) y_test.csv: vector of size M2 representing the colonized abudnace of invasion species used to test the prediction models.

