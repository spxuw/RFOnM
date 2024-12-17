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

(2) R code to obtain gene-wise p-values of dataset from GEO.

(2) Python code to run RFOnM.


# Data type for COP
(1) graph.csv: adjacency matrix of PPI.

(2) feature.csv: Z-transformed gene-wise p-values of each gene. Each column represents the p-values from a omics.


