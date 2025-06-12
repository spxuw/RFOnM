# RFOnM (A statistical physics approach to integrating multi-omics data for disease module detection)

This is a Pytorch and Sklearn implementation of COP, as described in our paper:

Xu-Wen Wang, Min Hyung Ryu, Michael H. Cho, Peter Castaldi, Craig P. Hersh, Edwin K. Silverman, Yang-Yu Liu. [A statistical physics approach to integrating multi-omics data for disease module detection].

## Contents

- [Installation](#Installation)
- [Data type for RFOnM](#Data-type-for-RFOnM)
- [How the use the RFOnM framework](#How-the-use-the-DKI-framework)

# Installation
<pre>
git clone https://github.com/spxuw/RFOnM.git
cd RFOnM
pip install -e .
</pre>

# Input data for RFOnM
<pre>
(1) graph.csv: adjacency matrix of PPI.

(2) feature.csv: Z-transformed gene-wise p-values of each gene. Each column represents the p-values from an omics.
</pre>

# Usage
<pre>
RFOnM --max_iterations 1000 --n_workers 2 --adj_path demo_graph.csv --feature_path demo_feature.csv --top_size 50  --result_path . 
</pre>

# Parameters
<pre>
max_iterations: the number of iterations for simulated annealing algorithm.
top_size and bottom_size: maximum and minimum module size for refining (500 and 1000 in manuscript).
</pre>

# Output
Output is a one-column csv file indicating the ids of genes in disease module.
