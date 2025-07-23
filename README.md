# RFOnM (A statistical physics approach to integrating multi-omics data for disease module detection)

This is a Python implementation of RFOnM, as described in our paper:

Xu-Wen Wang, Min Hyung Ryu, Michael H. Cho, Peter Castaldi, Craig P. Hersh, Edwin K. Silverman, Yang-Yu Liu. [A statistical physics approach to integrating multi-omics data for disease module detection].

## Contents

- [Installation](#Installation)
- [Data type for RFOnM](#Data-type-for-RFOnM)
- [How the use the RFOnM framework](#How-the-use-the-DKI-framework)

# Installation
<pre>
  conda create -n rfonm python=3.12
  conda activate test_rfonm
  conda install -c xuwenwang rfonm
</pre>

# Input data for RFOnM
<pre>
(1) graph.csv: adjacency matrix of PPI.

(2) feature.csv: Z-transformed gene-wise p-values of each gene. Each column represents the p-values from an omics.
</pre>

Full dataset can be downloaded from [Dropbox](https://www.dropbox.com/scl/fo/oha0igt23h15bw5ddhu06/AOAODCkDfPjQrX5epIfKVCM?rlkey=6op1kp30qedd54eg3jrvucax4&dl=0).

# Usage
<pre>
rfonm --max_iterations 1000 --n_workers 2 --adj_path ../data/graph_Alzheimer.csv --feature_path ../data/features_Alzheimer.csv --top_size 1000 --bottom_size 500 --result_path . 
</pre>

Results of RFOnM and other methods can be downloaded from [Dropbox](https://www.dropbox.com/scl/fo/6pedlejdl9m2q6pqk8gjj/APeqlZ96admY6uisHKcJjS8?rlkey=dxix574m844widitj8ixkgi8l&dl=0)
# Parameters
<pre>
max_iterations: the number of iterations for simulated annealing algorithm.
top_size and bottom_size: maximum and minimum module size for refining (500 and 1000 in manuscript).
</pre>

# Output
Outputa are a one-column csv file indicating the ids of genes in disease module and another file indicating the angle (theta) of each gene.
