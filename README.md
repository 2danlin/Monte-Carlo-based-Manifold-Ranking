# Preview

This package contains the source codes for all the algorithms related to Manifold Ranking-based Image Retrieval.
Please cite our paper below for usage:
```
@inproceedings{lin2019first,
  title={First Index-Free Manifold Ranking-Based Image Retrieval with Output Bound},
  author={Lin, Dandan and Wei, Victor Junqiu and Wong, Raymond Chi-Wing},
  booktitle={2019 IEEE International Conference on Data Mining (ICDM)},
  pages={1216--1221},
  year={2019},
  organization={IEEE}
}
```

This package contains all source codes for 
1. MCMR (ours)
2. MCMR+ (ours)
3. Power (for ground truth)
4. Mogul (VLDB'14)
5. Mogul-E (VLDB'14)


# Preparation before Usage (taking COIL for example)

## a. Dataset preparation (a large 500-NN graph precomputed by using deep learning methods)
1. Make a folder under "dataset" : 
	
```
mkdir large_knn
```
   
2. Put the "500-NN graph" file of each dataset under the new folder
	
e.g., an image dataset called "coil" and its "500-NN graph" file is with the name "coil.txt" 
	
(Format of the file could be found in Appendix)

## b. Generate the real "k-NN" graph with a given value of parameter $k$ (the parameters of the command could be found in Appendix)

```
./MCMR -d coil -algo BUILD_KNN_GRAPH -k 5
```
	
(the output is two files under folder "dataset": one is called "coil_5.txt" and the other is called "sym_coil_5.txt")

## c. Generate the file of query nodes (the parameters of the command could be found in Appendix)
```
./MCMR -d coil_5 -algo GEN_QUERY -n 100
```

## d. Generate the ground truth for each query node
1. Make a folder under "real_mr":
```
mkdir coil_5
```
  
2. Generate the ground truth for this dataset (the parameters of the command could be found in Appendix)
```
./MCMR -d coil_5 -algo GEN_GROUND_TRUTH -n 50
```

## e. Results preparation
1. Make a folder under "estimated_mr":
```
mkdir coil_5
```

## f. Index preparation
1. Make a folder under "index":
```
mkdir coil_5
```

# Usage Step (MCMR)

## a. Compilation
```
./compile.sh
```

## b. Execution of query phase (the parameters of the command could be found in Appendix)
```
./MCMR -d coil_5 -algo MCMR -n 50 -k 5 -e 0.5 -a 0.99 -c 10000
```

## c. Output Files (Format of those files can be found in Appendix)
1. "MCMR_5.txt" which captures the average time and precision statistics of the algorithm


# Usage Step (MCMR+)

## a. Compilation
```
./compile.sh
```

## b. Execution of query phase (the parameters of the command could be found in Appendix)
```
./MCMR -d coil_5 -algo MCMRPlus -e 0.5 -a 0.99 -n 50 -k 5
```

## c. Output Files (Format of those files can be found in Appendix)
1. "MCMRPlus_5.txt" which captures the average time and precision statistics of the algorithm

# Usage Step (Power)

## a. Compilation
```
./compile.sh
```

## b. Execution of query phase (the parameters of the command could be found in Appendix)
```
./MCMR -d coil_5 -algo POWER -n 50 -k 5 -t 1000 -a 0.99
```

## d. Output Files of query phase (Format of those files can be found in Appendix)
1. "POWER_5.txt" which captures the average time and precision statistics of the algorithm



# Usage Step (Mogul)

## a. Compilation
```
./compile.sh
```

## b. Preprocessing phase (the output files can be found in Appendix)
```
./MCMR -d coil_5 -algo PRE_MOGUL
```

## c. Execution of query phase (the parameters of the command could be found in Appendix)
```
./MCMR -d coil_5 -algo MOGUL -n 50 -k 5 -a 0.99
```

## d. Output Files of query phase (Format of those files can be found in Appendix)
1. "Mogul_5.txt" which captures the average time and precision statistics of the algorithm

# Usage Step (Mogul-E)

## a. Compilation
```
./compile.sh
```

## b. Preprocessing phase (the output files can be found in Appendix)
```
./MCMR -d coil_5 -algo PRE_MOGUL_E
```

## c. Execution of query phase (the parameters of the command could be found in Appendix)
```
./MCMR -d coil_5 -algo MOGUL_E -n 50 -k 5 -a 0.99
```

## d. Output Files of query phase (Format of those files can be found in Appendix)
1. "Mogul_E_5.txt" which captures the average time and precision statistics of the algorithm


# Appendix A. Parameters of this program
------------------------------------
This program contains 8 parameters.
The general command of this program is shown below:
```
./MCMR [parameters]
```

1.-d <dataset> : the name of dataset 
> e.g., "-d coil" or "-d coil_5"

2.-algo <algorithm> : the name of algorithm (all the supported algorithms could be found in Appendix)
> e.g., "-algo MCMR"

3.-n <node_count> : no. of query nodes
> e.g., "-n 1" (the default value is 1)

4.-alpha <alpha> : the value of the constant parameter in Manifold Ranking
> e.g., "-alpha 0.2" (the default value is 0.2)

5.-k <no. of top-k nodes> : the number of nodes to be returned for top-k query
> e.g., "-k 5" (the default value is 1) 

6.-e <epsilon> : the value of relative error to be guaranteed (for algorithms MCMR and MCMR+ only)
> e.g., "-e 0.5" (the default value is 0.5)

7.-t <no. of iterations> : the number of iterations used in Power (for algorithm POWER only)
> e.g., "-t 100" (the default value is 1000)

8.-c <no. of random walks> : the parameter for tuning the number of random walks used in MCMR only
> e.g., "-c 100000" (the default value is 100)
 
# Appendix B. Algorithms supported by this program

## 1. "BUILD_KNN_GRAPH": 

the algorithm for generating real k-nn graph by exploiting the precomputed large 500-NN graph where $k$ is given in the query phase

## 2. "GEN_QUERY" : 

the algorithm for generating the query nodes

## 3. "GEN_GROUND_TRUTH" : 

the algorithm for generating the exact MR scores for the query nodes

## 4. "MCMR" : 

the algorithm for top-k MR query (proposed by us)

## 5. "MCMRPlus" : 

the improved algorithm for top-k MR query (proposed by us)

## 6. "Power" : 

the iterative algorithm for top-k MR query (baseline)

## 7. "MOGUL" : 

the existing most efficient approximate algorithm for top-k MR query (VLDB'14)

## 8. "MOGUL_E" : 

the existing most efficient exact algorithm for top-k MR query (VLDB'14)



# Appendix C. Format of dataset file (taking the large 500-NN graph"coil.txt" or the real knn graph "coil_5.txt" for example)

## 1. The first line in the file is 
       - > no. of nodes in this graph

## 2. The following lines after the first line is
      
      - the edges of this graph which are represented by the following format:
      - 
        ```
	<From-node> <To-node> <weight> (separated by space)
	```



# Appendix D. Format of "MCMR_5.txt", "MCMRPlus_5.txt","POWER_5.txt", "MOGUL_5.txt" and "MOGUL_E_5.txt"

## the information of each line is given in the following:

1. "k: 5" --> the value of $k$ (i.e., the top-k query) 

2. "avg precision: " --> average precision of top-k query of this algorithm

3. "avg time: " -->average query time cost of the top-k query of this algorithm

# Appendix E. The output files of the preprocessing phase of "Mogul"

## 1. "cluster.txt": the results of clustering in their paper
	
	- 1.1. The first line in the file is 
 	- > no. of clusters in this graph

	- 1.2. The following lines after the first line are represented by the following format:
        - 
	  ```
	  <cluster-id> <node-id> (separated by space)
	  ```

## 2. "reordered_cluster.txt": the results of reordered clustering in their paper
	
	- 2.1. The first line in the file is 
 	- > no. of clusters in this graph

	- 2.2. The following lines after the first line are represented by the following format:
   	- 
	```
	<cluster-id> <node-id> (separated by space)
	```
	
## 3. "permutation.txt":  the permutation matrix based on the reordered clusters
	- 3.1. The first line in the file is
	- > no. of nodes in this graph

	- 3.2 Each line after the first line is represented by the following format:
	- 
	  ```
	  <new-node-id> <original-node-id> (separated by space)
	  ```

## 4. "incom_chol_D.txt": the decomposed matrices "D" in their paper
	- 4.1 The $i$-th line represents the value of D_ii
	- 
	  ```
	  <value>
	  ```

## 5. "incom_chol_L.txt": the decomposed matrices "L" in their paper
	- 5.1 Each line represents the value of L_ij where i is the id of row and j is the id of column
	- 
	```
	<row> <col> <value> (separated by space)
	```
	
## 6. "preprocessing.time" : the time file for this preprocessing
	- 6.1 "clustering time" : the time for doing the graph clustering
	- 6.2 "reorder time" : the time for reordering the clusters
	- 6.3 "permute time" : the time for getting the permutation matrix 
	- 6.4 "decomposition time" : the time for decomposing the matrix
	- 6.5 "total preprocessing time" :  the total preprocessing time of Mogul

# Appendix F. The output files of the preprocessing phase of "Mogul-E"

## 1. "cluster.txt": the results of clustering in their paper
	- 1.1. The first line in the file is 
 	- > no. of clusters in this graph

	- 1.2. The following lines after the first line are represented by the following format:
   	- 
	  ```
	  <cluster-id> <node-id> (separated by space)
	  ```

## 2. "reordered_cluster.txt": the results of reordered clustering in their paper
	- 2.1. The first line in the file is 
 	- >no. of clusters in this graph

	- 2.2. The following lines after the first line are represented by the following format:
   	- 
	 ```
	 <cluster-id> <node-id> (separated by space)
	 ```
	 
## 3. "permutation.txt":  the permutation matrix based on the reordered clusters
	
	- 3.1. The first line in the file is
	- > no. of nodes in this graph

	- 3.2 Each line after the first line is represented by the following format:
	- 
	 ```
	 <new-node-id> <original-node-id> (separated by space)
	 ```
	 
## 4. "chol_D.txt": the decomposed matrices "D" in their paper
	- 4.1 The $i$-th line represents the value of D_ii
	- 
	 ```
	 <value>
	 ```
	 
## 5. "chol_L.txt": the decomposed matrices "L" in their paper
	- 5.1 Each line represents the value of L_ij where i is the id of row and j is the id of column
	-
	 ```
	 <row> <col> <value> (separated by space)
	 ```
	 
## 6. "preprocessing_Mogul_E.time" : the time file for this preprocessing

	- 6.1 "clustering time" : the time for doing the graph clustering
	- 6.2 "reorder time" : the time for reordering the clusters
	- 6.3 "permute time" : the time for getting the permutation matrix 
	- 6.4 "decomposition time" : the time for decomposing the matrix
	- 6.5 "total preprocessing time" :  the total preprocessing time of Mogul


