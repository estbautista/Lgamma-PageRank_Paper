# <img src="https://latex.codecogs.com/svg.latex?\Huge&space;L^\gamma" title="" />-PageRank for Semi-Supervised Learning

This repository complements our manuscript entitled <img src="https://latex.codecogs.com/svg.latex?\normalsize&space;L^\gamma" title="" />-PageRank for Semi-Supervised Learning by Esteban Bautista, Patrice Abry and Paulo Gonçalves (Submitted to *Applied Network Science*). The repository provides both the results reported in the paper and the code to obtain them.

## Tools
Codes are written in [MATLAB](https://fr.mathworks.com).

The studied datasets are (as cited on the article): [MNIST](http://yann.lecun.com/exdb/mnist/) [26], [Gender images](http://cmp.felk.cvut.cz/~spacelib/faces/) [27], [BBC articles](http://mlg.ucd.ie/datasets/bbc.html) [28], [Phoneme](https://www.openml.org/d/1489) [29]. They are all included in the [Datasets folder](https://github.com/estbautista/Lgamma-PageRank_Paper/tree/master/Datasets).

## Experiments
#### Validation of Algorithm 1 in estimation of the optimal <img src="https://latex.codecogs.com/svg.latex?\Large&space;\gamma" title="" />
The [Algorithm_evaluation folder](https://github.com/estbautista/Lgamma-PageRank_Paper/tree/master/Algorithm_evaluation) contains:
* the evaluation of Algorithm 1 on the MNIST (reported in Table 1) 
* the code to generate Figure 2

#### Experiments on the Planted partition
The [PlantedPartition_experiment folder](https://github.com/estbautista/Lgamma-PageRank_Paper/tree/master/PlantedPartition_experiment) contains:
* the assessment of <img src="https://latex.codecogs.com/svg.latex?\normalsize&space;L^\gamma" title="" />-PageRank on the Planted Partition when partitions are retrieved by means of the sweep-cut (reported in Figure 3)
* the code to generate Figure 3

#### Experiments on real world datasets
The [RealWorldData_experiment folder](https://github.com/estbautista/Lgamma-PageRank_Paper/tree/master/RealWorldData_experiment) contains the assessment of Algorithm 1 and <img src="https://latex.codecogs.com/svg.latex?\normalsize&space;L^\gamma" title="" />-PageRank (with partitions via the sweep-cut) on the classification of real world datasets (both results reported in Table 2).


#### Experiments on unbalanced labeled data

The [UnbalancedLabels_experiment folder](https://github.com/estbautista/Lgamma-PageRank_Paper/tree/master/UnbalancedLabels_experiment) contains the performance assessment (multi-class approach) of <img src="https://latex.codecogs.com/svg.latex?\normalsize&space;L^\gamma" title="" />-PageRank in the presence of unbalanced labeled data (reported in Table 3).

###### Questions
Please contact [esteban.bautista-ruiz@ens-lyon.fr](mailto:esteban.bautista-ruiz@ens-lyon.fr) for any questions or problems.
