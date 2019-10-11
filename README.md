
# Segment clustering reference to other segment
Cluster a segment(domain) by superposing different part(domain) of a protein
structure. 

## Requirement:
* mod9.21
* python2 or python3
* Scikit-learn

## Experiment:
* Required distances were computed using [mod9.21](https://salilab.org/modeller/9.21/release.html). 
* Subsequent clustering was performed using sklearn.cluster.DBSCAN
  (Density-based spatial clustering of applications with noise ). 
* Precomputed distance matrix using mod9.21 was inputted in DBSCAN.

## Test:
Clusters computed using DBSCAN were further matched with Gromacs for 
testing purpose. Gromacs implements single linkage clustering algorithm by default. 
Clusters produced by both methods are statistically indistinguishable.

## Usage:
```
$ python ClusterDistMatrix.py -i ./unittest/rmsd.txt.tgz
```

## Reference:
* DBSCAN: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
* Gromacs: http://manual.gromacs.org/documentation/2018/onlinehelp/gmx-cluster.html
