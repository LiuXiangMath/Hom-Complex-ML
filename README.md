Hom Complex based machine learning for protein-protein interactions
====
    This manual is for the code implementation of paper "Hom complex based machine learning (HCML) for protein-protein binding affinity prediction"
    
****

# Software configuration
---
        Platform: Python>=3.6
        Packages needed: math, numpy>=1.18.1, scipy>=1.4.1, scikit-learn>=0.22.1

# Flow of HCML model
---
Protein-protein interaction complex  -->  Hom complex representation  -->   Feature generation  -->  Machine learning 
# Details about each step

## Protein-protein interaction complex
Before the representatoin, protein-protein complex needs to generated. More specifically, the atom coordinates are needed.

## Hom complex representatoin
For each protein-protein complex, eight atom sites are generated as in the paper. For each atom site, an alpha complex is generated. Then, the one-skeleton of the Alpha complex is extracted to form a filtration of graphs. Based on this filtration of graphs G, a sequence of nested Hom complexes $Hom(Z_2,G)$ are generated. The key functions are as follows:
```python
def get_alpha_one_skeleton(point,filtration)
    # this functoin get the one-skeleton graphs of the alpha complex from the point. 
    # parameter "point" is for the alpha complex generatoin
    # parameter "filtration" is the maximal edge we use in the alpha complex

def get_hom_complex_from_graph_new(V,E)
    # this function construct the Hom complex Hom(Z_2,G) from the graph G.
    # parameter "V" is the vertex-set of the graph
    # parameter "E" is a list of edges of the graph, each element of this list consists of three components: 
    # V1,V2, d representing that V1 and V2 have an edge and the length of this edge is d. 

def get_one_skeleton_hom_complex_from_graph(V,E)
    # if you only want the one-skeleton of the Hom complex, this function is a faster one
    # the parameters are same with the function above.

```
## Feature generation
Topological feature and auxiliary feature are used in our model. For the topological feature, persistent homology and persistent Euler characteristic are considered. The key functions are as follows:
```python
def zero_homology_of_a_complex_to_file(complex1,filename)
    # This functoin compute the zero dimensional persistent homology of the Hom complex, and write to file
    # parameter "complex1" is the Hom complex to be dealed with
    # parameter "filename" is the filename path used to store the persistent homology

def alpha_simplex_feature_to_file(typ,start,end,filtration,grid_size)
    # this function compute the persistent Euler characteristic
    # parameter "typ" needs to be set "euler"
    # parameter "start" and "end" need to be set 0, 645 and 0, 1131 for datasets S645 and SKEMPI S1131 respectively.
    # parameter "filtration" is 5 in our model
    # parameter "grid_size" is 0.1 in our model
```
## Machine learning
GradientBoostingRegressor is used to do the ten fold cross-validation in our model. Detailed parameters of GBT can be found in our paper. The corresponding python script is "10_cross_validation.py"

