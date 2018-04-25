# SmCCNet
A canonical correlation analysis based method for discovering (quantitative) trait-specific heterogeneous regulatory networks.

SmCCNet is a data integration method based on canonical correlation analysis (CCA). It allows the integration of multiple data types and a quantitaive phenotypic trait. This method also includes a subsampling scheme of features to improve the robustness of the canonical weights. 

Ther are three R scripts under the Code folder. The ModifiedPMA.R script contains the modified CCA functions, that were built based on the PMA package written by Daniella Witten. The SmCCNetSource.R script provides all necessary functions for performing SmCCNet. It also includes a plot function that utilizes the igraph package to allow visualization of identified module networks. The LocalCV.R script allows to run K-fold cross validations in parallel locally.

The ToyExample folder will include a small example to illustrate how to use SmCCNet. 
