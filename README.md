# SmCCNet
A canonical correlation analysis based method for discovering (quantitative) trait-specific heterogeneous regulatory networks.

SmCCNet is a data integration method based on canonical correlation analysis (CCA). It allows the integration of multiple data types and a quantitaive phenotypic trait. This method also includes a subsampling scheme of features to improve the robustness of the canonical weights. 

The main CCA function was build upon the PMA package written by Daniella Witten. The ModifiedPMA.R script contains the modified functions. The SmCCNetSource.R script provides all necessary functions for performing SmCCNet. It also includes a function based on the igraph package to allow visualization of identified modules. 
