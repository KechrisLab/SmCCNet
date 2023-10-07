# SmCCNet
SmCCNet is a canonical correlation analysis-based method for discovering (quantitative) trait(phenotype)-specific multi-omics networks. The method allows the integration of single/multiple data types and a quantitative/binary phenotypic trait. It incorporates a feature subsampling scheme to improve the robustness of the canonical weights. 

To install and use the package, you may download the directory or follow the instructions below.
```{r, install-and-example}
# Install package
if (!require("devtools")) install.packages("devtools")
devtools::install_github("KechrisLab/SmCCNet")

# Load package
library(SmCCNet)
```

In the **vignettes** folder, users can find documentation that illustrates how to implement SmCCNet with an example data set. The data file is included under the **data** folder. Details on all the functions included in the package are documented in the package manual under the **package** folder. Users may directly download the package tarball for all functions and example data, which is also under the **package** folder.

Visualization can be done in both ways: ShinyApp or Cytoscape. After the subnetwork .Rdata is obtained, the user can upload it to the ShinyApp (https://smccnet.shinyapps.io/smccnetnetwork/) and tune visualization parameters, it will return the subnetwork visualization and the user will be able to download it as a .pdf file.

Please report any issues at https://github.com/KechrisLab/SmCCNet/issues.


