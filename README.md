# SmCCNet
SmCCNet is a canonical correlation analysis based method for discovering (quantitative) trait-specific multi-omics networks. The method allows the integration of multiple data types and a quantitaive phenotypic trait. It incorporates a feature subsampling scheme of features to improve the robustness of the canonical weights. 

To install and use the package, you may download the directory or follow the instructions below.
```{r, install-and-example}
# Install package
if (!require("devtools")) install.packages("devtools")
devtools::install_github("KechrisLab/SmCCNet")

# Load package
library(SmCCNet)
```

In the **vignettes** folder, users can find a documentation that illustrates how to implement SmCCNet with an example data set. The data file is included under the **data** folder. Details on all the functions included in the package are documented in the package manual under the **package** folder. Users may directly download the package tarball for all functions and example data, which is also under the **package** folder.

Please report any issues at https://github.com/KechrisLab/SmCCNet/issues.
