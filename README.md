
# Factorization-based Inferred Expression Single-cell Transcriptomic Analysis (FIESTA)

1. In order to install this package: 

``` {r}
library(devtools)
install_github("elnazmirzaei/Fiesta")
library(Fiesta)
library(NMF)
library(mixtools)
library(fitdistrplus)
```

2. After installation you need to prepare your gene-expression dataset as a matrix in R such that rows represent genes and columns represent cells. Then, pick whether to use WNMF or sNMF. Default is sNMF. Based on the size of your data this step might take a while.

``` {r}
A = GeneExpressionData
A_imputed = impute( A , method="sNMF" )
```


3. Now it is time for scaling.

``` {r}
A_imputed_scaled = scale( A , A_imputed )
```

4. And at the end run the thresholding step.

``` {r}
A_imputed_scaled_thresholded = thresholding ( A , A_imputed_scaled )
```

#### Voila! Now you have a properly imputed gene-expression data.
