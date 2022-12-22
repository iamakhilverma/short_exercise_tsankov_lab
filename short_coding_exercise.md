---
title: "short_coding_exercise"
author: "Akhil Kumar"
date: "2022-12-23"
output: 
  html_document:
    code_folding: hide
    keep_md: true
    toc: true
    number_sections: true
    toc_depth: 3
    toc_float: true
    theme: readable
    highlight: default
---



# Libraries

## Installations


```r
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("GEOquery")
```

## Loading


```r
# library(Seurat)
library(dplyr)
library(patchwork)

library(GEOquery)

library(R.utils)
```

# Data retrieval


```r
setwd("/Users/iamakhilverma/Desktop/short_exercise_tsankov_lab")
getwd()
```

```
## [1] "/Users/iamakhilverma/Desktop/short_exercise_tsankov_lab"
```


GEO accession: `GSE125449`
A total of 19 tumors were profiled. Set 1 contains scRNA-seq data of twelve samples, i.e., S16_P10_LCP18, S02_P01_LCP21, S10_P05_LCP23, S09_P04_LCP25, S08_P03_LCP26, S07_P02_LCP28, S11_P06_LCP29, S12_P07_LCP30, S20_P12_LCP35, S21_P13_LCP37, S15_P09_LCP38, and S19_P11_LCP39. Set 2 includes scRNA-seq data of seven samples, i.e., S351_P10_LCP34, S355_P13_LCP42, S358_P16_LCP46, S305_P06_LCP56, S300_P02_LCP60, 364_P21_LCP65, and S365_P22_LCP66. Detailed information can be found in samples.txt file of each Set.

`GPL18573`	Illumina NextSeq 500 (Homo sapiens)
`GPL20301`	Illumina HiSeq 4000 (Homo sapiens)

We care about the data contained in Set 1, i.e., `GSE125449-GPL18573_series_matrix.txt.gz`.


```r
# gse <- getGEO('GSE125449',GSEMatrix=TRUE)
# show(gse)
```

```r
# options(timeout = max(300, getOption("timeout")))
# options(download.file.method.GEOquery = "wget")
```


```r
if (!(file.exists("/Users/iamakhilverma/Desktop/short_exercise_tsankov_lab/GSE125449/"))){
  filePaths <- getGEOSuppFiles('GSE125449')
  filePaths
}
```

```r
unzipper <- function(filename) {
  gunzip(filename, destname=paste(
  '/Users/iamakhilverma/Desktop/short_exercise_tsankov_lab/GSE125449/gse125449_set1/',
  gsub("^GSE125449_Set1_", "", tail(strsplit(gsub("[.]gz$", "", filename), '/')[[1]], n=1)),
  sep=''
), remove=FALSE)
}
```



```r
if (!(file.exists("/Users/iamakhilverma/Desktop/short_exercise_tsankov_lab/GSE125449/gse125449_set1"))){
  filenames <- list.files(path="/Users/iamakhilverma/Desktop/short_exercise_tsankov_lab/GSE125449/", pattern="GSE125449_Set1_[A-Za-z0-9.]+.gz", full.names=TRUE)
  filenames
  sapply(filenames, unzipper)
}
```


```r
# livcan.data <- 
```


# R Markdown

For Markdown setup options, see <https://bookdown.org/yihui/rmarkdown/html-document.html>. 

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
summary(cars)
```

```
##      speed           dist       
##  Min.   : 4.0   Min.   :  2.00  
##  1st Qu.:12.0   1st Qu.: 26.00  
##  Median :15.0   Median : 36.00  
##  Mean   :15.4   Mean   : 42.98  
##  3rd Qu.:19.0   3rd Qu.: 56.00  
##  Max.   :25.0   Max.   :120.00
```

# Including Plots

You can also embed plots, for example:

![](short_coding_exercise_files/figure-html/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.


# Session Information


```r
sessionInfo(package=NULL)
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] R.utils_2.12.2      R.oo_1.25.0         R.methodsS3_1.8.2  
## [4] GEOquery_2.62.2     Biobase_2.54.0      BiocGenerics_0.40.0
## [7] patchwork_1.1.2     dplyr_1.0.10       
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.0  xfun_0.36         bslib_0.4.2       purrr_1.0.0      
##  [5] colorspace_2.0-3  vctrs_0.5.1       generics_0.1.3    htmltools_0.5.4  
##  [9] yaml_2.3.6        utf8_1.2.2        rlang_1.0.6       jquerylib_0.1.4  
## [13] pillar_1.8.1      glue_1.6.2        DBI_1.1.3         lifecycle_1.0.3  
## [17] stringr_1.5.0     munsell_0.5.0     gtable_0.3.1      evaluate_0.19    
## [21] knitr_1.41        tzdb_0.3.0        fastmap_1.1.0     fansi_1.0.3      
## [25] highr_0.10        readr_2.1.3       scales_1.2.1      cachem_1.0.6     
## [29] limma_3.50.3      jsonlite_1.8.4    ggplot2_3.4.0     hms_1.1.2        
## [33] digest_0.6.31     stringi_1.7.8     grid_4.1.0        cli_3.5.0        
## [37] tools_4.1.0       magrittr_2.0.3    sass_0.4.4        tibble_3.1.8     
## [41] tidyr_1.2.1       pkgconfig_2.0.3   ellipsis_0.3.2    data.table_1.14.6
## [45] xml2_1.3.3        assertthat_0.2.1  rmarkdown_2.19    rstudioapi_0.14  
## [49] R6_2.5.1          compiler_4.1.0
```

