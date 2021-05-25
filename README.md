 APApathway
=============================

This R package is to obtain the dynamic APA usage difference score of the pathway path through the gene expression matrix, pathway pathway data, and PPI network.

Installing APApathway
=============
Mandatory 
---------

* R (>=3.6) is recommended.

Required R Packages
---------
[Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)

Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("BMILAB/APApathway")
library(APApathway)
```
reparations
====================

Gene expression count matrix
---------
The input to scNPF is matrix of gene expression count. The rows correspond to genes and the columns correspond to examples. The disease-related seed gene APA score we used in this study comes from the data set in the XIA paper. Including BRCA, BLCA, LUSC, LUAD, HNSC, UCEC and KIRC seven kinds of cancer disease related seed gene APA score.

```
##Loading gene expression count matrix
data("LUCS")
dim(LUCS)
```

Gene-gene interaction network
---------

The gene-gene interaction network (adjacency matrix) is used to randomly walk the obtained APA scores, so that every gene in the interaction network gets a score.If users want, they can provide gene co-expression network from publicly available database or a specific gene-gene network established on your own method. In this package, we provided a human gene-gene interaction networks [String](https://doi.org/10.1093/nar/gks1094).  

Pathway data set
---------
After each gene in the interaction network has a corresponding APA score, we obtain the score of the pathway according to the corresponding gene in the pathway gene set. Get the correlation between pathway and disease according to the score.User can choose to use their own pathway data or our data provided in the package.we provided a pathway data set [KEGG](http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP:KEGG).

Using APApathway
=============
In order to facilitate user understanding, we use the provided example LUCS,string,KEGGto illustrate the standard analysis work-flow of APApathway. APApathway randomly walks the scores of the cell expression matrix on the interaction network, and then obtains the pathway scores based on the genes in the pathway data gene set. The GetScore function contains five parameters, x is the gene expression matrix, network is the interaction network, pathway is the pathway data gene set, and gamma is the trade-off between prior   information and network diffusion, governing the distance that a signal is allowed to diffuse through the network during smoothing.The specific value of gamma has little effect on the results of network propagation over a sizable range.The default value is 0.5. The select parameter is a threshold, and the absolute value of the score is less than this value as an invalid value and will not participate in the calculation.he default value is 0.005.

```
##Using String network ,KEGG,LUCS data to obtain pathway score.
data("KEGG")
data("LUCS")
data("string")
PathwayScore <- APApathway(x=LUCS.data,network=string,pathway=KEGG)
```
