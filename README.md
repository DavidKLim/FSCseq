
# FSCseq

FSCseq is an R package for simultaneous feature selection and clustering
of RNA-seq gene expression data. It can also correct for differences in
sequencing depth using size factors from `DESeq2` ([Love et
al, 2014](https://doi.org/10.1186/s13059-014-0550-8)), as well as for
covariates such as batch. The main application is in delineating tumor
subtypes, but `FSCseq` can be used for other applications involving
discovery of subpopulations and identification of significant features.
Code to replicate the results from the FSCseq paper is available at
<https://github.com/DavidKLim/FSCseqPaper>.

## Installation

You can install the released version of FSCseq from this repository
with:

``` r
devtools::install_github("DavidKLim/FSCseq")
```

## Example

This example gives a brief overview of how to run `FSCseq` analysis with
simulated data. We show how to use a real RNA-seq read count dataset
instead of simulated data. The extension to a real dataset is
straight-forward

### Step 1a: Simulating data

To simulate data, use the `simulateData` function available in `FSCseq`.
Read count expression will be simulated from a finite mixture of
negative binomials. For ease of use, simulated data will be saved
automatically in the `save_dir` directory, and will be saved in object
`sim.dat`. In this example, one dataset is simulated (`nsims`) with
10000 genes (`G`) and 50 samples (`n`) from 1 batch (`B`) with 2
underlying clusters (`K`), baseline
![\\log\_2](https://latex.codecogs.com/png.latex?%5Clog_2 "\\log_2")
mean of 12 (`beta0`), and overdispersion of 0.35 (`phi`).

``` r
B=1; g=10000; K=2; n=50; LFCg=1; pDEg=0.05; beta0=12; phi0=0.35
set.seed(999)
sim.dat = FSCseq::simulateData(B=B, g=g, K=K, n=n, LFCg=LFCg, pDEg=pDEg,
             beta0=beta0, phi0=phi0, nsims=1, save_file=FALSE)[[1]]
# for save_file=TRUE, can input custom save_dir and save_prefix for parallelization of downstream analyses
```

The `simulateData` function outputs a list of length `nsims` with a
`sim.dat` list object for each simulation. In the above example, we
subset to just the 1st element (since `nsims=1`) of the output list, but
for `nsims>1`, other simulated datasets can be accessed by changing the
index value. The contents of each `sim.dat` object can be accessed, and
counts and true cluster labels can be extracted as follows:

``` r
str(sim.dat)
#> List of 9
#>  $ cts       : num [1:10000, 1:50] 7494 19478 3728 3721 10575 ...
#>  $ cts_pred  : num [1:10000, 1:25] 3436 12903 2449 4670 2733 ...
#>  $ cls       : int [1:50] 1 1 1 1 2 1 1 2 2 1 ...
#>  $ cls_pred  : int [1:25] 1 2 1 1 2 2 2 1 2 2 ...
#>  $ batch     : num [1:50] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ SF        : num [1:50] 0.93 0.672 1.199 1.068 0.931 ...
#>  $ SF_pred   : num [1:25] 0.971 0.839 1.436 1.092 0.983 ...
#>  $ DEg_ID    : logi [1:10000] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ sim_params:List of 18
#>   ..$ K       : num 2
#>   ..$ B       : num 1
#>   ..$ g       : num 10000
#>   ..$ n       : num 50
#>   ..$ n_pred  : num 25
#>   ..$ pK      : num [1:2] 0.5 0.5
#>   ..$ pB      : num 1
#>   ..$ LFCg    : num 1
#>   ..$ pDEg    : num 0.05
#>   ..$ sigma_g : num 0.1
#>   ..$ LFCb    : num 0
#>   ..$ pDEb    : num 0.5
#>   ..$ sigma_b : num 0
#>   ..$ beta    : num [1:10000, 1:2] 12 12 12 12 12 12 12 12 12 12 ...
#>   ..$ phi     : num [1:10000] 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35 ...
#>   ..$ disp    : chr "gene"
#>   ..$ LFCg_mat: num [1:10000, 1:2] 1 1 0 1 1 1 1 1 1 1 ...
#>   ..$ DEb_ID  : logi [1:10000] FALSE TRUE FALSE TRUE TRUE TRUE ...
cts=sim.dat$cts; true_cls=sim.dat$cls
```

### Step 1b: Analyzing custom data

To perform analysis on your own data, download the read counts and load
it. This example shows how to acquire the TCGA Breast Cancer Dataset
available on the NCI [GDCPortal](https://portal.gdc.cancer.gov/), using
the `TCGAbiolinks` package. Warning: this query contains 1215 files with
a total of about 1.84 GB, and will take a long time to download:

``` r
# library(devtools)
# devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
query1 = GDCquery(project="TCGA-BRCA",
                data.category = "Gene expression",
                data.type = "Gene expression quantification",
                platform = "Illumina HiSeq",
                file.type  = "results",
                experimental.strategy = "RNA-Seq",
                legacy = TRUE)
GDCdownload(query1)
GDCprepare(query = query1, save = TRUE, save.filename = "TCGA_BRCA_exp.rda")
```

Then, read the saved data into the R environment

``` r
load(file="TCGA_BRCA_exp.rda")
library(SummarizedExperiment)
cts <- round(assay(data),0)
cts <- cts[!duplicated(cts[,1:ncol(cts)]),]
anno <- colData(data)@listData
```

Optionally, you may want to pre-filter out genes with low FPKM values.
Subtype information for the TCGA BRCA dataset used in our paper can be
obtained using the `TCGAquery_subtype()` function, and used as the true
cluster labels. These true cluster labels are optional, but useful to
track diagnostics in FSCseq:

``` r
BRCA_tab = TCGAquery_subtype("BRCA")
match_ids = match(anno$patient,BRCA_tab$patient) # match patients
anno$subtypes=BRCA_tab$BRCA_Subtype_PAM50[match_ids]

true_cls = as.numeric(as.factor(anno$subtypes))
```

Then, you can proceed with the subsequent steps with the `cts` matrix
and true cluster labels `true_cls`, as in the simulated data. Details of
the processing steps and analyses on the TCGA BRCA dataset performed in
our paper can be found [here](https://github.com/DavidKLim/FSCseqPaper).
In the subsequent steps, we walk through just the simulated dataset, but
analysis can be done on your own data using the same steps.

### Step 2: Performing clustering and feature selection

Input the simulated (or custom) `cts` matrix into `FSCseq_workflow`.
Default search grids for tuning parameters are preset. For brevity of
illustration, we go through FSCseq analysis on the previously simulated
dataset with a much smaller grid of values of tuning parameters (takes
about 7-8 minutes). Note that `dir_name` should be unique, in order to
avoid utilizing saved results from a previously analyzed dataset:

``` r
cts=sim.dat$cts; true_cls=sim.dat$cls
t0 = as.numeric(Sys.time())
FSCseq_results = FSCseq::FSCseq_workflow(cts=cts,K_search=c(2:3),lambda_search=c(1.0, 1.5),
                                         alpha_search=c(0.1, 0.2),dir_name="~/test/Saved_Results")
#> Warning in FSCseq::FSCseq_workflow(cts = cts, K_search = c(2:3), lambda_search = c(1, : No input batch. Assuming all samples from same batch
#> converting counts to integer mode
t1 = as.numeric(Sys.time())
print(paste("time elapsed:",t1-t0))
#> [1] "time elapsed: 460.350342035294"
```

Note that we did not simulate batch in this case. If batch was
simulated, an additional argument `batch=...` can be input to adjust for
these batch effects

### Step 3: Summarizing and visualizing results

We can now summarize our clustering results. `FSCseq_workflow` outputs
the processed data after pre-filtering and normalizing for differences
in sequencing depth, as well as the results from FSCseq analysis. Store
included genes in `idx` to compare FSCseq results `res` with simulated
data:

``` r
res = FSCseq_results$results; processed.dat = FSCseq_results$processed.dat
idx = processed.dat$idx  # IDs of genes that are included in analysis after pre-filtering
library(mclust)
#> Package 'mclust' version 5.4.5
#> Type 'citation("mclust")' for citing this R package in publications.
print(paste("True K:",K,"Optimal K:",length(unique(res$cls))))
#> [1] "True K: 2 Optimal K: 2"
print(paste("ARI:",adjustedRandIndex(true_cls,res$cls)))
#> [1] "ARI: 1"
table(true_cls,res$cls)
#>         
#> true_cls  1  2
#>        1 23  0
#>        2  0 27
```

We can also show the true positive rate (TPR) and false positive rate
(FPR) of feature selection in discovering cluster-discriminatory genes:

``` r
true_disc = sim.dat$DEg_ID[idx];
FSC_disc = res$discriminatory

print(paste("TPR: ",sum(true_disc & FSC_disc)/sum(true_disc)))
#> [1] "TPR:  0.746231155778894"
print(paste("FPR: ",sum(!true_disc & FSC_disc)/sum(!true_disc)))
#> [1] "FPR:  0.00217296827466319"
```

We can visualize the expression patterns by plotting a heatmap, with
column annotations denoting cluster membership (red/black)

``` r
norm_y = processed.dat$norm_y
heatmap(log(norm_y[sim.dat$DEg_ID,]+0.1),scale="row",ColSideColors = as.character(res$cls),xlab="Samples",ylab="Genes",main="Heatmap of cluster-discriminatory genes")
```

<img src="man/figures/README-FSCseqHM-1.png" width="100%" />

### Step 4 (optional): Predicting on new data

`simulateData` additionally simulates a test set with the same simulated
parameters, in order to perform prediction after fitting the FSCseq
model on the training set. Input the FSCseq fitted object `res$fit` into
`FSCseq_predict_workflow`, along with the count matrix of the test set.
The count matrix of the training set is also required to use as a
pseudo-reference for the calculation of size factors in the test set.
Input `idx` to narrow down list of genes to those included in FSCseq
analyses.

``` r
cts_pred = sim.dat$cts_pred
true_cls_pred = sim.dat$cls_pred
fit_pred = FSCseq::FSCseq_predict_workflow(fit=res$fit,cts=cts,cts_pred=cts_pred,idx=idx)
#> converting counts to integer mode
```

``` r
res_pred = fit_pred$results
print(paste("pARI: ",adjustedRandIndex(true_cls_pred,res_pred$clusters)))
#> [1] "pARI:  1"
```

This anaysis can be easily generalized to real data by replacing
`cts_pred` with a separate test dataset, after fitting the FSCseq model
on the training dataset.
