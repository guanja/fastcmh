# Genetic Heterogeneity Discovery with FastCMH

Implementation in FastCMH algorithm in C and R. The R package `fastcmh` is [on CRAN](https://CRAN.R-project.org/package=fastcmh). The paper can be found [here](https://goo.gl/2QN2La).

## Installation

Either the R package can be used, or the C code can be used directly.

### R package

The easiest way to install the R package is via CRAN, by using the R console:

```
install.packages("Rcpp")    #dependency  
install.packages("bindata") #dependency  
install.packages("fastcmh")
```

or simply:

```
install.packages("fastcmh", dependencies=TRUE)
```

### Data format

Suppose you have a data set with 
* `n` samples (e.g. patients), 
* `L` features (e.g. SNPs) that have a binary encoding,
* each sample has two possible labels (e.g. case/control), 
* each sample one of `K` classes for the categorical covariate on which you wish to condition (e.g. country: Spain, France, Germany...). 

In order to run FastCMH, three files are needed:
1. `data.txt`: a file containing `L` rows and `n` columns, each (_i_, _j_ )th space-separated entry either `0` or `1` corresponding to the value of the binary feature for the _i_ th feature of the _j_ th sample.
2. `label.txt`: a file containing `n` rows, each row containing a single entry that is either `0` or `1`, where the _j_ th row gives the label for the _j_ th sample.
3. `cov.txt`: a file containing `K` rows, each row containing a single positive integer, where the _m_ th row indicates the number of samples that have the _m_ th value of the categorical covariate.

**Note:** The rows in `data.txt` and `label.txt` need to be ordered so that the first `n_1` rows have covariate class 1, the next `n_2` rows have covariate class 2, ..., the final `n_K` rows have covariate class `n_K`.

**Note:** The files need not have the default file names above, but for the rest of this description we shall assume that this is the case.


## Example: sample data

This example shows a minimal synthetic dataset in order to illustrate the format of the data. It can be generated in R using the package `fastcmh` using the following commands:

```
library(fastcmh)
makefastcmhdata(folder="./", L=20, n=50, K=2, tau1=5, taulength1=4, tau2=12, taulength2=4, seednum=3)
```

This will create three files `data.txt`, `label.txt`, `cov.txt` in the current working directory. As the code suggests, this dataset has `L=20` features, `n=50` samples and `K=2` classes for the categorical covariate. The true significant interval starts at `tau1=5` and has length `4` (i.e. the interval `[5,8]`), while a confounded interval is created starting at `tau2=12`, and also with length `4` (i.e. the interval `[12, 15]`) The random seed is set to be `3`.

**Note:** While the FastCMH algorithm will handle _any number_ `K` for the number of classes of the categorical covariate, the R script `makefastcmhdata` will only generate synthetic data for `K=2`.


The contents of `data.txt` are (`L=20` features/rows, `n=50` samples/columns):

```
0 1 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 1 1 1 0 0 1 0 0 0 1 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0
0 0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 1 0 0 0 0 0 0 0 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 0 1
0 0 0 0 0 0 0 1 0 1 0 0 0 1 1 0 0 0 1 0 0 1 1 0 0 1 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 0 0 1
1 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 1 0 1 1 0 1 0 0 1 0 0 1 0 0
0 0 1 1 0 1 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1 0 0 1 0 1 0
0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 1 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 1 0 1 0 0 1
0 0 1 1 0 1 0 1 1 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 1 0 1 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0
0 0 0 1 1 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 1 1 0 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0
1 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 1 0 1 1 1 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 1 0 1 0 1 0 0 0 0 0 0 1 1 0 1 1 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0 1 1 0
0 0 0 0 0 1 0 0 0 0 1 0 1 1 1 0 0 1 1 1 0 1 0 1 1 1 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1
0 0 0 0 1 0 0 1 0 1 0 0 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 1 1 1 0 1 1 0 0 0 0 0
1 1 0 0 0 1 1 0 1 0 0 1 0 1 0 1 1 0 0 0 0 0 1 0 0 0 0 1 1 0 1 1 1 0 1 0 0 1 0 0 0 0 0 1 0 1 0 1 1 0
0 1 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0 1 0 0
```

The `label.txt` file has `n=50` rows (only first ten are shown):

```
0
0
0
1
1
1
1
1
1
0
[etc, truncated]
```

The `cov.txt` file has `K=2` rows (each row with one entry, and the rows sum up to `n=50`):

```
22
28
```

The FastCMH algorithm can be run in R using the following code:

```
df <- runfastcmh(folder="./")
```

And the (filtered) significant intervals can be obtained simply using:

```
df$sig
```

which will return the following data frame:

```
  start end       pvalue
  1     5   8 0.0003525223
```

This concludes the description of the data format.


## Compiling and running the C code

To compile the C version of the code, enter the `C` folder and run `make`:

```
cd C
make
```

In the `C` folder is a shell script which shows how to run the C code on the sample data:

```
data="../../sampledata/data.txt"
label="../../sampledata/label.txt"
cov="../../sampledata/cov.txt"

alpha=0.05
L_max=0

outputfolder="./output/"

basefile="fastcmh"
basefilename=$outputfolder$basefile

postprocessing_folder="../postprocessing/"

pval="allpval.txt"
pval_file=$outputfolder$pval

mkdir -p $outputfolder

./significant_interval_search_meta_cmh $data $label $cov $alpha $L_max $basefilename -postprocessing_folder $postprocessing_folder -pval_file $pval_file
```


## Contact

Any questions can be directed to:  
* Felipe Llinares Lopez (C code): felipe.llinares [at] bsse.ethz.ch  
* Dean Bodenham (R package): dean.bodenham [at] bsse.ethz.ch 
