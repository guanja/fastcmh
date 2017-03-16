data="../../sampledata/data.txt"
label="../../sampledata/label.txt"
cov="../../sampledata/cov.txt"
alpha=0.05
Lmax=0
outputFolder="./output"
basename="/fastcmh"
basefilename=$outputFolder$basename
postprocessingfolder="../postprocessing"
pval="/pval"
pvalfile=$outputFolder$pval

mkdir $outputFolder

#run FastCMH
#./significant_interval_search_meta_cmh $data $label $cov $alpha $Lmax $basefilename -postprocessingfolder $postprocessingfolder -pval_file $pvalfile
./significant_interval_search_meta_cmh $data $label $cov $alpha $Lmax $basefilename -postprocessingfolder
# ./significant_interval_search_meta_cmh $data $label $cov $alpha $Lmax $basefilename 
