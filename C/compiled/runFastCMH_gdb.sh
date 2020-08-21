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

gdb --args ./significant_interval_search_meta_cmh $data $label $cov $alpha $L_max $basefilename -postprocessing_folder $postprocessing_folder -pval_file $pval_file
