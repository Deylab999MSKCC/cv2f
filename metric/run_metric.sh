module load gcc/10.2.0
module load R/R-4.2.0
conda activate cv2f

snpcell="../data"
featurecell="../data/feature_tables"
mafpath="../data/MAF_features_Aug032022.txt"
bimpath="../data/1000G_BIMS_hg38/1000G.EUR.QC."
ldblockspath="../data/LAVA_LDblocks_published.txt"
outputcell="../results/comparisons_metric"


pos_train_file="positive_set.txt"
neg_train_file="negative_set.txt"
pos_test_file="positive_set.txt"
neg_test_file="negative_set.txt"
feature_file="baseline."
metric_file="mwe_baseline.metrics"
echo $pos_train_file $neg_train_file $pos_test_file $neg_test_file $metric_file
if [ ! -f $outputcell/$metric_file ]
then
cmd="Rscript run_cv2f_MAF5LDmatch_metrics.R --positive_set_train $snpcell/$pos_train_file --negative_set_train $snpcell/$neg_train_file --positive_set_test $snpcell/$pos_test_file --negative_set_test $snpcell/$neg_test_file --feature_file $featurecell/$feature_file --mafpath $mafpath  --bimpath $bimpath --ldblockspath $ldblockspath --output_metrics $outputcell/$metric_file"
bsub -W 360 -R "rusage[mem=10]" -e ../results/metric/mwe_baseline.metrics.cv2f.err -o ../results/metric/mwe_baseline.metrics.cv2f.out -n 1 "$cmd"
fi

