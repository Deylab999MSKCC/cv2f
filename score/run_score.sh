module load gcc/10.2.0
module load R/R-4.2.0
conda activate cv2f

snpcell="../data"
featurecell="../data/feature_tables"
mafpath="..//MAF_features_Aug032022.txt"
bimpath="../data/1000G_BIMS_hg38/1000G.EUR.QC."
ldblockspath="../data/LAVA_LDblocks_published.txt"
outputcell="../results/score"
chrm=22

pos_file="positive_set.txt"
neg_file="negative_set.txt"
feature_file="baseline."
score_file="mwe_baseline.22.cv2f"
echo $pos_file $neg_file $metric_file
if [ ! -f $outputcell/$metric_file ]
then
cmd="Rscript run_cv2f_MAF5LDmatch_cv2f_scores.R --positive_set $snpcell/$pos_file --negative_set $snpcell/$neg_file --feature_file $featurecell/$feature_file --mafpath $mafpath  --bimpath $bimpath --ldblockspath $ldblockspath --output_cv2f $outputcell/$score_file"
bsub -W 360 -R "rusage[mem=10]" -e ../results/score/mwe_baseline.22.cv2f.err -o ../results/score/mwe_baseline.22.cv2f.out -n 1 "$cmd"
fi

