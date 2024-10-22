module load gcc/10.2.0
module load R/R-4.2.0

tissue="ALL"
trait="FEV1FVC"
rsid="rs17168118"
outputcell="../results/plots"
cv2f_dir=""


cmd="Rscript cv2f_plot.R --tissue $tissue --trait $trait --rsid $rsid --outputcell $outputcell"
bsub -W 2440 -R "rusage[mem=20]" -e ${outputcell}/${tissue}_${trait}_${rsid}.err -o ${outputcell}/${tissue}_${trait}_${rsid}.out -n 1 "$cmd"

