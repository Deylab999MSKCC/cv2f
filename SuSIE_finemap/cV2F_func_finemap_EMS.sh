module load gcc/10.2.0
module load R/R-4.2.0
TASKFILE=/data/deyk/GWAS/Finucane_2023/traitnames2.txt

IFS="
"

for line in `cat $TASKFILE | uniq`;
do
    annot=`echo $line | awk '{print $1}'`
    echo $annot
    cmd="Rscript cV2F_func_finemap_EMS.R --trait $annot"
    bsub -W 100 -R "rusage[mem=35]" -e chunks.err -o chunks.out -n 1 "$cmd"
done

