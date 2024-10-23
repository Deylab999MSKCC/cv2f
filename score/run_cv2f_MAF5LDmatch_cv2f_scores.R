#modified by Thahmina A. Ali

required_packages <- c("data.table", "xgboost", "pROC", "PRROC", "optparse")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')


library(data.table)
library(R.utils)
library(optparse)
library(xgboost)
library(pROC)
library(PRROC)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--positive_set", type="character", default = "../data/positive_set.txt",
              help="rsIDs of SNPs in the positive set (celltype relevant trait)"),
  make_option("--negative_set", type="character", default = "../data/negative_set.txt",
              help="rsIDs of SNPs in the negative set (non-celltype relevant traits)"),
  make_option("--feature_file", type="character", default = "../data/feature_tables/baseline.",
              help="Data frame of feature tables related to a cell type"),
  make_option("--bimpath", type="character", default = "../data/1000G_BIMS_hg38/1000G.EUR.QC.",
	      help="Path and prefix of the bimfile of the BIMFILE"),
  make_option("--mafpath", type="character", default = "../data/MAF_features_Aug032022.txt",
	      help="Path and prefix of the frequency file for all SNPs"),
  make_option("--ldblockspath", type="character", default = "../data/LAVA_LDblocks_published.txt", 
	      help="Path and prefix of the LD blocks file"),
  make_option("--output_cv2f", type="character", default = "../results/score/mwe_baseline.22.cv2f",
              help=" Path of the output celltype cV2F file"),
  make_option("--chrm", type="integer", default = 22,
              help=" chromosome number to produce scores for")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

# opt = list()
# opt$positive_set =  "/data/deyk/GWAS/UKBiobank/finemapping/Finemap_SNPsets/ALL_combined_PIP0_75.txt"
# opt$negative_set =  "/data/deyk/GWAS/UKBiobank/finemapping/Finemap_SNPsets/Notfinemapped_PIP0_01.txt"
# opt$feature_file = "/data/deyk/kushal/cV2F/data/Deliverables/Feature_Tables/baselineLD/feature_table_perchr."
# opt$mafpath = "/data/deyk/extras/MAF_features_Aug032022.txt"
# opt$bimpath = "/data/deyk/extras/1000G_BIMS_hg38/1000G.EUR.QC."
# opt$ldblockspath = "/data/deyk/extras/loci_files/LAVA_LDblocks_published.txt"

pre_positive_snps = read.table(opt$positive_set, header=F)[,1]
pre_negative_snps = read.table(opt$negative_set, header=F)[,1]
maf_tabb = data.frame(fread(opt$mafpath))
LDblocks_tabb = data.frame(fread(opt$ldblockspath))

##########################################################################################################

cat("Define positive and LD and MAF matched negative set of variants")
cat("\n")

bimpooltabb = c()
for(numchr in 1:22){
  bimdf = data.frame(fread(paste0(opt$bimpath, numchr, ".bim")))
  bimpooltabb = rbind(bimpooltabb, bimdf)
  cat("Processing SNPs from chr:", numchr, "\n")
}

colnames(bimpooltabb) = c("CHR", "SNP", "CM", "BP", "A1", "A2")

positive_snps = intersect(pre_positive_snps, maf_tabb$SNP)

pre_negative_snps = setdiff(intersect(pre_negative_snps, maf_tabb$SNP), pre_positive_snps)

negative_snps = c()
for(numl in 1:nrow(LDblocks_tabb)){
  snps_in_block = bimpooltabb$SNP[which(bimpooltabb$CHR == LDblocks_tabb$CHR[numl] &
                                  bimpooltabb$BP > LDblocks_tabb$START[numl] &
                                  bimpooltabb$BP < LDblocks_tabb$STOP[numl]) ]
  positives_in_block = intersect(positive_snps, snps_in_block)
  negatives_in_block = intersect(snps_in_block, pre_negative_snps)

  if(length(positives_in_block) > 0){
    mafs_positive_block = maf_tabb$MAF[match(positives_in_block, maf_tabb$SNP)]
    mafs_negative_block = maf_tabb$MAF[match(negatives_in_block, maf_tabb$SNP)]
    snps_temp = c()
    for(tt in 1:length(mafs_positive_block)){
      num_snps = pmin(5, length(mafs_negative_block))
      snps_temp = c(snps_temp, negatives_in_block[order(abs(mafs_positive_block[tt] - mafs_negative_block), decreasing = F)[1:num_snps]])
    }
    negative_snps = c(negative_snps, snps_temp)
  }
  cat("Matching for LD-block:", numl, "\n")
}

negative_snps = negative_snps[!is.na(negative_snps)]


positive_snps2 = c()
negative_snps2 = c()
pos_feature_tabb_list = c()
neg_feature_tabb_list = c()

for(numchr in 1:22){
  dff = data.frame(fread(paste0(opt$feature_file, numchr, ".txt")))
  pos_feature_tabb_list[[numchr]] = dff[match(intersect(positive_snps, dff$SNP), dff$SNP), -(1:4)]
  neg_feature_tabb_list[[numchr]] = dff[match(intersect(negative_snps, dff$SNP), dff$SNP), -(1:4)]
  positive_snps2 = c(positive_snps2, intersect(positive_snps, dff$SNP))
  negative_snps2 = c(negative_snps2, intersect(negative_snps, dff$SNP))
  cat("We are at chr:", numchr, "\n")
}

pos_feature_tabb = do.call(rbind, pos_feature_tabb_list)
neg_feature_tabb = do.call(rbind, neg_feature_tabb_list)
rownames(pos_feature_tabb) = positive_snps2
rownames(neg_feature_tabb) = negative_snps2


##########################################################################################################

combined_feature_tabb = rbind(pos_feature_tabb, neg_feature_tabb)
rownames(combined_feature_tabb) = c(positive_snps2, negative_snps2)
labels = c(rep(1, nrow(pos_feature_tabb)), rep(0, nrow(neg_feature_tabb)))
annotation_feature_tabb = combined_feature_tabb


cat("Create cV2F")
cat("\n")

cv2f_list = c()

numchr = opt$chrm

bimtabb = data.frame(fread(paste0(opt$bimpath, numchr, ".bim")))
temp_positive_snps = setdiff(positive_snps, bimtabb[,2])
temp_negative_snps = setdiff(negative_snps, bimtabb[,2])
dff = data.frame(fread(paste0(opt$feature_file, numchr, ".txt")))
test_data = dff[match(bimtabb[,2], dff$SNP), -(1:4)]

preds_tabb = c()
for(nboot in 1:10){
sample_positive_snps = temp_positive_snps
sample_negative_snps = sample(temp_negative_snps, length(temp_positive_snps), replace=F)

temp_pos_feature_tabb = annotation_feature_tabb[match(sample_positive_snps, rownames(annotation_feature_tabb)), ]
temp_neg_feature_tabb = annotation_feature_tabb[match(sample_negative_snps, rownames(annotation_feature_tabb)), ]
temp_combined_feature_tabb = rbind(temp_pos_feature_tabb, temp_neg_feature_tabb)
temp_labels = c(rep(1, nrow(temp_pos_feature_tabb)), rep(0, nrow(temp_neg_feature_tabb)))

training_data = temp_combined_feature_tabb
training_labels = temp_labels

NUMITER_XGBoost=1000
bstSparse <-  xgboost(data = as.matrix(training_data),
                      label = training_labels,
                      n_estimators = c(25, 40, 50),
                      max_depth = c(10, 15, 25),
                      learning_rate = 0.05,
                      gamma = 10,
                      min_child_weight = c(6, 8, 10),
                      nthread = 2,
                      scale_pos_weight = 1,
                      subsample = c(0.6, 0.8, 1),
                      nrounds = NUMITER_XGBoost,
                      objective = "binary:logistic",
                      eval_metric = "auc")

predict_labels_sample = predict(bstSparse, as.matrix(test_data))
preds_tabb = cbind(preds_tabb, predict_labels_sample)
}
predict_labels = rowMeans(preds_tabb)

cv2f_list[[numchr]] = cbind.data.frame(bimtabb[,1], bimtabb[,4], bimtabb[,2], bimtabb[,3], predict_labels)
cat("We are at chr:", numchr, "\n")


cv2f_df = do.call(rbind, cv2f_list)
colnames(cv2f_df) = c("CHR", "BP", "SNP", "CM", "cV2F")
fwrite(cv2f_df, file = paste0(opt$output_cv2f), row.names=F, col.names=T, sep = "\t", quote=F)
