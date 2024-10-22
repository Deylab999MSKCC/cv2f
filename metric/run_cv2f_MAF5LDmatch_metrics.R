#modified by Thahmina A. Ali

required_packages <- c("data.table", "xgboost", "pROC", "PRROC", "optparse", "SHAPforxgboost")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')


library(data.table)
library(R.utils)
library(optparse)
library(xgboost)
library(pROC)
library(PRROC)
require(SHAPforxgboost)
require(plyr)
require(dplyr)
library(stringr)


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
  make_option("--bimpath", type="character", help="Path and prefix of the bimfile of the BIMFILE"),
  make_option("--mafpath", type="character", help="Path and prefix of the frequency file for all SNPs"),
  make_option("--ldblockspath", type="character", help="Path and prefix of the LD blocks file"),
  make_option("--output_metrics", type="character", default = "../data/out",
              help=" Path of the output AUROC/AUPRC metrics file")
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


cat("Compute performance metrics to assess classification accuracy")
cat("\n")

all_shap = list()
all_training = list()

aucvec = c()
prcvec = c()
for(num_iter in 1:5){
  idx1 = sample(1:nrow(pos_feature_tabb), floor(0.70*nrow(pos_feature_tabb)), replace = F)
  idx2 = sample(1:nrow(neg_feature_tabb), floor(0.70*nrow(pos_feature_tabb)), replace =	F)
  training_idx = c(idx1, nrow(pos_feature_tabb) + idx2)

  tidx1 = setdiff(1:nrow(pos_feature_tabb), idx1)
  tidx2 = sample(setdiff(1:nrow(neg_feature_tabb), idx2), length(tidx1), replace=F)

  test_idx = c(tidx1, nrow(pos_feature_tabb) + tidx2)

  training_labels = labels[training_idx]
  test_labels = labels[test_idx]

  training_data = annotation_feature_tabb[training_idx, ]
  test_data = annotation_feature_tabb[test_idx, ]

  all_training[[num_iter]] = training_data

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
  predict_labels = predict(bstSparse, as.matrix(test_data))
  roc_test <- roc(test_labels, predict_labels, algorithm = 2)
  probs1 = predict_labels[test_labels == 1]
  probs2 = predict_labels[test_labels == 0]
  pr <- pr.curve(scores.class0 = probs1, scores.class1 = probs2, curve = T)
  aucvec = c(aucvec, auc(roc_test))
  prcvec = c(prcvec, pr$auc.integral)

  shap_values <- shap.values(xgb_model = bstSparse, X_train = as.matrix(training_data))
  all_shap[[num_iter]] <- shap_values$shap_score
  
  cat("We are at iter:", num_iter, "\n")
}

new_shap <- bind_rows(all_shap)
new_training <- bind_rows(all_training)
shap_long <- shap.prep(shap_contrib = new_shap, X_train = as.matrix(new_training))
write.csv(shap_long, paste0(sub('.[^.]*$', '', opt$output_metrics), ".shap_values.csv"), quote=FALSE, row.names = FALSE)
shap_long <- shap.prep(shap_contrib = new_shap, X_train = as.matrix(new_training), top_n = 25)
png(filename=paste0(sub('.[^.]*$', '', opt$output_metrics), ".png"), width = 225 , height = 175, units='mm', res = 900, type="cairo")
shap.plot.summary(shap_long, min_color_bound = "#008afb", max_color_bound = "#f70068",)
dev.off()

dff_roc = data.frame("AUROC.mean" = mean(aucvec),
                 "AUROC.sd" = sd(aucvec),
                 "AUPRC.mean" = mean(prcvec),
                 "AUPRC.sd" = sd(prcvec))
write.table(dff_roc, file = paste0(opt$output_metrics), col.names=TRUE, row.names=FALSE, sep = "\t", quote=F)
