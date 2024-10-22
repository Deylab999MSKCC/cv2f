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
  make_option("--positive_set_train", type="character", default = "../data/positive_set.txt",
              help="rsIDs of SNPs in the positive set for training (celltype relevant trait)"),
  make_option("--negative_set_train", type="character", default = "../data/negative_set.txt",
              help="rsIDs of SNPs in the negative set for training (non-celltype relevant traits)"),
  make_option("--positive_set_test", type="character", default = "../data/positive_set.txt",
              help="rsIDs of SNPs in the positive set for validation (celltype relevant trait)"),
  make_option("--negative_set_test", type="character", default = "../data/negative_set.txt",
              help="rsIDs of SNPs in the negative set for validation (non-celltype relevant traits)"),
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

pre_positive_train_snps = read.table(opt$positive_set_train, header=F)[,1]
pre_negative_train_snps = read.table(opt$negative_set_train, header=F)[,1]
pre_positive_test_snps = read.table(opt$positive_set_test, header=F)[,1]
pre_negative_test_snps = read.table(opt$negative_set_test, header=F)[,1]
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

positive_train_snps = intersect(pre_positive_train_snps, maf_tabb$SNP)
positive_test_snps = intersect(pre_positive_test_snps, maf_tabb$SNP)

pre_negative_train_snps = setdiff(intersect(pre_negative_train_snps, maf_tabb$SNP), pre_positive_train_snps)
pre_negative_test_snps = setdiff(intersect(pre_negative_test_snps, maf_tabb$SNP), pre_positive_test_snps)

negative_train_snps = c()
negative_test_snps = c()
for(numl in 1:nrow(LDblocks_tabb)){
  snps_in_block = bimpooltabb$SNP[which(bimpooltabb$CHR == LDblocks_tabb$CHR[numl] &
                                  bimpooltabb$BP > LDblocks_tabb$START[numl] &
                                  bimpooltabb$BP < LDblocks_tabb$STOP[numl]) ]
  positives_train_in_block = intersect(positive_train_snps, snps_in_block)
  negatives_train_in_block = intersect(snps_in_block, pre_negative_train_snps)
  positives_test_in_block = intersect(positive_test_snps, snps_in_block)
  negatives_test_in_block = intersect(snps_in_block, pre_negative_test_snps)

  if(length(positives_train_in_block) > 0){
    mafs_positive_train_block = maf_tabb$MAF[match(positives_train_in_block, maf_tabb$SNP)]
    mafs_negative_train_block = maf_tabb$MAF[match(negatives_train_in_block, maf_tabb$SNP)]
    snps_temp = c()
    for(tt in 1:length(mafs_positive_train_block)){
      num_snps = pmin(5, length(mafs_negative_train_block))
      snps_temp = c(snps_temp, negatives_train_in_block[order(abs(mafs_positive_train_block[tt] - mafs_negative_train_block), decreasing = F)[1:num_snps]])
    }
    negative_train_snps = c(negative_train_snps, snps_temp)
  }
    
  if(length(positives_test_in_block) > 0){
    mafs_positive_test_block = maf_tabb$MAF[match(positives_test_in_block, maf_tabb$SNP)]
    mafs_negative_test_block = maf_tabb$MAF[match(negatives_test_in_block, maf_tabb$SNP)]
    snps_temp = c()
    for(tt in 1:length(mafs_positive_test_block)){
      num_snps = pmin(5, length(mafs_negative_test_block))
      snps_temp = c(snps_temp, negatives_test_in_block[order(abs(mafs_positive_test_block[tt] - mafs_negative_test_block), decreasing = F)[1:num_snps]])
    }
    negative_test_snps = c(negative_test_snps, snps_temp)
  }
  cat("Matching for LD-block:", numl, "\n")
}

negative_train_snps = negative_train_snps[!is.na(negative_train_snps)]
negative_test_snps = negative_test_snps[!is.na(negative_test_snps)]


positive_train_snps2 = c()
negative_train_snps2 = c()
pos_feature_train_tabb_list = c()
neg_feature_train_tabb_list = c()

positive_test_snps2 = c()
negative_test_snps2 = c()
pos_feature_test_tabb_list = c()
neg_feature_test_tabb_list = c()

temp_dff = c()
for(numchr in 1:22){
  temp_dff[[numchr]] = data.frame(fread(paste0(opt$feature_file, numchr, ".txt")))
}

##########################################################################################################

#combined_feature_tabb = rbind(pos_feature_tabb, neg_feature_tabb)
#rownames(combined_feature_tabb) = c(positive_snps2, negative_snps2)
#labels = c(rep(1, nrow(pos_feature_tabb)), rep(0, nrow(neg_feature_tabb)))
#annotation_feature_tabb = combined_feature_tabb


cat("Compute performance metrics to assess classification accuracy")
cat("\n")

all_training = list()

aucvec = c()
prcvec = c()
for(numchr in 1:22){
    positive_train_snps2 = c()
    negative_train_snps2 = c()
    pos_feature_train_tabb_list = c()
    neg_feature_train_tabb_list = c()

    positive_test_snps2 = c()
    negative_test_snps2 = c()
    pos_feature_test_tabb_list = c()
    neg_feature_test_tabb_list = c()

  for (other_numchr in 1:22){
      if (numchr != other_numchr) {
          dff = temp_dff[[other_numchr]]
          pos_feature_train_tabb_list[[other_numchr]] = dff[match(intersect(positive_train_snps, dff$SNP), dff$SNP), -(1:4)]
          neg_feature_train_tabb_list[[other_numchr]] = dff[match(intersect(negative_train_snps, dff$SNP), dff$SNP), -(1:4)]
          positive_train_snps2 = c(positive_train_snps2, intersect(positive_train_snps, dff$SNP))
          negative_train_snps2 = c(negative_train_snps2, intersect(negative_train_snps, dff$SNP))
          cat("We are at chr:", other_numchr, "\n")
          
      }
  }

    cat("now appending annotations\n")
  pos_feature_train_tabb = do.call(rbind, pos_feature_train_tabb_list)
  neg_feature_train_tabb = do.call(rbind, neg_feature_train_tabb_list)
  rownames(pos_feature_train_tabb) = positive_train_snps2
  rownames(neg_feature_train_tabb) = negative_train_snps2

  combined_feature_train_tabb = rbind(pos_feature_train_tabb, neg_feature_train_tabb)
  rownames(combined_feature_train_tabb) = c(positive_train_snps2, negative_train_snps2)
  labels_train = c(rep(1, nrow(pos_feature_train_tabb)), rep(0, nrow(neg_feature_train_tabb)))
  annotation_feature_train_tabb = combined_feature_train_tabb

    cat("now creating test set\n")
  dff = temp_dff[[numchr]]
  pos_feature_test_tabb_list[[numchr]] = dff[match(intersect(positive_test_snps, dff$SNP), dff$SNP), -(1:4)]
  neg_feature_test_tabb_list[[numchr]] = dff[match(intersect(negative_test_snps, dff$SNP), dff$SNP), -(1:4)]
  positive_test_snps2 = c(positive_test_snps2, intersect(positive_test_snps, dff$SNP))
  negative_test_snps2 = c(negative_test_snps2, intersect(negative_test_snps, dff$SNP))
  cat("We are at chr:", numchr, "\n")

  pos_feature_test_tabb = do.call(rbind, pos_feature_test_tabb_list)
  neg_feature_test_tabb = do.call(rbind, neg_feature_test_tabb_list)
  rownames(pos_feature_test_tabb) = positive_test_snps2
  rownames(neg_feature_test_tabb) = negative_test_snps2

    cat("is intersection\n")
  combined_feature_test_tabb = rbind(pos_feature_test_tabb, neg_feature_test_tabb)
  rownames(combined_feature_test_tabb) = c(positive_test_snps2, negative_test_snps2)
  labels_test = c(rep(1, nrow(pos_feature_test_tabb)), rep(0, nrow(neg_feature_test_tabb)))
  annotation_feature_test_tabb = combined_feature_test_tabb
    
    cat("separating out the stuff\n")
    if (nrow(pos_feature_train_tabb) < nrow(neg_feature_train_tabb)) {
      idx1 = seq(1, nrow(pos_feature_train_tabb))
      idx2 = sample(1:nrow(neg_feature_train_tabb), length(idx1), replace =	F)
    } else {
        idx2 = seq(1, nrow(neg_feature_train_tabb))
        idx1 = sample(1:nrow(pos_feature_train_tabb), length(idx2), replace =	F) 
    }
  training_idx = c(idx1, nrow(pos_feature_train_tabb) + idx2)
    
    if (nrow(pos_feature_test_tabb) < nrow(neg_feature_test_tabb)) {
      tidx1 = seq(1, nrow(pos_feature_test_tabb))
      tidx2 = sample(1:nrow(neg_feature_test_tabb), length(tidx1), replace=F)
    } else {
      tidx2 = seq(1, nrow(neg_feature_test_tabb))
      tidx1 = sample(1:nrow(pos_feature_test_tabb), length(tidx2), replace=F) 
    }
  test_idx = c(tidx1, nrow(pos_feature_test_tabb) + tidx2)

  training_labels = labels_train[training_idx]
  test_labels = labels_test[test_idx]

  training_data = annotation_feature_train_tabb[training_idx, ]
  test_data = annotation_feature_test_tabb[test_idx, ]

  all_training[[numchr]] = training_data

    cat("training")
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
  
  cat("We are at iter:", numchr, "\n")
}

dff_roc = data.frame("AUROC.mean" = mean(aucvec),
                 "AUROC.sd" = sd(aucvec),
                 "AUPRC.mean" = mean(prcvec),
                 "AUPRC.sd" = sd(prcvec))
write.table(dff_roc, file = paste0(opt$output_metrics), col.names=TRUE, row.names=FALSE, sep = "\t", quote=F)
