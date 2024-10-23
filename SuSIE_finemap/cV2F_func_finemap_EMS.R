library(data.table)
library(R.utils)
library(optparse)
library(zoo)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--trait", type="character", default = "../data/biosamples.txt", help="Biosample names"),
  make_option("--finemappath", type="character", default = "../data/ukbb-finemapping",
	      help="Path to directory of finemapped traits"),
  make_option("--bimpath", type="character", default = "../data/1000G_BIMS_hg38/1000G.EUR.QC.",
	      help="Path and prefix of the bimfile of the BIMFILE"),
  make_option("--cv2fpath", type="character", default = "../results/score/mwe_baseline.",
	      help="Path and prefix of cv2f scores files")
  make_option("--output_finemap", type="character", default = "../results/SuSIE_finemap/EMS_all_cV2F",
	      help="Output path to directory of results")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

traitname = opt$trait
dff = read.delim2(paste0(opt$finemappath, "/", traitname, ".tsv.bgz"))

pooled_tabb_list = list()
for(numchr in 1:22){
  tabb = data.frame(fread(paste0(opt$bimpath, numchr, ".bim")))
  pooled_tabb_list[[numchr]] = tabb
}
pooled_tabb = do.call(rbind, pooled_tabb_list)


pre_cv2f_tabb = c()
for(numchr in 1:22){
  cv2f_scores = data.frame(fread(paste0(opt$cv2fpath, numchr, ".cv2f.txt")))
  pre_cv2f_tabb = rbind(pre_cv2f_tabb, cv2f_scores)
  cat("We are at chr:", numchr, "\n")
}
snps_intersect = intersect(pre_cv2f_tabb$SNP, pooled_tabb[,2])
pre_cv2f_tabb$BP_hg38 = pre_cv2f_tabb$BP
tt = rep(0.10, nrow(pre_cv2f_tabb))
tt[pre_cv2f_tabb$cV2F>0.75] = 1
pre_cv2f_tabb$cV2F = tt
cv2f_tabb = pre_cv2f_tabb[match(snps_intersect, pre_cv2f_tabb$SNP), ]
cv2f_tabb$BP = pooled_tabb[match(snps_intersect, pooled_tabb[,2]), 4]


u_regions = unique(dff$region)

out_merged = c()
for(numu in 1:length(u_regions)){
  dff_temp = dff[which(dff$region == u_regions[numu]), ]
  dff_temp2 = dff_temp[which(dff_temp$cs_id != -1), ]
  if(nrow(dff_temp2) == 0){
    cat("We are at region:", numu, "\n")
    next
  }
  xx = cbind.data.frame(as.numeric(dff_temp2$alpha1),
                        as.numeric(dff_temp2$alpha2),
                        as.numeric(dff_temp2$alpha3),
                        as.numeric(dff_temp2$alpha4),
                        as.numeric(dff_temp2$alpha5),
                        as.numeric(dff_temp2$alpha6),
                        as.numeric(dff_temp2$alpha7),
                        as.numeric(dff_temp2$alpha8),
                        as.numeric(dff_temp2$alpha9),
                        as.numeric(dff_temp2$alpha10))

  cv2f_temp = cv2f_tabb$cV2F[match(dff_temp2$rsid, cv2f_tabb$SNP)]
  Cz = zoo(cv2f_temp)
  Cz_approx2 <- na.approx(Cz, na.rm=FALSE, rule=2)
  Cz_approx = Cz_approx2
  # if(length(Cz_approx2) ==1){
  #   Cz_approx = (exp(Cz_approx2))/max(exp(Cz_approx2))
  # }else{
  #   Cz_approx = (exp(Cz_approx2)-min(exp(Cz_approx2)))/(max(exp(Cz_approx2)) - min(exp(Cz_approx2)))
  # }

  xx2 = cbind.data.frame(as.numeric(dff_temp2$alpha1)*Cz_approx,
                        as.numeric(dff_temp2$alpha2)*Cz_approx,
                        as.numeric(dff_temp2$alpha3)*Cz_approx,
                        as.numeric(dff_temp2$alpha4)*Cz_approx,
                        as.numeric(dff_temp2$alpha5)*Cz_approx,
                        as.numeric(dff_temp2$alpha6)*Cz_approx,
                        as.numeric(dff_temp2$alpha7)*Cz_approx,
                        as.numeric(dff_temp2$alpha8)*Cz_approx,
                        as.numeric(dff_temp2$alpha9)*Cz_approx,
                        as.numeric(dff_temp2$alpha10)*Cz_approx)

  ucs = unique(dff_temp2$cs_id)
  for(num_cs in 1:length(ucs)){

    idx= which(dff_temp2$cs_id == ucs[num_cs])
    xx_filtered = xx[idx, ]
    pip2 = apply(xx_filtered, 1, function(z) return(1 - prod(1-as.numeric(z))))

    effects_to_choose = as.numeric(which(colSums(xx_filtered) > 0.90))
    xx3 = xx2[idx, ]
    xx4 = sweep(xx3, 2, colSums(xx3), "/")
    xx5 = cbind.data.frame(xx4[,effects_to_choose])
    pip3 = apply(xx5, 1, function(z) return(1 - prod(1-as.numeric(z))))

    cv2f_out = Cz_approx[idx]
    dff_temp3 = cbind.data.frame(dff_temp2[idx, ], pip2, pip3, cv2f_out)

    numm = length(which(cumsum(sort(xx4[,effects_to_choose], decreasing = T)) < 0.95))+1
    denn = length(which(cumsum(sort(xx[,effects_to_choose], decreasing = T)) < 0.95))+1
    cs_shrink = (1 - numm/denn)
    dff_temp3$cs_shrink_percentage = cs_shrink

    colnames(dff_temp3) = c(colnames(dff_temp2), "pip", "cV2F.pip",
                            "cV2F", "CS.shrink.percent")
    out_merged = rbind(out_merged, dff_temp3)
  }
  cat("We are at region:", numu, "\n")
}

write.table(out_merged, file = paste0(output_finemap, "/", traitname, ".cv2f_ukbb.txt"),
            row.names=F, col.names=T, sep = "\t", quote=F)
