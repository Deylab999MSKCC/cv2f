library(data.table)
library(locuszoomr)
library("EnsDb.Hsapiens.v86")
library(R.utils)
library(optparse)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--tissue", type="character", default = "ALL",
	      help="Whether we are using tissue/tissue-agnostic finemapped results. Currently only supports ALL, BLOOD, and LIVER"),
  make_option("--trait", type="character", default = "FEV1FVC",
	      help="The complex trait from summary statistic"),
  make_option("--rsid", type="character", default = "rs17168118",
	      help="The rsid of the variant to center"),
  make_option("--outputcell", type="character", default = "../results/plots",
	      help="Output folder for the plot"),
  make_option("--cv2f_dir", type="character", default = "../results/SuSIE_finemap/EMS_all_cV2F",
	      help="Directory of variants finemapped with cV2F scores"),
  make_option("--trait_cv2f_dir", type="character", default = "../results/SuSIE_finemap/EMS_all_Liver_cV2F",
              help="Directory of variants finemapped with trait specific cV2F score. Only required if not tissue-agnostic"),
  make_option("--sumstats_dir", type="character", default = "../data/sumstats/sumstats_Dey",
	      help="Directory of GWAS summary stats"),
  make_option("--bimpath", type="character", default = "../data/1000G_BIMS_hg38/1000G.EUR.QC.",
	      help="Path and prefix of the bimfile of the BIMFILE")
)


opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

tissue = opt$tissue
trait = opt$trait
rsid = opt$rsid
outputcell = opt$outputcell
cv2f_dir = opt$cv2f_dir
trait_cv2f_dir = opt$trait_cv2f_dir
sumstats_dir = opt$sumstats_dir
bimpath = opt$bimpath


cv2f = read.table(paste0(cv2f_dir, "/UKBB.", trait, ".SuSiE.cv2f_ukbb.txt"), sep="\t", header=TRUE)
cv2f = cv2f[which(cv2f$chromosome!="X" & cv2f$chromosome!="Y"),]

if (tissue == "ALL") {
    data = cv2f[which(cv2f$rsid == rsid),]
    chrm = data$chromosome
    pos = data$position
    left_lim = max(0, pos - 100000)
    right_lim = pos + 100000
} else {
    trait_cv2f = read.table(paste0(trait_cv2f_dir, "/UKBB.", trait, ".SuSiE.cv2f_ukbb.txt"), sep="\t", header=TRUE)
    trait_cv2f = trait_cv2f[which(trait_cv2f$chromosome!="X" & trait_cv2f$chromosome!="Y"),]

    data = trait_cv2f[which(trait_cv2f$rsid == rsid),]
    chrm = data$chromosome
    pos = data$position
    left_lim = max(0, pos - 100000)
    right_lim = pos + 100000
}

bim = read.table(paste0(bimpath, chrm, ".bim"), sep="\t")
if (file.exists(paste0(sumstats_dir, "/UKB.", trait, ".BOLT.sumstats.gz"))) {
    sumstats = read.table(paste0(sumstats_dir, "/UKB.", trait, ".BOLT.sumstats.gz"), sep="\t", header=TRUE)
} else {
    sumstats = read.table(paste0(sumstats_dir, "/UKB.", trait, ".SAIGE.sumstats.gz"), sep="\t", header=TRUE)
}
colnames(bim) <- c("chromosome", "rsid", ".", "position", "ref", "alt")
#colnames(sumstats) <- c("rsid", "A1", "A2", "Z", "N")
names(sumstats)[names(sumstats) == 'SNP'] <- 'rsid'

bim = merge(bim, sumstats, sort=FALSE)
bim$P = pnorm(abs(bim$Z), 0, 1, lower.tail=FALSE)
#bim = bim[which(bim$P>10**(-20)),]

general_data = bim[which(bim$rsid == rsid),]
general_chrm = general_data$chromosome
general_pos = general_data$position
general_left_lim = max(0, general_pos - 100000)
general_right_lim = general_pos + 100000

pdf(paste0(outputcell,"/",tissue,"_",trait,"_",rsid,"_plot.pdf"), width = 6, height = 6, pointsize=8)

if (tissue != "ALL") {
    print("here")
    loc <- locus(data = bim, seqname=general_chrm, xrange=c(general_left_lim, general_right_lim), ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid", p="P")
    loc_pip <- locus(data = trait_cv2f, seqname=chrm, xrange=c(left_lim, right_lim), yvar = "pip",
                        ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid")
    loc_cv2f <- locus(data = trait_cv2f, seqname=chrm, xrange=c(left_lim, right_lim), yvar = "cV2F", 
                        ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid")
    loc_cv2fpip <- locus(data = trait_cv2f, seqname=chrm, xrange=c(left_lim, right_lim), yvar = "cV2F.pip",
                        ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid")
    loc_cv2fpip_general <- locus(data = cv2f, seqname=chrm, xrange=c(left_lim, right_lim), yvar = "cV2F.pip",
                        ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid")
    pf_bim <- quote({
        abline(v = general_pos, lty=2, col = "blue")
    })
    pf <- quote({
        abline(h = 0.9, lty=2, col = "orange")
        abline(v = pos, lty=2, col = "blue")
    })
    
    if (tissue == "LIVER") {
        pl <- quote({
            legend("bottomright", legend = c("LIVER cV2F.pip", "cV2F.pip", "pip"),
                 pch = c(23, 23, 21), cex=1.8, pt.bg = c("green", "blue", "gray"), bty = "n")
        })
        oldpar <- set_layers(3)
        scatter_plot(loc, xticks = FALSE, cex=1, lwd=0.5, bg="green", legend_pos = 'topleft', panel.first=pf_bim, align = FALSE)
        
        scatter_plot(loc_cv2fpip_general, xticks = FALSE, ylim = c(0,1), cex=2.5, pch=23, bg="blue", lwd=0.5, cex.axis=1.3)
        scatter_plot(loc_cv2fpip, xticks = FALSE, ylim = c(0,1), cex=2.5, pch=23, bg="green", lwd=0.5, add = TRUE)
        scatter_plot(loc_pip, xticks = FALSE, ylim = c(0,1), cex=2.2, pch=21, bg="gray", lwd=0.5, panel.first = pf, panel.last = pl, add = TRUE)
        
    } else {
        pl <- quote({
            legend("bottomright", legend = c("BLOOD cV2F.pip", "cV2F.pip", "pip"),
                 pch = c(23, 23, 21), cex=1.8, pt.bg = c("red", "blue", "gray"), bty = "n")
        })
        oldpar <- set_layers(3)
        scatter_plot(loc, xticks = FALSE, cex=1, lwd=0.5, bg="green", legend_pos = 'topleft', panel.first=pf_bim, align = FALSE)
        
        scatter_plot(loc_cv2fpip_general, xticks = FALSE, ylim = c(0,1), cex=2.5, pch=23, bg="blue", lwd=0.5, cex.axis=1.3)
        scatter_plot(loc_cv2fpip, xticks = FALSE, ylim = c(0,1), cex=2.5, pch=23, bg="red", lwd=0.5, add = TRUE)
        scatter_plot(loc_pip, xticks = FALSE, ylim = c(0,1), cex=2.2, pch=21, bg="gray", lwd=0.5, panel.first = pf, panel.last = pl, add = TRUE)
    }

    genetracks(loc, filter_gene_biotype = "protein_coding")
    
} else {
    loc <- locus(data = bim, seqname=general_chrm, xrange=c(general_left_lim, general_right_lim), ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid", p="P")
    loc_pip <- locus(data = cv2f, seqname=chrm, xrange=c(left_lim, right_lim), yvar = "pip",
                        ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid")
    loc_cv2f <- locus(data = cv2f, seqname=chrm, xrange=c(left_lim, right_lim), yvar = "cV2F", 
                        ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid")
    loc_cv2fpip_general <- locus(data = cv2f, seqname=chrm, xrange=c(left_lim, right_lim), yvar = "cV2F.pip",
                        ens_db = "EnsDb.Hsapiens.v86", chrom="chromosome", pos="position", labs="rsid")

    pf_bim <- quote({
        abline(v = general_pos, lty=2, col = "blue")
    })
    pf <- quote({
        abline(h = 0.9, lty=2, col = "orange")
        abline(v = pos, lty=2, col = "blue")
    })
    pl <- quote({
        legend("bottomright", legend = c("cV2F.pip", "pip"),
             pch = c(23, 21), cex=1.8, pt.bg = c("blue", "gray"), bty = "n")
    })
    oldpar <- set_layers(3)
    scatter_plot(loc, xticks = FALSE, cex=1, lwd=0.5, bg="green", legend_pos = 'topleft', panel.first=pf_bim, align = FALSE)
    
    scatter_plot(loc_cv2fpip_general, xticks = FALSE, ylim = c(0,1), cex=2.5, pch=23, bg="blue", lwd=0.5, cex.axis=1.3)
    scatter_plot(loc_pip, xticks = FALSE, ylim = c(0,1), cex=2.2, pch=21, bg="gray", lwd=0.5, panel.first = pf, panel.last = pl, add = TRUE)
    genetracks(loc, filter_gene_biotype = "protein_coding")
}

par(oldpar)  # revert par() settings
dev.off()

