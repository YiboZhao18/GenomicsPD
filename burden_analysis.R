# burden analysis
# before analysis, prepare dir "burden" contains: .bed, .fam, .bim file from plink and b37_map.txt for gene-SNP mapping
# 1. Prepare SNP-gene map
gene.posi <- read.delim("b37_map.txt", header = F, stringsAsFactors = F)
gene.names <- gene.posi$V4
gene.start <- gene.posi$V2
gene.end <- gene.posi$V3

outfile <- "NeuroX_LRRK2.int_rm.noRS_RSconvert_coding_merge3_rm.palindromic_geno_hwe_Caucasian_pheno.onset_QC"
list.files()
bim.file <- paste(outfile, "bim", sep = ".")
SNP.bim <- read.delim(bim.file, stringsAsFactors = F, header = F)
SNP.id <- SNP.bim$V2
SNP.start <- SNP.bim$V4

i = 1
total <- length(SNP.id)
SetID <- integer()
for (i in 1:total){
        j = 1
        all.genes <- length(gene.names)
        for (j in 1:all.genes){
                if ((SNP.start[i] > gene.start[j]) & (SNP.start[i] < gene.end[j])) {
                        flag <- "match found!"
                        SetID[i] <- gene.names[j]
                }
                
        }
}
SetID.file <- paste(outfile, "SetID", sep = ".")
write.table(cbind(SetID, SNP.id), SetID.file, sep = "\t", row.names = F, quote = F)

# Generate a SNP set data file (SSD) from binary plink data files using user specified SNP sets
library(SKAT)

all.files <- list.files()
all.files

outSSD <- paste(outfile, "SSD", sep = ".")
outInfo <- paste(outfile, "Info", sep = ".")
# (Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
Generate_SSD_SetID(all.files[6], all.files[7], all.files[8], all.files[11], outSSD, outInfo)

# Generate a genotype matrix (from plink)
SSD <- Open_SSD(outSSD, outInfo)
# get genotype matrix from plink
# /share/apps/genomics/plink-1.9/plink --bfile IMMUNO_rm.unlifted_liftOver.b37_sorted_LRRK2int_rm.noRS_RSconvert_coding_rm.palindromic_geno_hwe_pheno1_caucasian --recode A-transpose --out genotype_matrix
genotype_matrix <- read.delim("NeuroX_merge_genotype_matrix.traw", stringsAsFactors = F)
# for SKAT, a genotype matrix has samples as rows while SNPs as columns, no headers
row.names(genotype_matrix) <- genotype_matrix$SNP
Z.t <- genotype_matrix[,-c(1:6)]
Z <- t(Z.t)

# Generate a parameter file (with covariates)
### 1. Prepare a data file like SKAT.example
###### Z: genotype matrix; X: covaritate; y: phenotype (0/1)
fam.file <- paste(outfile, "fam", sep = ".")
in.fam <- read.table(fam.file, stringsAsFactors = F, header = F)
y <- in.fam[,6]
i = 1
total.y <- length(y)
for (i in 1:total.y){
        if (y[i] == 2){
                y[i] = 0
        }
}
# get covariate matrix generated from PCA on single SNP analyses
in.covariate <- read.table("sex_PC1.5", header = T, stringsAsFactors = F)
X <- as.matrix(in.covariate[,-c(1,2)])
obj <- SKAT_Null_Model(y~X, out_type = "D")

i = 1
total <- length(SSD$SetInfo$SetIndex)
pval <- integer()
for (i in 1:total) {
        res <- SKAT.SSD.OneSet_SetIndex(SSD, i, obj, kernel = "linear.weighted")
        pval[i] <- res$p.value
}
adj.p <- p.adjust(pval, "bonferroni", total)
output <- cbind(SSD$SetInfo$SetID, pval, adj.p)
output
write.table(output, "burden_res_NeuroX_merge.txt", row.names = F, quote = F, sep = "\t")
# summary(pval<0.05)
# SSD$SetInfo$SetID[pval<0.05]

