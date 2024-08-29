
#########################################
##### 1. Merge NeuroX1 & 2 Datasets #####
#########################################

###### 1.1 SNP QC before merge
###### 1.1.1 DATASET1
/share/apps/genomics/plink-1.9/plink --bfile NEUROX --missing --out dataset1_missing
# 267607 variants loaded from .bim file.
# Warning: 5465 het. haploid genotypes present (see dataset1_missing.hh ) // Note: het. haploid genotypes may be due to genotyping errors or due to the pseudo-autosomal regions on the X chromosome, where males have homologous sequences on both the X and Y chromosomes. This region is coded as 25 in PLINK

## TO CHECK WHETHER THESE SNPS ARE ON THE X CHROMOSOME
awk '{print $3}' dataset1_missing.hh | sort | uniq > hhSNPs.txt

/share/apps/genomics/plink-1.9/plink --bfile NEUROX --extract hhSNPs.txt --make-bed --out hhSNPs
# --extract: 183 variants remaining.

awk '{print $1}' hhSNPs.bim | sort | uniq -c
# None of these SNPs are on Chr25

## REMOVE THE HHSNPS
/share/apps/genomics/plink-1.9/plink --bfile NEUROX --exclude hhSNPs.txt --make-bed --out dataset1_rmhh
# 267424 variants and 619 people pass filters and QC.

###### 1.1.2 DATASET2
/share/apps/genomics/plink-1.9/plink --bfile BINARY_toLoni --missing --out dataset2_missing
# 990768 variants loaded from .bim file.
# Warning: 25272 het. haploid genotypes present (see dataset2_missing.hh )

awk '{print $3}' dataset2_missing.hh | sort | uniq > hhSNPs2.txt

/share/apps/genomics/plink-1.9/plink --bfile BINARY_toLoni --extract hhSNPs2.txt --make-bed --out hhSNPs2
# --extract: 662 variants remaining.

awk '{print $1}' hhSNPs2.bim | sort | uniq -c
# None of these SNPs are on Chr25

## REMOVE THE HHSNPS
/share/apps/genomics/plink-1.9/plink --bfile BINARY_toLoni --exclude hhSNPs2.txt --make-bed --out dataset2_rmhh
# 990106 variants and 222 people pass filters and QC.

###### 1.2 Subject QC before merge
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh --mind 0.05 --geno 0.05 --make-bed --out dataset1_rmhh_QC
# 13 people removed due to missing genotype data (--mind).
# 12440 variants removed due to missing genotype data (--geno).
# 254984 variants and 606 people pass filters and QC.

/share/apps/genomics/plink-1.9/plink --bfile dataset2_rmhh --mind 0.05 --geno 0.05 --make-bed --out dataset2_rmhh_QC
# 0 people removed due to missing genotype data (--mind).
# 35286 variants removed due to missing genotype data (--geno).
# 954820 variants and 222 people pass filters and QC.

###### 2. Merge datasets
###### 2.1 Convert SNP ID to RSID
# code for ID conversion please see RSID_convert.R, using dbSNP mapping
# output files: 1. *.map (SNP name & RS ID); 2. *_exclude (list of SNPs with no RS ID record)
awk '{print $7"\t"$2}' dataset1_rmhh_QC.bim_RSID.map  > dataset1_RSID.txt 
awk '{print $2}' dataset1_rmhh_QC.bim_noRS_exclude > dataset1_noRSID.txt

awk '{print $7"\t"$2}' dataset2_rmhh_QC.bim_RSID.map  > dataset2_RSID.txt 
awk '{print $2}' dataset2_rmhh_QC.bim_noRS_exclude > dataset2_noRSID.txt

## Exclude SNPs with no RSID converted
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC --exclude dataset1_noRSID.txt --make-bed --out dataset1_rmhh_QC_rm.noRS --keep-allele-order
# --exclude: 242050 variants remaining.
# 242050 variants and 606 people pass filters and QC.

/share/apps/genomics/plink-1.9/plink --bfile dataset2_rmhh_QC --exclude dataset2_noRSID.txt --make-bed --out dataset2_rmhh_QC_rm.noRS --keep-allele-order
# --exclude: 902177 variants remaining.
# 902177 variants and 222 people pass filters and QC.

## Update SNP names in dataset1 and dataset2 list
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS --update-name dataset1_RSID.txt 1 2 --make-bed --out dataset1_rmhh_QC_rm.noRS_RSconvert --keep-allele-order
# --update-name: 242050 values updated.

/share/apps/genomics/plink-1.9/plink --bfile dataset2_rmhh_QC_rm.noRS --update-name dataset2_RSID.txt 1 2 --make-bed --out dataset2_rmhh_QC_rm.noRS_RSconvert --keep-allele-order
# --update-name: 902177 values updated.

###### 2.2 Check for overlapping SNPs
awk '{print $2}' dataset1_rmhh_QC_rm.noRS_RSconvert.bim | sort > sorted_SNPset1.txt
awk '{print $2}' dataset2_rmhh_QC_rm.noRS_RSconvert.bim | sort > sorted_SNPset2.txt

comm -12 sorted_SNPset1.txt sorted_SNPset2.txt > comm_SNPset.txt
wc comm_SNPset.txt
#  223378  223378 2588561 comm_SNPset.txt

###### 2.3 Extract overlapping SNPs in dataset1 and 2
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert --extract comm_SNPset.txt --make-bed --out dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP
# 223378 variants and 606 people pass filters and QC.

/share/apps/genomics/plink-1.9/plink --bfile dataset2_rmhh_QC_rm.noRS_RSconvert --extract comm_SNPset.txt --make-bed --out dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP
# 223378 variants and 222 people pass filters and QC.

###### 2.4 Merge1
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP --bmerge  dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP --make-bed --out merge1
# Error: 43100 variants with 3+ alleles present. // NOTE: check these SNPs whether it is because of strand inconsistency

## There are some alleles have 0 or - on A1 // NOTE: not sure what happened but better to remove them
awk '$5 == "0" {print $2}' dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP.bim > snps_to_remove1.txt
wc snps_to_remove1.txt
# 124364  124364 1479873 snps_to_remove1.txt

awk '$5 == "-" {print $2}' dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP.bim > snps_to_remove2.txt
wc snps_to_remove2.txt
# 159817  159817 1894351 snps_to_remove2.txt

cat snps_to_remove1.txt snps_to_remove2.txt > combined_snps_to_remove.txt
sort combined_snps_to_remove.txt | uniq > unique_combined_list.txt
wc unique_combined_list.txt
# 167438  167438 1984681 unique_combined_list.txt

/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP --exclude unique_combined_list.txt --make-bed --out dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele
# 55940 variants and 606 people pass filters and QC.

/share/apps/genomics/plink-1.9/plink --bfile dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP --exclude unique_combined_list.txt --make-bed --out dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele
# 55940 variants and 222 people pass filters and QC.

###### 2.5 Merge2
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele --bmerge dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele --make-bed --out merge2
# Error: 21 variants with 3+ alleles present. // NOTE: these may be due to strand inconsistency

## FLIP
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele --flip merge2-merge.missnp --make-bed --out dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele_flip

/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele_flip --bmerge dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele --make-bed --out merge3
# Error: 1 variant with 3+ alleles present.
# this SNP has "I" on A1 //NOTE: exclude this SNP

/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele_flip --exclude merge3-merge.missnp --make-bed --out dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele_flip1
# 55939 variants and 606 people pass filters and QC.

/share/apps/genomics/plink-1.9/plink --bfile dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele --exclude merge3-merge.missnp --make-bed --out dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele1
# 55939 variants and 222 people pass filters and QC.

###### 2.6 Merge4
/share/apps/genomics/plink-1.9/plink --bfile dataset1_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele_flip1 --bmerge dataset2_rmhh_QC_rm.noRS_RSconvert_commSNP_rm3allele1 --make-bed --out merge4
# 55939 variants and 828 people pass filters and QC.


