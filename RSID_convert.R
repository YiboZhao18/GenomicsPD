# convert to RSID
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library (SNPlocs.Hsapiens.dbSNP144.GRCh37)

setwd("/home/yibozhao/mydata/IL12")

get_RSID <- function(file){
        mydata = read.table(file,stringsAsFactors=F,header=F)
        names(mydata) <- c("Chr", "Marker", "zero", "Bp", "A1", "A2")
        
        mydata$rsid = rep("norsid",nrow(mydata))
        mydata$key = paste0(mydata$Chr,":",mydata$Bp)
        chrs = unique(mydata$Chr)
        snps = SNPlocs.Hsapiens.dbSNP144.GRCh37 # SNP map for RSID
        snpcount(snps)
        cat("Well try to get RSID from",nrow(mydata),"variants\n")
        for(chr in chrs){
                mydata.local = mydata[mydata$Chr == chr,]
                cat("Converting to RSID at chromosome",chr,nrow(mydata.local),"snps\n")
                result = snpsByOverlaps(snps,GRanges(seqnames=rep(paste0(chr),nrow(mydata.local)), ranges=IRanges(start=mydata.local$Bp,end=mydata.local$Bp)))
                result = cbind(start(result),elementMetadata(result)$RefSNP_id)
                mask = !is.na(match(paste0(chr,":",result[,1]),mydata$key))
                mydata$rsid[match(paste0(chr,":",result[,1]),mydata$key)[mask]] = result[,2][mask]
                cat("We converted to RSID at chromosome",chr,sum(mask),"snps\n")
        }
        
        cat("We got converted",sum(mydata$rsid != "norsid"),"variants\n")
        RSID.map = mydata[mydata$rsid != "norsid",]
        norsid = mydata[mydata$rsid == "norsid",]
        outfile1 <- paste(file, "_RSID.map", sep = "")
        outfile2 <- paste(file, "_noRS_exclude", sep = "")
        
        write.table(RSID.map,outfile1,col.names=F,quote=F,row.names=F)
        write.table(norsid, outfile2, col.names = F, quote = F, row.names = F)
}

get_RSID("dataset1_rmhh_QC.bim")
get_RSID("dataset2_rmhh_QC.bim")