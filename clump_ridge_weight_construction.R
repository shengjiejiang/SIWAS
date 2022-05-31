slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job <- as.numeric(slurm_arrayid)
chr <- job

library(data.table)
library(matlib)
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(MASS)

suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
source("dist_support.R")
source("ACAT.R")

# This function (allele.qc) is downloaded from TWAS (http://gusevlab.org/projects/fusion/#typical-analysis-and-output)

allele.qc = function(a1,a2,ref1,ref2) {
    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip
    
    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;
    
    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    return(snp)
}

JointRidge = function(B,S,N,XX=diag(1,nrow=1),lambda= 0.1){ #x first
    
    B = as.matrix(B)
    S = as.matrix(S)
    N = as.matrix(N)
    n=max(N)
    nrx = dim(B)[1]
    yy1 = NULL
    
    for(rx in 1:nrx){
        Dj = XX[rx,rx]
        Sj = S[rx]
        Bj = B[rx]
        yy1 = c(yy1, N[rx,1] / n * (n-1) * Dj * Sj^2 + Dj * Bj^2)
    }
    yy = median(yy1)
    xxd = diag(XX)
    XY1 = xxd * B      #vector
    B = c(XY1)
    sA = ginv( (XX + diag(lambda,dim(XX)[1],dim(XX)[2])) )
    beta = sA %*% B
    eign = eigen(XX)$values
    df = sum(eign/(eign + lambda))
    sigma2 = (yy - t(beta) %*% B ) / (n- df )
    se = sigma2[1,1] * sA %*% XX %*% sA   #cov
    pvalue = 2 * pnorm(abs(beta/sqrt(diag(se))) ,lower.tail=FALSE)
    return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2,yy = yy,xy = XY1,df = df))
}


ref_ld <- "/gpfs/home/sj14m/summerproject/project2/1000G/1000G.EUR.QC.CHR"
genos <- read_plink(paste(ref_ld,job,sep=""), impute="avg")
VOL.file <- as.data.frame(read.table(paste("/gpfs/home/sj14m/summerproject/project3/VOL_chr", job, ".txt", sep=""), header=TRUE))
dim(VOL.file)
snp.used <- match(VOL.file$RSID, genos$bim[,2])
keep <- !is.na(snp.used)
genos.bim <- genos$bim[snp.used[keep],]
genos.bed <- genos$bed[,snp.used[keep]]
dim(genos.bim)
dim(genos.bed)

# QC / allele-flip the input and output
qc <- allele.qc( VOL.file$A1 , VOL.file$A2 , genos.bim[,5] , genos.bim[,6] )

# Flip Z-scores for mismatching alleles
VOL.file[ c(which(qc$flip=="TRUE")), c(paste("BETA", 1:58, sep="")) ] = -1 * VOL.file[ c(which(qc$flip=="TRUE")), c(paste("BETA", 1:58, sep="")) ]
VOL.file$A1[ qc$flip ] = genos.bim[qc$flip,5]
VOL.file$A2[ qc$flip ] = genos.bim[qc$flip,6]
VOL.file$POS[qc$keep] = genos.bim[qc$keep, 4]

# Remove strand ambiguous SNPs (if any)
# which (qc$keep=="FLASE")  #integer(0)
if ( sum(!qc$keep) > 0 ) {
    genos$bim = genos$bim[qc$keep,]
    genos$bed = genos$bed[,qc$keep]
    VOL.file = VOL.file[qc$keep,]
}

gene2000chr <- as.data.frame(read.table(paste("gene2000chr", job, ".txt", sep=""), header=TRUE))
head(gene2000chr)

n <- nrow(VOL.file)
d <- 1
VOL.file_new <- VOL.file[1,]
for (i in 1:n){
	print(i)
    	A <- (VOL.file[i,2] >= gene2000chr$start) & (VOL.file[i,2] <= gene2000chr$end)
    	k <-which(A)
    	B <- all(!A)
    	m <- sum(A)
    	for(l in 1:m){
        	VOL.file_new[d, ] <- VOL.file[i,]
        	VOL.file_new$GeneID[d] <- ifelse(B, "NULL", as.character(gene2000chr[k[l],2]))
        	d <- d+1
    	}
}

filename <- paste("/gpfs/home/sj14m/summerproject/project3/simulation/VOL.file_new_",job,".rds", sep="")
saveRDS(VOL.file_new, filename)

VOL.file_new <- readRDS(filename)

VOL.file_new <- VOL.file_new[order(VOL.file_new$GeneID, VOL.file_new$POS), ]
head(VOL.file_new)
dim(VOL.file_new)  #134514*123


VOL.file_new2 <- VOL.file_new[VOL.file_new$GeneID != "NULL",]
z <- split(VOL.file_new2, VOL.file_new2$GeneID)
length(z)

dim(VOL.file_new2) #52798 * 123
names(VOL.file_new2)[123] <- "Gene.ENSG"
names(VOL.file_new2)[2]   <- "BP"
names(VOL.file_new2)[3]   <- "SNP"
names(VOL.file_new2)[1]   <- "chr"


y <- split(VOL.file_new2, VOL.file_new2$Gene.ENSG)
m <- length(y)

wl <- as.data.frame(matrix(NA,0,3))
colnames(wl) <- c("GeneID", "CHR", "file.path")

l<-1

for (i in 1:2){  #m
	print(i)
	k <- z[[i]][,1:6]
    for (j in 1:2){ #58
		ss <- y[[i]][, c("SNP", "Gene.ENSG", "chr", "BP", paste("BETA", j, sep=""), paste("SEBETA", j, sep="") )]
		ss$Zscore <- ss[,5]/ss[,6]
		
		tmp.outdir = "/gpfs/home/sj14m/summerproject/project3/simulation/clumps/"
		tmp.snpfile = paste(tmp.outdir, ss[1, "Gene.ENSG"],".txt", sep="")
		write.table(ss[,1], file =tmp.snpfile, quote = F, col.names =F, row.names = F)
		system(paste0("plink --bfile ", "1000G.EUR.ALLSNP.QC.CHR", chr, " --chr ", chr, "  --extract ", tmp.snpfile, " --make-bed --out ", tmp.outdir, ss[1, "Gene.ENSG"]))
		ss$pval = 2*pnorm(abs(ss$Zscore), lower.tail=F)
		tmp.ssfile = paste(tmp.outdir, ss[1, "Gene.ENSG"], "_ss", ".txt", sep="")
		write.table(ss, file=tmp.ssfile, quote =F, col.names = T, row.names = F)
		command <- paste("plink --bfile ", tmp.outdir, ss[1, "Gene.ENSG"], " --clump ", tmp.ssfile, " --clump-p1 1 --clump-kb 250 --clump-r2 0.95  --clump-snp-field SNP --clump-field pval --out ", tmp.outdir, ss[1, "Gene.ENSG"], sep = "")
		system(command)
		
		snps.keep <- read.table(paste(ss[1, "Gene.ENSG"], ".clumped", sep=""), header=TRUE)
		nr <- length(snps.keep$SNP)
		if (nr<=1){
			next
		}
		N <- rep(8428, nr)
		used.snpx <- snps.keep$SNP
		cur.genos <- genos.bed[, used.snpx]
		cur.genos <- cur.genos-matrix(rep(colMeans(cur.genos), each=dim(cur.genos)[1]), dim(cur.genos)[1], dim(cur.genos)[2])
		XX <- cov(cur.genos)
		B <- y[[i]][y[[i]]$SNP %in% used.snpx, c(paste("BETA", j, sep=""))]
		S <- y[[i]][y[[i]]$SNP %in% used.snpx, c(paste("SEBETA", j, sep=""))]
		A <- JointRidge(B,S,N,XX,lambda= 0.1)
		z[[i]][, c(paste("BETA", j, sep=""))] <- 0
		z[[i]][y[[i]]$SNP %in% used.snpx, c(paste("BETA", j, sep=""))] <- A$beta
		z[[i]][, c(paste("SEBETA", j, sep=""))] <- 0
		z[[i]][y[[i]]$SNP %in% used.snpx, c(paste("SEBETA", j, sep=""))] <- sqrt(diag(A$cov))
		
		a <- unique(z[[i]]$GeneID)
    		file.out<- paste("/gpfs/home/sj14m/summerproject/project3/simulation/clumps/weights/VOL_gene2000chr",job, "_", a, ".wgt.rds", sep="")
    		saveRDS(z[[i]], file.out)
    		wl[l,1] <- as.character(a)
    		wl[l,2] <- job
    		wl[l,3] <- paste("/gpfs/home/sj14m/summerproject/project3/simulation/clumps/weights/VOL_gene2000chr",job, "_", a, ".wgt.rds", sep="")
    		l <- l+1
		print(l)
	}
}
dim(wl)
head(wl)
saveRDS(wl, paste("/gpfs/home/sj14m/summerproject/project3/simulation/clumps/weights/VOL_gene2000_chr",job, "_weight.info.rds", sep=""))





