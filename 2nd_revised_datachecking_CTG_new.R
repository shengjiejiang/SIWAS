
library(data.table)
library(matlib)
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(MASS)

suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))

#suppressMessages(library('aSPU2'))
source("dist_support.R")
source("ACAT.R")

col <- 60
row <- 5
ind <- rep(0, 3*row*col)
content <- rep(0, row*col)
a1 <- c()
a2 <- c()
a3 <- c()
a4 <- c()
a5 <- c()
d1 <- c()
d2 <- c()
d3 <- c()
d4 <- c()
d5 <- c()
e1 <- c()
e2 <- c()
e3 <- c()
e4 <- c()
e5 <- c()
out <- matrix(NA, 0, 11)
list1 <- matrix(NA, 0, 12)
list2 <- matrix(NA, 0, 17)
list3 <- matrix(NA, 0, 12)
list4 <- matrix(NA, 0, 17)
list5 <- matrix(NA, 0, 12)
list6 <- matrix(NA, 0, 17)
list7 <- matrix(NA, 0, 12)
list8 <- matrix(NA, 0, 17)
list9 <- matrix(NA, 0, 12)
list10<- matrix(NA, 0, 17)
DD <- matrix(NA, 0, 0)

for (i in 1:22){
    print(i)    # for chromosome 1
    file.name1 = paste("/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/2nd_revised_VOL_gene2000_output_new",i, "_1.txt", sep="")
    file.name2 = paste("/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/2nd_revised_VOL_gene2000_summary_pVal_", i, "_1.rds", sep="")
    file.name3 = paste("/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/2nd_revised_VOL_gene2000_summary_pVal_", i, "_1.rds", sep="")
    #file.name2 <- paste("summary_best_cis_CTG.rds", sep="")
    #file.name3 <- paste("summary_best_cis_AD_Jansene.rds", sep="")    
	tmp <- read.table("/gpfs/home/sj14m/summerproject/project3/ENSEMBL_GRch37_gene_list.txt", header=TRUE)
	dim(tmp)   #24,811*7
	head(tmp)  #1.ensembl_gene_id 2.gene_biotype 3. hgnc_symbol 4. chromosome 5. start_position 6. end_position 7. uniqind
    A <- read.table(file.name1, header=TRUE) #125,198*9
	AA <- A[which(!is.na(A$gene)),]
	dim(AA)   # 136,290*9
	head(AA)  # 1. gene 2. CHR 3. nonzero_SNPs 4. weight_mod 5. TWAS_asy 6. SSU_asy 7. minP 8. ACAT1 9. ACAT2
	
m <- match(AA$gene, tmp$ensembl_gene_id)
	AA$GeneName <- tmp[m, c("hgnc_symbol")]
	dim(AA)   #  136,290*10
	head(AA)   # 1. gene 2. CHR 3. CHR 4. weight_mod 5. TWAS_asy 6. SSU_asy 7. minP 8. ACAT1 9. ACAT2 10. GeneName
	
	B <- readRDS(file.name2)
	dim(B)   #2,567*5
	head(B)  #1. CHR 2. P0 3. P1 4. GeneName 5. GeneID 6. pVal
	
	C <- readRDS(file.name3)
	dim(C)   # 2567*5
	head(C)  #1. CHR 2. P0 3. P1, 4. GeneName 5. GeneID  6. pVal
    
    for ( l in 1:nrow(AA)){   #nrow(A)=123,769
print(l)
       if (length(which(B$GeneName==AA[l,10]))==0){
           AA[l,11] <- "NA"
       } else {
           AA[l,11] <- min(B[which(B$GeneName == AA[l,10]), 6])  # put the smallest p value to 11th line of A
       }
       if (length(which(C$GeneName==AA[l,10]))==0){
           AA[l,12] <- "NA"
       } else {
           AA[l,12] <- min(unique(C[which(C$GeneName ==AA[l,10]), 6]))   # put the smallest p value to 12th line of A
       }     
    
   }
    names(AA)[11]<- c("pVal1")
    names(AA)[12]<- c("pVal2")
    dim(AA)   #123,769*12

    
    E <- AA[AA$weight_mod == 59, ]    # write Unweighted results in E   2,098*12
    D <- AA[AA$weight_mod < 59, ]    # write 1-58 results in D  121,671*12
    a <- unique(D$gene)                  # unique name of 2098 genes
    b <- length(a)          # 2098 number of valid genes
    for (d in 1:b){                      # calculate the cauchy combination results of 1-58
        if ( nrow(D[D$gene==a[d],])!= 58){
	    D <- D[-(D$gene ==a[d]),]
            next
        }else{
        D[D$gene==a[d], 13]<- ACAT(c(D[D$gene==a[d], 5])[1:58], c(rep(1/(col-2),(col-2))))
        colnames(D)[13] <- "ACAT_SUM"
        D[D$gene==a[d], 14]<- ACAT(c(D[D$gene==a[d], 6])[1:58], c(rep(1/(col-2), (col-2))))
        colnames(D)[14] <- "ACAT_SSU"
        D[D$gene==a[d], 15]<- ACAT(c(D[D$gene==a[d], 7])[1:58], c(rep(1/(col-2), (col-2))))
        colnames(D)[15] <- "ACAT_minP"
        D[D$gene==a[d], 16]<- ACAT(c(D[D$gene==a[d], 8])[1:58], c(rep(1/(col-2), (col-2))))
        colnames(D)[16] <- "ACAT_ACAT1"
        D[D$gene==a[d], 17]<- ACAT(c(D[D$gene==a[d], 9])[1:58], c(rep(1/(col-2), (col-2))))
        colnames(D)[17] <- "ACAT_ACAT2"
        DD <- rbind(DD, D[D$gene==a[d],][1,])
        }
    }

      for (n in 0:(row-1)){    
		for (m in 1:col){
            if (n==0){
                A1 <- D[which(D$TWAS_asy<(0.05/20000)), ]        #A1 sum test significant 
                a1 <- append(a1, as.character(A1$gene))		   #a1 is the sum test significant gene list
           							  
                E1 <- E[which(E$TWAS_asy<(0.05/20000)), ]          #Sum test without weight sum test significant
                e1 <- append(e1, as.character(E1$gene))   	   #e1 is the gene list of without weight significant 								   	   #sum test
         							
                D1 <- D[which(D$ACAT_SUM<0.05/20000), ]		   #ACAT_SUM list
                d1 <- append(d1, as.character(D1$gene))            #

		    list1 <- rbind(list1, E1)
		    list2 <- rbind(list2, D1)
                
            } else if (n==1){
                A1 <- D[which(D$SSU_asy<(0.05/20000)), ]       
                a2 <- append(a2, as.character(A1$gene))
                E2 <- E[which(E$SSU_asy<(0.05/20000)), ]
                e2<- append(e2, as.character(E1$gene))
                D2 <- D[which(D$ACAT_SSU<(0.05/20000)),]
                d2 <- append(d2, as.character(D1$gene))
		    list3 <- rbind(list3, E2)
		    list4 <- rbind(list4, D2)
		    D1 <- D2
                
            } else if (n==2){
                A1 <- D[which(D$minP<(0.05/20000)), ]
                a3 <- append(a3, as.character(A1$gene))
                E3 <- E[which(E$minP<(0.05/20000)), ]
                e3<- append(e3, as.character(E1$gene))
                D3 <- D[which(D$ACAT_minP<(0.05/20000)),]
                d3 <- append(d3, as.character(D1$gene))
		    list5 <- rbind(list5, E3)
		    list6 <- rbind(list6, D3)
		    D1 <- D3
                
            } else if (n==3){
                A1 <- D[which(D$ACAT1<(0.05/20000)), ]
                a4 <- append(a4, as.character(A1$gene))
                E4 <- E[which(E$ACAT1<(0.05/20000)), ]
                e4 <- append(e4, as.character(E1$gene))
                D4 <- D[which(D$ACAT_ACAT1<(0.05/20000)),]
                d4 <- append(d4, as.character(D1$gene))
		    list7 <- rbind(list7, E4)
		    list8 <- rbind(list8, D4)
		    D1 <- D4
                
            } else {
                A1 <- D[which(D$ACAT2<(0.05/20000)), ]
                a5 <- append(a5, as.character(A1$gene))
                E5 <- E[which(E$ACAT2<(0.05/20000)), ]
                e5<- append(e5, as.character(E1$gene))
                D5 <- D[which(D$ACAT_ACAT2<(0.05/20000)),]
                d5 <- append(d5, as.character(D1$gene))
		    list9 <- rbind(list9, E5)
		    list10 <- rbind(list10, D5)
		    D1 <- D5
            }

            k <- col*n*3+(m-1)*3+1
    if (m<col){
            ind[k] <- ind[k]+nrow(A1[A1$weight_mod == m,])
            ind[k+1] <- ind[k+1]+nrow(A1[A1$weight_mod == m & as.numeric(A1$pVal1) < 5e-8,])
            ind[k+2] <- ind[k+2]+nrow(A1[A1$weight_mod == m & as.numeric(A1$pVal2) < 5e-8,])
    }else{  
		
            ind[k] <- ind[k]+(nrow(D1)/(col-2))
            ind[k+1] <- ind[k+1]+(sum(as.numeric(!is.na(D1$pVal1))< 5e-8)/(col-2))
            ind[k+2] <- ind[k+2]+(sum(as.numeric(!is.na(D1$pVal2))< 5e-8)/(col-2))
    }

    out <- rbind(out, A1)
}
}
}
	write.table(out, file="/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/ResultwithpVal3.txt", col.names=TRUE)

for (j in 1:(row*col)){
	l<- 3*j-2
	m<- 3*j-1
	n<- 3*j
	content[j] <- as.character(paste(ind[l],"/", ind[m],"/",ind[n], sep=''))
}
summary <-matrix(content, nrow=row, ncol=col, byrow=TRUE, dimnames=list(c("SUM", "SSU", "minP", "ACAT1", "ACAT2"), c(paste("BETAVOL", 1:(col-2), sep=""), "No_Weight", "ACAT1_58")))
summary



##########################################
#####     CTG.rds data summary     #####
##########################################


saveRDS(content, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/contentCTG.rds")
saveRDS(list1, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list1_CTG.rds")
saveRDS(list2, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list2_CTG.rds")
saveRDS(list3, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list3_CTG.rds")
saveRDS(list4, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list4_CTG.rds")
saveRDS(list5, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list5_CTG.rds")
saveRDS(list6, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list6_CTG.rds")
saveRDS(list7, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list7_CTG.rds")
saveRDS(list8, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list8_CTG.rds")
saveRDS(list9, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list9_CTG.rds")
saveRDS(list10, "/gpfs/home/sj14m/summerproject/project3_revised/result_1/CTG/list10_CTG.rds")


#library(VennDiagram)
#grid.newpage()
#draw.triple.venn(area1 = 3, area2 = 32, area3= 3, n12 = 3, n23 = 0, n13 = 0,
#n123 =0, category = c("No_weight",  "SIWAS", "TWAS"), lty = "blank",
#fill = c("green", "coral", "cornflower blue"))


