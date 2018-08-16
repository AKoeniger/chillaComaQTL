# QTL mapping of cold tolerance in D. ananassae from Bangkok
# skript modified from S. Arif ("CheckvariantcallDist.R")
# markers generated with Stacks v. 1.45 ("QTL_Mapping_pipeline4")

                                 #### Version 6 ####

                                 
                                 
#SNPs: second corrections applied; filtered for heterozygozity                                 
                                 
                                 
setwd("/media/koeniger/Elements1/RAD-Seq/QTL_Mapping_pipeline_4/stacks_out/cor_stacks_2")

#Step 1 Filtering on expected genotype and allele frequencies

#We will use markers that were called in 80% (75/94) of the progeny
snps6 <- read.table("batch_1.genotypes_75.onemap.txt", skip=3, na.strings=c("-","--"), stringsAsFactors = F)

#lowercase h is the original genotype call for heterozygote and captial H is the correction applied by the genotypes utility in stacks
#convert everything to uppercase
numsnps6= nrow(snps6)  #3092
numind6=ncol(snps6)-1  #94

library("stringr")
snps6[,2:ncol(snps6)] <- data.frame(lapply(snps6[,2:ncol(snps6)], function(v) {
  return(as.character(toupper(v)))
}))

head(snps6)

#make a dataframe to each marker name, heterozygity (total het call over all calls) and allele frequency drift (absoulte different in a-b allele fre)
marker_summ6 <- data.frame(markers=snps6$V1, het=numeric(length=numsnps6), drift=numeric(length=numsnps6), PropMiss=numeric(length=numsnps6),
                           numAA=numeric(length=numsnps6), numAB=numeric(length=numsnps6), numBB=numeric(length=numsnps6))


for (i in 1:numsnps6){
  #do the genotype counts
  x <- table(t(snps6[i,2:ncol(snps6)]), useNA = "always")
  #calculate heterozygote frequencies
  marker_summ6$het[i] <- x[3]/sum(x[1:3])
  #calculate the allele frequency drift
  #1. calculate allele frequency of one allele
  a <- ((2*x[1])+x[3])/(2*sum(x[1:3]))
  #2. then the absolute differenece of it's frequency from 0.5
  marker_summ6$drift[i] <- abs(0.5-a)
  #finally addinformation on proportion of missing indivduals for the tag
  marker_summ6$PropMiss[i] <- x[4]/numind6
  marker_summ6$numAA[i] <- x[1]
  marker_summ6$numAB[i] <- x[3]
  marker_summ6$numBB[i] <- x[2]
}

head(marker_summ6)

# noow use the criter for het >= 0.15 and <= 0.35 and drift <=0.1 to eliminate crazy marker
# experiment with marker het between 0.1-0.4 and drift upto 0.2
marker_clean6 <- marker_summ6[which((marker_summ6$het >= 0.15 & marker_summ6$het <= 0.35) & marker_summ6$drift <= 0.1),]
nrow(marker_clean6) # 1400
head(marker_clean6)
#Alternatively using segregation distortion tests
#The following can be use Instead to check against a segregation distortion of 0.375:0.25:0.375
#Calculate all the p-values first
#marker_seg6 <- data.frame(markers=snps6$V1, pvals=numeric(length=numsnps6), adjPvals=numeric(length=numsnps6))
#exp_p <- c(0.375,0.25,0.375)
#for (i in 1:numsnps6){
  #avoid missing genotypes for markers
#  if(any(is.na(c(marker_summ6$numAA[i],marker_summ6$numAB[i],marker_summ6$numBB[i])))) {
#    marker_seg6$pvals[i] = NA
#  }
#  else{
    #print(c(marker_summ$numAA[i],marker_summ$numAB[i],marker_summ$numBB[i]))
#    test_results6 <-chisq.test(c(marker_summ6$numAA[i],marker_summ6$numAB[i],marker_summ6$numBB[i]), p=exp_p  )
#    marker_seg3$pvals[i] <- test_results3$p.value
#  }
#}

#Adjust the p-values usign Bejnaminin-Hochberg, could use bonferroni to be more conservative
#marker_seg3$adjPvals <- p.adjust(marker_seg3$pvals, method="BH")

#Filter markers based on NAs and padj <0.05  
#marker_clean_seg3 <- marker_summ6[which(marker_seg3$adjPvals > 0.05),]
#nrow(marker_clean_seg3)


# ----------------------------------------------------------------------------------------------------------

#Step 2: Prepare an input file for R/QTL from the input file

#we'll need the phenotypes
pheno <- read.csv("../../Analysis_in_R/pheno_file.csv")

#construct the first column, the pheno colun, for the read.cross file
pheno_cross <- c("chill_coma", " ", " ", pheno$pheno)

#now we need scaffold and positions information, this is stored "batch_1.catalog.tags.tsv.gz"
# output file from the stacks genotypes run
#system("gunzip batch_1.catalog.tags.tsv.gz")

mark_info6 <- read.table("batch_1.catalog.tags.tsv", header=F)
head(mark_info6)
str(mark_info6)
#catalog id from the one map file is in the 3rd column (it is the row number), we want to recover the 4th and the 5th col (scaffold position in kb)
#if filtered above
cat_ids6 <- as.numeric(str_sub(marker_clean6$markers, 2, -1))
#if NOTfiltered
#cat_ids6 <- as.numeric(str_sub(snps6$V1, 2, -1))

marker_header6 <- mark_info6[cat_ids6, 3:5]
str(marker_header6)
#convert all genotypes calls to same case - already done on lines 24-26
#recoded <- apply(snps[2:ncol(snps)],2,tolower)

#replace the marker columns in the snps file with the marker name, scaffold and position info
#if filtered
new_snps6 <- cbind(marker_header6, snps6[which(snps6$V1 %in% marker_clean6$markers), 2:ncol(snps6)])
#if NOT filtered
#new_snps6 <- cbind(marker_header6, snps6[, 2:ncol(snps6)])
head(new_snps6)
#reorder the columns by scaffold than positions
new_snps6 <- new_snps6[with(new_snps6, order(V4, V5)), ]
#transpose the new_snps and combine with the phen0o_cross
#save the file as .csv
crossfile6 <- cbind(pheno_cross, t(new_snps6))
str(crossfile6)
#getrid of the rownames and colnames
rownames(crossfile6)=NULL
colnames(crossfile6)=NULL
#save the file
write.table(crossfile6, "../../Analysis_in_R/Cold_r75_crossfile6.csv", quote=F, sep=",",col.names=F, row.names =F)
#Replace the long "scaffolds_" with "s" in the file 
system("sed -i 's/scaffold_/s/g' ../../Analysis_in_R/Cold_r75_crossfile6.csv")

#-------------------------------------------------------------------------------------------------------
#  Step 3: Playing around in R/Qtl

#Note this assumes the cross file is in the same directory and the working directory is the the one where the script is stored
#change accordingly

cold6 <- read.cross(format="csv", ".", "Cold_r75_crossfile6.csv", genotypes=c("A",  "H", "B"), estimate.map = F)
#I have not specified F.Gen as i will not estimating a genetic map
#you may want to repeat the analysis with cleaned version below and genetic distances estimated assuming F.gen=5, for comparison
#you can convert your cross object to F.gen =t with estimated map as follows

#cold6<- convert2bcsft(cold6, F.gen = 5)

summary(cold6)
#rescale map from kb to mb
cold6 <- rescalemap(cold6, 1e-6) 

#----------------------------------------------------------------------------------------------------------------------------
#Check genotyping quality after filtering
geno.image(cold6)
#in an ideal world this should only be contiguous blocks of red,green (homozygous) and some blue (heterozygous)

#We still need to clean this up a bit
#R/qtl does deal with genotyping error but as should be obvious from the plot above
#there is still quite a bit of error here.

#We could do two things:
#(i)get rid of problematic indivdiuals
#(ii) get get rid od problematic markers in singl individuals.

#In this cleaning approach both filtering methods rely on detecting unexpected double crossover events

#The following is adapted from pages 33-39 from this tutorial
#http://www.rqtl.org/tutorials/geneticmaps.pdf

#Lets remove problematic individuals first as they will make detecting genotyping errors more tenable later

#Plot the number of double crossovers per individual
plot(countXO(cold6), ylab="Number of crossovers")
#most lines have double crossovers betweem 30-150 but ehre are clear outliers. Lets remove these for now
cold6 <- subset(cold6, ind=(countXO(cold6) < 200))
#We lose 8 lines
nind(cold6)   #87
#Did this qualitatively improve our genotyping?
geno.image(cold6)
#I'm not sure if retaining the indvidual with missing data is helpful or not, it is not helpful to have the individual for the nexy
#step (identifying potential genotyping errors ) so I will remove that individual
cold6 <- subset(cold6, ind=(ntyped(cold6)>700))   #68

#Now can We identify potentially erroneous double crossovers assumign a genotyping error rate of 1%
cold6 <- calc.errorlod(cold6, error.prob=0.01)

#We will remove any marker with an abritrary large LOD score, for e.g. greater than 5 (cutoff)
nrow(toperr <- top.errorlod(cold6, cutoff=5))

#remove markers with LOD scores greater than five and store it as a new cross object
cold6c <-cold6
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  cold6c$geno[[chr]]$data[id, mar] <- NA
}

#NOTE: I would recommend looking at the help page for calc.errorlod to find out what is happening and also looking at the function
#clean.geno() as an alternative/complimentary

summary(cold6)
summary(cold6c)
#We haven't lost a lot of markers but we have lost 9 indiviuals #8

#Has the genotyping qualtiatively improved?
geno.image(cold6c)



#------------------------------------------------------------------------------------------------------------------------------------------
#YOU COULD IGNORE THE FOLLOWING IF YOU LIKE, it is not directly related to above it is just checking the genetic map
#estimating the map using the est.map function
newmap <- est.map(cold6c, n.cluster=6)
plotMap(cold6c, newmap, alternate.chrid=T)
plotMap(cold6c, newmap, alternate.chrid=T)
#There seem to be some potentially problematic markers on s13337  (and maybe the terminal one on s13117)
#Lets calculate pairwise recombination frequency
#Can we remove the problematic (unlinked markers) on s13337
cold6c <- est.rf(cold6c)
#and plot it 
plotRF(cold6c, alternate.chrid=T) 
checkAlleles(cold6c)
#No apparent problems.

#cold6<-drop.markers(cold6, c("19153", "19632"))
#cold6 <- est.rf(cold6)
#plotRF(cold6, alternate.chrid=T) 
#PLEASE STOP IGNORING NOW
#----------------------------------------------------------------------------------------------------------------------------------------------

#You could not continue with your QTL mapping approach as before but try with the cold6.clean cross
#Note I've put step to 0 so that it does not estimate any intermarker distances as this does not make much sense
#if you haven't estimated a genetic map. I would recommend that in additional to the analysis with no estimated genetic map
#you could do one with a genetic map estimated assuming F.gen=5, you can covert the current cross to F5 and simulatenously estimate the map as follows:

cold6.c.F5<- convert2bcsft(cold6c, F.gen = 5)
#Have a look at the genetic map before proceeding? Does it look reasonable?
plotMap(cold6c.F5, alternate.chrid = T)
summaryMap(cold6c.F5)
#         n.mar length ave.spacing max.spacing
# overall  1400  962.0         0.7        55.5
lg<-formLinkageGroups(cold6c.F5, max.rf = 0.35, min.lod = 6)
table(lg[,2])
#  1   2   3   4   5   6   7   8   9 
#784 200 180 125  79  27   2   2   1 

############################ scanone

cold6c.F5 <- calc.genoprob(cold6c.F5, step = 1)   # set back to 1 since we are using the genetic map now
result_em_cold6c.F5 <- scanone(cold6c.F5, pheno.col=1, method= "em")
result_hk_cold6c.F5 <- scanone(cold6c.F5, pheno.col=1, method= "hk")
result_ehk_cold6c.F5 <- scanone(cold6c.F5, pheno.col=1, method= "ehk")
# find significance threshold
perm_em_cold6c.F5<-scanone(cold6c.F5, method="em", n.perm = 1000)
perm_hk_cold6c.F5<-scanone(cold6c.F5, method="hk", n.perm = 1000)
perm_ehk_cold6c.F5<-scanone(cold6c.F5, method="ehk", n.perm = 1000)
summary(perm_em_cold6c.F5)
#lod
#5%  4.66
#10% 4.18
summary(perm_hk_cold6c.F5)
#lod
#5%  4.00
#10% 3.58
summary(perm_ehk_cold6c.F5)  
#lod
#5%  153
#10% 153
summary(result_em_cold6c.F5, threshold=4.1, df=TRUE)
#                 chr     pos  lod
#22061         s13337  0.0839 6.81
#cs13340.loc80 s13340 80.0531 5.13
summary(result_hk_cold6c.F5, threshold=3.5, df=TRUE)
#               chr     pos  lod
#22061         s13337  0.0839 6.82
#cs13340.loc80 s13340 80.0531 5.27
#23306        s13340 11.450 5.86


################### scantwo

cold6c.F5 <- calc.genoprob(cold6c.F5, step = 1) 
result_scan2_hk_cold6c.F5 <- scantwo(cold6c.F5, pheno.col=1, method= "hk", clean.output=TRUE)       
perm_scan2_hk_cold6c.F5<- scantwopermhk(cold6c.F5, n.perm = 1000)  # takes FOREVER
scan2_hk_cold6c.F5_summary <- summary(result_scan2_hk_cold6c.F5, perms = perm_scan2_hk_cold6c.F5, alpha = 0.2, pvalues=T)

#                pos1f pos2f lod.full pval lod.fv1 pval lod.int pval     pos1a pos2a lod.add pval lod.av1 pval
#cs12916:cs13340 92.16  74.1     10.1 0.03    4.86 0.99    1.63    1     92.16  80.1     8.5    0    3.23 0.19
#cs13337:cs13340  7.08  30.1     12.6 0.00    5.77 0.71    1.24    1      7.08  30.1    11.4    0    4.54 0.02

write.table(scan2_hk_cold6_summary, "../../Analysis_in_R/version6c.F5_scan2.csv", row.names = F)

summary(perm_scan2_hk_cold6c.F5, alpha = 0.05)
#   full  fv1  int  add  av1  one
#5% 9.65 7.36 6.33 6.61 3.68 3.86

########################### find markers

find.marker(cold6c.F5, chr = "s13340", pos =80.1 )    #QTL3
#"23306"
find.marker(cold6c.F5, chr = "s13340", pos =30.1)     #QTL2
#"27008"
find.marker(cold6c.F5, chr = "s12916", pos =92.16)    #QTL4
#" 2791"
find.marker(cold6c.F5, chr = "s13337", pos =7.08)     #QTL1
# "20419"


############################### Plots
# scanone
plot(result_em_cold6c.F5, main="Scanone EM Algorithm",col="#08519c", alternate.chrid=TRUE)
abline(h=summary(perm_em_cold6c.F5, alpha=0.05))
#scantwo
plot(result_scan2_hk_cold6c.F5, main="Test for interaction vs add model", alternate.chrid=TRUE, upper="int", lower = "add")



######## MQM by imputation

#Multiple QTL model:
#Multiple QTL model:
cold6c.F5<- sim.geno(cold6c.F5, step = 1, n.draws = 32, error.prob = 0.0001)
qtl_model6c <- sim.geno(cold6c.F5, step = 1, n.draws = 32, error.prob = 0.0001)
qtl6c <- makeqtl(qtl_model6c,chr=c("s13337", "s13340","s13340", "s12916"), pos=c(7.08, 30.1,80.1,92.16)) # candidate regions
qtl6c                                                                                                    

#QTL object containing imputed genotypes, with 32 imputations. 
#
#name    chr     pos n.gen
#Q1  s13337@7.1 s13337  7.0839     3
#Q2 s13340@30.1 s13340 30.1412     3
#Q3 s13340@80.1 s13340 80.0531     3
#Q4 s12916@92.2 s12916 92.1596     3

# NOTE:  QTL3: using the "additive position" at 80.1 cM gives better results then the "full position" at 74,1cM

plot(qtl6, alternate.chrid=T, main= "Genetic map version 6c.F5") # Map of the chromosomes with candidate loci

# possible model without interaction:

out.fqa6c <- fitqtl(qtl_model6c, qtl=qtl6c, get.ests=T, pheno.col=1,formula=y~Q1+Q2+Q3+Q4)
mqm_imp_cold6c_summary <- summary(out.fqa6c)
mqm_imp_cold6c_summary
# fitqtl summary
# 
# Method: multiple imputation 
# Model:  normal phenotype
# Number of observations : 86 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 + Q2 + Q3 + Q4 
# 
# df       SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model  8 10740.64 1342.57957 16.31011 58.24629 4.687362e-13 5.758949e-12
# Error 77  7699.40   99.99221                                            
# Total 85 18440.04                                                       
# 
# 
# Drop one QTL at a time ANOVA table: 
#   ----------------------------------  
#   df Type III SS   LOD   %var F value Pvalue(Chi2) Pvalue(F)    
# s13337@7.1   2        2768 5.735 15.008  13.839        0.000  7.34e-06 ***
#   s13340@30.1  2        1651 3.629  8.956   8.258        0.000  0.000563 ***
#   s13340@80.1  2        1246 2.801  6.757   6.230        0.002  0.003106 ** 
#   s12916@92.2  2        1090 2.473  5.911   5.450        0.003  0.006113 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# Estimated effects:
#   -----------------
#   est       SE      t
# Intercept     46.2261   1.4858 31.112
# s13337@7.1a    7.0847   1.3692  5.175
# s13337@7.1d   -1.7591   2.4670 -0.713
# s13340@30.1a  -0.0689   1.4843 -0.046
# s13340@30.1d -10.6973   2.6662 -4.012
# s13340@80.1a   5.1059   1.4984  3.408
# s13340@80.1d  -1.2126   2.9654 -0.409
# s12916@92.2a  -3.9349   1.3576 -2.899
# s12916@92.2d  -4.1797   2.6314 -1.588

mqm_imp_cold6c.F5_tab <- mqm_imp_cold6c.F5_summary$result.drop 
write.table(mqm_imp_cold6c.F5_tab, file="mqm_imp_cold6c.F5.csv", sep="\t", quote=FALSE, row.names=TRUE) 

