install.packages("bigrquery")
install.packages("httpuv")
install.packages("ggplot2")
install.packages("reshape")
install.packages("scales")
install.packages("dplyr")
BiocManager::install('gargle')
library(bigrquery)
library(httpuv)
library(ggplot2)
library(reshape)
library(scales)
library(dplyr)
BiocManager::install('gdsfmt')
library(gdsfmt)
library(gargle)
install.packages('readtext')
library("BiocParallel")
register(MulticoreParam(40)) #Use all the available cores, its a large genome
library(SNPRelate)
library(ggplot2)
library(gtools)
BiocManager::install('rjson')
BiocManager::install('bigsnpr')
devtools::install_github("privefl/bigsnpr")
library(bigsnpr)
require(gridExtra)
library(pcadapt)
#library(tm)

list_datasets('isb-cgc')
list_tables('isb-cgc','genome_reference')
pr1='isb-cgc-psm'
sql1="select * from 'isb-cgc.genome_reference.ClinVar_20180401_GRCh38'"
tbl1=query_exec(pr1,sql1,use_legacy_sql = F)
setwd('/shire/zillur/exome_ancestry/captured_regions/')
sql3="SELECT Chromosome, Start_Position, End_Position, Strand, Variant_Type FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation`"
tbl3=query_exec(sql3,pr1,use_legacy_sql = F,max_pages = Inf)

sql4="SELECT Chromosome,Start_Position,End_Position,Genome_Change FROM `isb-cgc.ccle_201602_alpha.Mutation_calls`"
tbl4=query_exec(sql4,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl4,'ccle_201602_alpha_mutation.csv')
sql5="SELECT Chromosome,Start_Position,End_Position,Genome_Change  FROM `isb-cgc.tcga_201607_beta.Somatic_Mutation_calls`"
tbl5=query_exec(sql5,pr1,use_legacy_sql = F,max_pages = Inf)
sql6="SELECT chromosome, start_pos, end_pos,num_probes FROM `isb-cgc.TCGA_hg19_data_v0.Copy_Number_Segment_Masked`"
tbl6=query_exec(sql6,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl6,row.names = F,'tcga_hg19_CNS.csv')
sql7="SELECT Chromosome, Start_Position,End_Position,Genome_Change FROM `isb-cgc.TCGA_hg19_data_v0.Somatic_Mutation_DCC`"
tbl7=query_exec(sql7,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl7,row.names = F,'tcga_hg19_somatic_dcc.csv')
sql8="SELECT Chromosome,Start_Position,End_Position, Variant_Type FROM `isb-cgc.TCGA_hg19_data_v0.Somatic_Mutation_MC3`"
tbl8=query_exec(sql8,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl8,row.names = F,'tcga_hg19_somatic_mc3.csv')
sql9="SELECT chromosome,start_pos,end_pos,num_probes FROM `isb-cgc.TCGA_hg38_data_v0.Copy_Number_Segment_Masked`"
tbl9=query_exec(sql9,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl9,row.names = F,'tcga_hg38_CNS.csv')
sql10="SELECT chromosome,start_pos,end_pos,num_probes FROM `isb-cgc.TCGA_hg38_data_v0.Copy_Number_Segment_Masked_r14`"
tbl10=query_exec(sql10,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl10,row.names = F,'tcga_hg38_CNS_r14.csv')
sql11="SELECT Chromosome,Start_Position,End_Position,Strand FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR10`"
tbl11=query_exec(sql11,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl11,row.names = F,'tcga_hg38_somatic_dr10.csv')
sql12="SELECT Chromosome,Start_Position,End_Position,Strand FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR6`"
tbl12=query_exec(sql12,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl12,row.names = F,'tcga_hg38_somatic_dr6.csv')
sql13="SELECT Chromosome,Start_Position,End_Position,Strand FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR7`"
tbl13=query_exec(sql13,pr1,use_legacy_sql = F,max_pages = Inf)
write.csv(tbl13,row.names = F,'tcga_hg38_somatic_dr7.csv')
sql14="SELECT  reference_name,start_position,reference_bases   FROM `bigquery-public-data.human_variant_annotation.ensembl_esp6500_aa_hg19_release93` "
tbl14=query_exec(sql14,pr1,use_legacy_sql = F,max_pages = Inf)
sql15="SELECT reference_name,start_position,end_position,reference_bases FROM `bigquery-public-data.human_genome_variants.simons_genome_diversity_project_sample_variants` "
tbl15=query_exec(sql15,pr1,use_legacy_sql = F,max_pages = Inf)
sql16="SELECT reference_name,start_position,end_position,reference_bases FROM `bigquery-public-data.human_genome_variants.platinum_genomes_deepvariant_variants_20180823`"
tbl16=query_exec(sql16,pr1,use_legacy_sql = F,max_pages = Inf)


sq_sm1="SELECT Chromosome,Start_Position,End_Position,Genome_Change  FROM `isb-cgc.TCGA_hg19_data_v0.Somatic_Mutation_DCC` WHERE Mutation_Status='Somatic'"
t_sm1=query_exec(sq_sm1,pr1,use_legacy_sql = F,max_pages = Inf)
sq_sm2="SELECT Chromosome,Start_Position,End_Position,STRAND  FROM `isb-cgc.TCGA_hg19_data_v0.Somatic_Mutation_MC3` WHERE SOMATIC='1'"
t_sm2=query_exec(sq_sm2,pr1,use_legacy_sql = F,max_pages = Inf)
sq_sm3="SELECT Chromosome,Start_Position,End_Position,Strand  FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation` WHERE Mutation_Status='Somatic'"
t_sm3=query_exec(sq_sm3,pr1,use_legacy_sql = F,max_pages = Inf)
sq_sm4="SELECT Chromosome,Start_Position,End_Position,Strand  FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR10` WHERE Mutation_Status='Somatic'"
t_sm4=query_exec(sq_sm4,pr1,use_legacy_sql = F,max_pages = Inf)
sq_sm5="SELECT Chromosome,Start_Position,End_Position,Strand FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR6` WHERE Mutation_Status='Somatic'"
t_sm5=query_exec(sq_sm5,pr1,use_legacy_sql = F,max_pages = Inf)
sq_sm6="SELECT Chromosome,Start_Position,End_Position,Strand FROM `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR7` WHERE Mutation_Status='Somatic'"
t_sm6=query_exec(sq_sm6,pr1,use_legacy_sql = F,max_pages = Inf)



b1=read.csv('1321.maf5.indep0.1.bim',sep = '\t',header = F)

fun1=function(df,snps2){
df=data.frame(lapply(df, as.character),stringsAsFactors = F)
df=df[(nchar(df$V5)==1 & nchar(df$V6)==1),]
df$snps2=paste(df$V5,df$V6)
return(df)
}
bs1=list(b10,b15,b8,b6,b3,b5,b1,b7,b11,b4,b9,b13,b2,b16,b12,b14)
for (i in 1:16) {bs1[[i]] <- fun1(bs1[[i]],snps2)}
pn1=c('gwas','tcga hg38 somatic mutation','exome','clinvar','bq ccle mutation','ccle hc','1321','broca','cpg hg38','rand1_1000','fone','rand1_300','aims','tst','rand1_150','rand1_50')

png('snps_plot/snp_bar_7.png',width = 1600,height = 1600)
par(mfrow=c(4,4))
for (i in 1:16){barplot(table(bs1[[i]]$snps2)/sum(table(bs1[[i]]$snps2)),col = rainbow(12),main=paste(nrow(bs1[[i]]), "SNPs in",pn1[i], sep=" "))}
dev.off()

cv1=read.csv('all.cv.error2',sep = ':',header = F)
cv2=cve1
cv2$V2=gsub("[^0-9.]", "", cv2$V2)
cv2$V1=gsub("log[0-9.]", "log", cv2$V1)
cv2$V1=gsub('0.out','.out',cv2$V1)

dfs=split(cv2,cv2$V1)
for (i in 1:16){dfs[[i]]$V2=as.numeric(as.character(dfs[[i]]$V2))}
for (i in 1:16){dfs[[i]]=dfs[[i]][order(dfs[[i]]$V2),]}
for (i in 1:16){dfs[[i]]$V2=sort(dfs[[i]]$V2)}
str(dfs)
cv3=rbind_list(dfs)
colnames(cv3)=c('Panels','K','CV_error')
png('snps_plot/cv_err_2.png')
ggplot(data=cv3, aes(x=K, y=CV_error, group=Panels, shape=Panels, color=Panels)) +
  geom_line() +
  geom_point()
dev.off()

png('snps_plot/test_pca_af_19.png',width = 1661,height = 2332)
x1321 <- read.pcadapt(bed_r1321, type = "bed")
x1321 <- pcadapt(input = x1321, K = 20)
xaims <- read.pcadapt(bed_aims, type = "bed")
xaims <- pcadapt(input = xaims, K = 20)
xbqccle <- read.pcadapt(bed_bq_ccle_mutation, type = "bed")
xbqccle <- pcadapt(input = xbqccle, K = 20)
xbroca <- read.pcadapt(bed_broca, type = "bed")
xbroca <- pcadapt(input = xbroca, K = 20)
xccle <- read.pcadapt(bed_ccle, type = "bed")
xccle <- pcadapt(input = xccle, K = 20)
xclinvar <- read.pcadapt(bed_clinvar, type = "bed")
xclinvar <- pcadapt(input = xclinvar,K=20)
xcpg <- read.pcadapt(bed_cpghg38, type = "bed")
xcpg <- pcadapt(input = xcpg, K = 20)
xexome <- read.pcadapt(bed_exome, type = "bed")
xexome <- pcadapt(input = xexome, K = 20)
xfone <- read.pcadapt(bed_fone, type = "bed")
xfone <- pcadapt(input = xfone, K = 20)
xgwas <- read.pcadapt(bed_gwas, type = "bed")
xgwas <- pcadapt(input = xgwas, K = 20)
xrand1_1000 <- read.pcadapt(bed_rand1_1000, type = "bed")
xrand1_1000 <- pcadapt(input = xrand1_1000, K = 20)
xrand1_150 <- read.pcadapt(bed_rand1_150, type = "bed")
xrand1_150 <- pcadapt(input = xrand1_150, K = 20)
xrand1_300 <- read.pcadapt(bed_rand1_300, type = "bed")
xrand1_300 <- pcadapt(input = xrand1_300, K = 20)
xrand1_50 <- read.pcadapt(bed_rand1_50, type = "bed")
xrand1_50 <- pcadapt(input = xrand1_50, K = 20)
xtcga_somatic <- read.pcadapt(bed_tcga_somatic, type = "bed")
xtcga_somatic <- pcadapt(input = xtcga_somatic, K = 20)
xtst <- read.pcadapt(bed_tst, type = "bed")
xtst <- pcadapt(input = xtst, K = 20)
xtcga_hg19_somatic_dcc2=read.pcadapt('tcga_hg19_somatic_dcc.maf5.indep0.1.bed',type = 'bed')
xtcga_hg19_somatic_dcc2=pcadapt(input = xtcga_hg19_somatic_dcc2,K=20)
xtcga_hg19_somatic_mc32=read.pcadapt('tcga_hg19_somatic_mc3.maf5.indep.1.bed',type = 'bed')
xtcga_hg19_somatic_mc32=pcadapt(input = xtcga_hg19_somatic_mc32,K=20)
xtcga_hg38_somatic_dr102=read.pcadapt('tcga_hg38_somatic_dr10.maf5.indep0.1.bed',type = 'bed')
xtcga_hg38_somatic_dr102=pcadapt(input = xtcga_hg38_somatic_dr102,K=20)
xtcga_hg38_somatic_dr62=read.pcadapt('tcga_hg38__somatic_dr6.maf5.indep0.1.bed',type = 'bed')
xtcga_hg38_somatic_dr62=pcadapt(input = xtcga_hg38_somatic_dr62,K=20)
xtcga_hg38_somatic_dr72=read.pcadapt('tcga_hg38_somatic_dr7.maf5.indep0.1.bed',type = 'bed')
xtcga_hg38_somatic_dr72=pcadapt(input = xtcga_hg38_somatic_dr72,K=20)

pl45=plot(xtcga_hg19_somatic_dcc2, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl46=qplot(xtcga_hg19_somatic_dcc2$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle(paste('hg19_somatic_dcc', length(xtcga_hg19_somatic_dcc2$af), 'SNPs'))+theme(plot.title = element_text(size = 22))
pl47=plot(xtcga_hg38_somatic_dr62, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl48=qplot(xtcga_hg38_somatic_dr62$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle(paste('hg38_somatic_dr6', length(xtcga_hg38_somatic_dr62$af), 'SNPs'))+theme(plot.title = element_text(size = 22))
pl49=plot(xtcga_hg19_somatic_mc32, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl50=qplot(xtcga_hg19_somatic_mc32$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle(paste('hg19_somatic_mc3', length(xtcga_hg19_somatic_mc32$af), 'SNPs'))+theme(plot.title = element_text(size = 22))
pl51=plot(xtcga_hg38_somatic_dr102, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl52=qplot(xtcga_hg38_somatic_dr102$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle(paste('hg38_somatic_dr10', length(xtcga_hg38_somatic_dr102$af), 'SNPs'))+theme(plot.title = element_text(size = 22))
pl53=plot(xtcga_hg38_somatic_dr72, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl54=qplot(xtcga_hg38_somatic_dr72$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle(paste('hg38_somatic_dr7', length(xtcga_hg38_somatic_dr72$af), 'SNPs'))+theme(plot.title = element_text(size = 22))




pl1=plot(x1321, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl2=plot(xaims, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl3=plot(xbqccle, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl4=plot(xbroca, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl5=plot(xccle, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl6=plot(xclinvar, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl7=plot(xcpg, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl8=plot(xexome, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl9=plot(xfone, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl10=plot(xgwas, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl11=plot(xrand1_1000, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl12=plot(xrand1_150, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl13=plot(xrand1_300, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl14=plot(xrand1_50, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl15=plot(xtcga_hg38_somatic_dr10, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_somatic_dr10', length(xtcga_hg38_somatic_dr10$af), 'SNPs'))
pl16=plot(xtst, option = "scores", pop = pop_code)+ggtitle('')+theme(plot.title = element_text(size = 0))
pl17=plot(xtcga_hg19_cns_maf, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg19_CNS', length(xtcga_hg19_cns_maf$af), 'SNPs'))
pl18=plot(xtcga_hg38_cns_maf, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_CNS', length(xtcga_hg38_cns_maf$af), 'SNPs'))
pl19=plot(xtcga_hg19_somatic_dcc, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg19_somatic_dcc', length(xtcga_hg19_somatic_dcc$af), 'SNPs'))
pl20=plot(xtcga_hg19_somatic_mc3, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg19_somatic_mc3', length(xtcga_hg19_somatic_mc3$af), 'SNPs'))
pl21=plot(xtcga_hg38_somatic_dr6, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_somatic_dr6', length(xtcga_hg38_somatic_dr6$af), 'SNPs'))
pl22=plot(xtcga_hg38_somatic_dr7, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_somatic_dr7', length(xtcga_hg38_somatic_dr7$af), 'SNPs'))
pl23=qplot(xtcga_hg19_cns_maf$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+theme(plot.title = element_text(size = 22))
pl24=qplot(xtcga_hg38_cns_maf$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+theme(plot.title = element_text(size = 22))
pl25=qplot(xgwas$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('gwas 525540 SNPs')+theme(plot.title = element_text(size = 22))
pl26=qplot(xtcga_hg19_somatic_dcc$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+theme(plot.title = element_text(size = 22))
pl27=qplot(xtcga_hg19_somatic_mc3$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+theme(plot.title = element_text(size = 22))
pl28=qplot(xtcga_hg38_somatic_dr6$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+theme(plot.title = element_text(size = 22))
pl29=qplot(xtcga_hg38_somatic_dr10$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+theme(plot.title = element_text(size = 22))
pl30=qplot(xtcga_hg38_somatic_dr7$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+theme(plot.title = element_text(size = 22))
pl31=qplot(xexome$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('exome 46156 SNPs')+theme(plot.title = element_text(size = 22))
pl32=qplot(xclinvar$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('clinvar hg38 17076 SNPs')+theme(plot.title = element_text(size = 22))
pl33=qplot(xbqccle$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('bq ccle 5443 SNPs')+theme(plot.title = element_text(size = 22))
pl34=qplot(xccle$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('ccle 5262 SNPs')+theme(plot.title = element_text(size = 22))
pl35=qplot(x1321$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('1321 4052 SNPs')+theme(plot.title = element_text(size = 22))
pl36=qplot(xbroca$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('broca 2495 SNPs')+theme(plot.title = element_text(size = 22))
pl37=qplot(xcpg$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('cpg hg38 2380 SNPs')+theme(plot.title = element_text(size = 22))
pl38=qplot(xrand1_1000$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('rand1_1000 2177 SNPs')+theme(plot.title = element_text(size = 22))
pl39=qplot(xfone$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('fone 736 SNPs')+theme(plot.title = element_text(size = 22))
pl40=qplot(xrand1_300$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('rand1_300 718 SNPs')+theme(plot.title = element_text(size = 22))
pl41=qplot(xaims$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('AIMS 445 SNPs')+theme(plot.title = element_text(size = 22))
pl42=qplot(xtst$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('tst 365 SNPs')+theme(plot.title = element_text(size = 22))
pl43=qplot(xrand1_150$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('rand1_150 326 SNPs')+theme(plot.title = element_text(size = 22))
pl44=qplot(xrand1_50$af,geom = 'histogram',binwidth=0.01,xlab = 'allele frequency',ylab = 'count')+ggtitle('rand1_50 108 SNPs')+theme(plot.title = element_text(size = 22))
grid.arrange(pl17,pl23,pl18,pl24,pl10,pl25,pl19,pl26,pl20,pl27,pl21,pl28,pl15,pl29,pl22,pl30,pl8,pl31,pl6,pl32,pl3,pl33,pl5,
             pl34,pl1,pl35,pl4,pl36,pl7,pl37,pl11,pl38,pl9,pl39,pl13,pl40,pl2,pl41,pl16,pl42,pl12,pl43,pl14,pl44, ncol=4)
dev.off()


grid.arrange(pl10,pl25,pl45,pl46,pl8,pl31,pl6,pl32,pl47,pl48,pl49,pl50,pl51,pl52,pl53,pl54,pl3,pl33,pl5,
             pl34,pl1,pl35,pl4,pl36,pl7,pl37,pl11,pl38,pl9,pl39,pl13,pl40,pl2,pl41,pl16,pl42,pl12,pl43,pl14,pl44, ncol=4)
dev.off()




xtcga_hg19_cns_maf <- read.pcadapt('tcga_hg19_CNS.maf.indep.bed', type = "bed")
xtcga_hg19_cns_maf <- pcadapt(input = xtcga_hg19_cns_maf, K = 20)
plot(xtcga_hg19_cns_maf, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg19_CNS', length(xtcga_hg19_cns_maf$af), 'SNPs'))


xtcga_hg19_somatic_dcc <- read.pcadapt('tcga_hg19_', type = "bed")
xtcga_hg19_somatic_dcc <- pcadapt(input = xtcga_hg19_somatic_dcc, K = 20)
plot(xtcga_hg19_somatic_dcc, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg19_somatic_dcc', length(xtcga_hg19_somatic_dcc$af), 'SNPs'))

xtcga_hg19_somatic_mc3 <- read.pcadapt('tcga_hg19_somatic_mc3.bed', type = "bed")
xtcga_hg19_somatic_mc3 <- pcadapt(input = xtcga_hg19_somatic_mc3, K = 20)
plot(xtcga_hg19_somatic_mc3, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg19_somatic_mc3', length(xtcga_hg19_somatic_mc3$af), 'SNPs'))

xtcga_hg38_cns_maf <- read.pcadapt('tcga_hg38_CNS.maf.indep.bed', type = "bed")
xtcga_hg38_cns_maf <- pcadapt(input = xtcga_hg38_cns_maf, K = 20)
plot(xtcga_hg38_cns_maf, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_CNS', length(xtcga_hg38_cns_maf$af), 'SNPs'))

xtcga_hg38_cns_r14_maf <- read.pcadapt('tcga_hg38_CNS_r', type = "bed")
xtcga_hg38_cns_r14_maf <- pcadapt(input = xtcga_hg38_cns_maf, K = 20)
plot(xtcga_hg38_cns_maf, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_CNS', length(xtcga_hg38_cns_maf$af), 'SNPs'))


xtcga_hg38_somatic_dr10 <- read.pcadapt('tcga_hg38_somatic_dr10.bed', type = "bed")
xtcga_hg38_somatic_dr10 <- pcadapt(input = xtcga_hg38_somatic_dr10, K = 20)
plot(xtcga_hg38_somatic_dr10, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_somatic_dr10', length(xtcga_hg38_somatic_dr10$af), 'SNPs'))

xtcga_hg38_somatic_dr7 <- read.pcadapt('tcga_hg38_somatic_dr7.bed', type = "bed")
xtcga_hg38_somatic_dr7 <- pcadapt(input = xtcga_hg38_somatic_dr7, K = 20)
plot(xtcga_hg38_somatic_dr7, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_somatic_dr7', length(xtcga_hg38_somatic_dr7$af), 'SNPs'))

xtcga_hg38_somatic_dr6 <- read.pcadapt('tcga_hg38_somatic_dr6.bed', type = "bed")
xtcga_hg38_somatic_dr6 <- pcadapt(input = xtcga_hg38_somatic_dr6, K = 20)
plot(xtcga_hg38_somatic_dr6, option = "scores", pop = pop_code)+ggtitle(paste('tcga_hg38_somatic_dr6', length(xtcga_hg38_somatic_dr6$af), 'SNPs'))



dev.off()

png('snps_plot/manhattan_all_2.png',width = 850,height = 850)
pl1=plot(x1321, option = "manhattan")+ggtitle('P_1321 4052 SNPs')
pl2=plot(xaims, option = "manhattan")+ggtitle('P_AIMS 445 SNPs')
pl3=plot(xbqccle, option = "manhattan")+ggtitle('P_bq ccle mutation 5443 SNPs')
pl4=plot(xbroca, option = "manhattan")+ggtitle('P_broca 2495 SNPs')
pl5=plot(xccle, option = "manhattan")+ggtitle('P_ccle 5262 SNPs')
pl6=plot(xclinvar, option = "manhattan")+ggtitle('P_clinvar hg38 17076 SNPs')
pl7=plot(xcpg, option = "manhattan")+ggtitle('P_cpg hg38 2380 SNPs')
pl8=plot(xexome, option = "manhattan")+ggtitle('P_exome 46156 SNPs')
pl9=plot(xfone, option = "manhattan")+ggtitle('P_fone 736 SNPs')
pl10=plot(xgwas, option = "manhattan")+ggtitle('P_gwas 525540 SNPs')
pl11=plot(xrand1_1000, option = "manhattan")+ggtitle('P_rand1_1000 2177 SNPs')
pl12=plot(xrand1_150, option = "manhattan")+ggtitle('P_rand1_150 326 SNPs')
pl13=plot(xrand1_300, option = "manhattan")+ggtitle('P_rand1_300 718 SNPs')
pl14=plot(xrand1_50, option = "manhattan")+ggtitle('P_rand1_50 108 SNPs')
pl15=plot(xtcga_somatic, option = "manhattan")+ggtitle('P_tcga somatic 96284 SNPs')
pl16=plot(xtst, option = "manhattan")+ggtitle('P_tst 365 SNPs')
grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,pl9,pl10,pl11,pl12,pl13,pl14,pl15,pl16, ncol=4)
dev.off()

png('snps_plot/qq_all_2.png',width = 850,height = 850)
pl1=plot(x1321, option = "qqplot", threshold = 0.1)+ggtitle('P_1321 4052 SNPs')
pl2=plot(xaims, option = "qqplot", threshold = 0.1)+ggtitle('P_AIMS 445 SNPs')
pl3=plot(xbqccle, option = "qqplot", threshold = 0.1)+ggtitle('P_bq ccle mutation 5443 SNPs')
pl4=plot(xbroca, option = "qqplot", threshold = 0.1)+ggtitle('P_broca 2495 SNPs')
pl5=plot(xccle, option = "qqplot", threshold = 0.1)+ggtitle('P_ccle 5262 SNPs')
pl6=plot(xclinvar, option = "qqplot", threshold = 0.1)+ggtitle('P_clinvar hg38 17076 SNPs')
pl7=plot(xcpg, option = "qqplot", threshold = 0.1)+ggtitle('P_cpg hg38 2380 SNPs')
pl8=plot(xexome, option = "qqplot", threshold = 0.1)+ggtitle('P_exome 46156 SNPs')
pl9=plot(xfone, option = "qqplot", threshold = 0.1)+ggtitle('P_fone 736 SNPs')
pl10=plot(xgwas, option = "qqplot", threshold = 0.1)+ggtitle('P_gwas 525540 SNPs')
pl11=plot(xrand1_1000, option = "qqplot", threshold = 0.1)+ggtitle('P_rand1_1000 2177 SNPs')
pl12=plot(xrand1_150, option = "qqplot", threshold = 0.1)+ggtitle('P_rand1_150 326 SNPs')
pl13=plot(xrand1_300, option = "qqplot", threshold = 0.1)+ggtitle('P_rand1_300 718 SNPs')
pl14=plot(xrand1_50, option = "qqplot", threshold = 0.1)+ggtitle('P_rand1_50 108 SNPs')
pl15=plot(xtcga_somatic, option = "qqplot", threshold = 0.1)+ggtitle('P_tcga somatic 96284 SNPs')
pl16=plot(xtst, option = "qqplot", threshold = 0.1)+ggtitle('P_tst 365 SNPs')
grid.arrange(pl10,pl15,pl8,pl6,pl3,pl5,pl1,pl7,pl11,pl4,pl9,pl13,pl2,pl16,pl12,pl14, ncol=4)
dev.off()

png('snps_plot/stat_dist_all_3.png',width = 850,height = 850)
pl1=plot(x1321, option = "stat.distribution", threshold = 0.1)+ggtitle('P_1321 4052 SNPs')
pl2=plot(xaims, option = "stat.distribution", threshold = 0.1)+ggtitle('P_AIMS 445 SNPs')
pl3=plot(xbqccle, option = "stat.distribution", threshold = 0.1)+ggtitle('P_bq ccle mutation 5443 SNPs')
pl4=plot(xbroca, option = "stat.distribution", threshold = 0.1)+ggtitle('P_broca 2495 SNPs')
pl5=plot(xccle, option = "stat.distribution", threshold = 0.1)+ggtitle('P_ccle 5262 SNPs')
pl6=plot(xclinvar, option = "stat.distribution", threshold = 0.1)+ggtitle('P_clinvar 17076 SNPs')
pl7=plot(xcpg, option = "stat.distribution", threshold = 0.1)+ggtitle('P_cpg hg38 2380 SNPs')
pl8=plot(xexome, option = "stat.distribution", threshold = 0.1)+ggtitle('P_exome 46156 SNPs')
pl9=plot(xfone, option = "stat.distribution", threshold = 0.1)+ggtitle('P_fone 736 SNPs')
pl10=plot(xgwas, option = "stat.distribution", threshold = 0.1)+ggtitle('P_gwas 525540 SNPs')
pl11=plot(xrand1_1000, option = "stat.distribution", threshold = 0.1)+ggtitle('P_rand1_1000 2177 SNPs')
pl12=plot(xrand1_150, option = "stat.distribution", threshold = 0.1)+ggtitle('P_rand1_150 326 SNPs')
pl13=plot(xrand1_300, option = "stat.distribution", threshold = 0.1)+ggtitle('P_rand1_300 718 SNPs')
pl14=plot(xrand1_50, option = "stat.distribution", threshold = 0.1)+ggtitle('P_rand1_50 108 SNPs')
pl15=plot(xtcga_somatic, option = "stat.distribution", threshold = 0.1)+ggtitle('P_tcga somatic 96284 SNPs')
pl16=plot(xtst, option = "stat.distribution", threshold = 0.1)+ggtitle('P_tst 365 SNPs')
grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,pl9,pl10,pl11,pl12,pl13,pl14,pl15,pl16, ncol=4)
dev.off()


png('snps_plot/stat_hist_all_8.png',width = 950,height = 850)
par(mfrow=c(4,4))
bs2=list(xgwas,xtcga_somatic,xexome,xclinvar,xbqccle,xccle,x1321,xbroca,xcpg,xrand1_1000,xfone,xrand1_300,xaims,xtst,xrand1_150,xrand1_50)
for (i in 1:16){hist(bs2[[i]]$pvalues,col = 'black',xlab='p-values',breaks=40,main=paste(length(bs2[[i]]$pvalues), "SNPs in",pn1[i], sep=" "))}
dev.off()

png('snps_plot/stat_af_all_4.png',width = 950,height = 850)
par(mfrow=c(4,4))
for (i in 1:16){hist(bs2[[i]]$af,col = 'black',xlab='allele frequency',breaks=50,main=paste(length(bs2[[i]]$af), "SNPs in",pn1[i], sep=" "))}
dev.off()

png('snps_plot/stat_maf_all_4.png',width = 950,height = 850)
par(mfrow=c(4,4))
for (i in 1:16){hist(bs2[[i]]$maf,col = 'black',xlab='minor allele frequency',breaks=50,main=paste(length(bs2[[i]]$maf), "SNPs in",pn1[i], sep=" "))}
dev.off()

png('snps_plot/stat_stat_all_3.png',width = 950,height = 850)
par(mfrow=c(4,4))
for (i in 1:16){hist(bs2[[i]]$stat,col = rainbow(22),xlab='stat',breaks=50,main=paste(length(bs2[[i]]$stat), "SNPs in",pn1[i], sep=" "))}
dev.off()

png('snps_plot/stat_chi2.stat_3.png',width = 950,height = 850)
par(mfrow=c(4,4))
for (i in 1:16){hist(bs2[[i]]$chi2.stat,col = rainbow(22),xlab='chi2.stat',breaks=50,main=paste(length(bs2[[i]]$chi2.stat), "SNPs in",pn1[i], sep=" "))}
dev.off()

par(mfrow=c(4,4))
hist(x1321$pvalues, xlab = "p-values", main = 'P_1321 4052 SNPs', breaks = 50, col = "orange")
hist(xaims$pvalues, xlab = "p-values", main = 'P_AIMS 445 SNPs', breaks = 50, col = "orange")
hist(xbqccle$pvalues, xlab = "p-values", main = 'P_bq ccle mutation 5443 SNPs', breaks = 50, col = "orange")
hist(xbroca$pvalues, xlab = "p-values", main = 'P_broca 2495 SNPs', breaks = 50, col = "orange")
hist(xccle$pvalues, xlab = "p-values", main = 'P_ccle 5262 SNPs', breaks = 50, col = "orange")
hist(xclinvar$pvalues, xlab = "p-values", main = 'P_clinvar hg38 17076 SNPs', breaks = 50, col = "orange")
hist(xcpg$pvalues, xlab = "p-values", main = 'P_cpg hg38 2380 SNPs', breaks = 50, col = "orange")
hist(xexome$pvalues, xlab = "p-values", main = 'P_exome 46156 SNPs', breaks = 50, col = "orange")
hist(xfone$pvalues, xlab = "p-values", main = 'P_fone 736 SNPs', breaks = 50, col = "orange")
hist(xgwas$pvalues, xlab = "p-values", main = 'P_gwas 525540 SNPs', breaks = 50, col = "orange")
hist(xrand1_1000$pvalues, xlab = "p-values", main = 'P_rand1_1000 2177 SNPs', breaks = 50, col = "orange")
hist(xrand1_150$pvalues, xlab = "p-values", main = 'P_rand1_150 326 SNPs', breaks = 50, col = "orange")
hist(xrand1_300$pvalues, xlab = "p-values", main = 'P_rand1_300 718 SNPs', breaks = 50, col = "orange")
hist(xrand1_50$pvalues, xlab = "p-values", main = 'P_rand1_50 108 SNPs', breaks = 50, col = "orange")
hist(xtcga_somatic$pvalues, xlab = "p-values", main = 'P_tcga somatic 96284 SNPs', breaks = 50, col = "orange")
hist(xtst$pvalues, xlab = "p-values", main = 'P_tst 365 SNPs', breaks = 50, col = "orange")
#grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,pl9,pl10,pl11,pl12,pl13,pl14,pl15,pl16, ncol=4)
dev.off()

png('snps_plot/scree_all_1.png',width = 850,height = 850)
pl1=plot(x1321, option = "screeplot", threshold = 0.1)+ggtitle('P_1321 4052 SNPs')
pl2=plot(xaims, option = "screeplot", threshold = 0.1)+ggtitle('P_AIMS 445 SNPs')
pl3=plot(xbqccle, option = "screeplot", threshold = 0.1)+ggtitle('P_bq ccle mutation 5443 SNPs')
pl4=plot(xbroca, option = "screeplot", threshold = 0.1)+ggtitle('P_broca 2495 SNPs')
pl5=plot(xccle, option = "screeplot", threshold = 0.1)+ggtitle('P_ccle 5262 SNPs')
pl6=plot(xclinvar, option = "screeplot", threshold = 0.1)+ggtitle('P_clinvar 17076 SNPs')
pl7=plot(xcpg, option = "screeplot", threshold = 0.1)+ggtitle('P_cpg hg38 2380 SNPs')
pl8=plot(xexome, option = "screeplot", threshold = 0.1)+ggtitle('P_exome 46156 SNPs')
pl9=plot(xfone, option = "screeplot", threshold = 0.1)+ggtitle('P_fone 736 SNPs')
pl10=plot(xgwas, option = "screeplot", threshold = 0.1)+ggtitle('P_gwas 525540 SNPs')
pl11=plot(xrand1_1000, option = "screeplot", threshold = 0.1)+ggtitle('P_rand1_1000 2177 SNPs')
pl12=plot(xrand1_150, option = "screeplot", threshold = 0.1)+ggtitle('P_rand1_150 326 SNPs')
pl13=plot(xrand1_300, option = "screeplot", threshold = 0.1)+ggtitle('P_rand1_300 718 SNPs')
pl14=plot(xrand1_50, option = "screeplot", threshold = 0.1)+ggtitle('P_rand1_50 108 SNPs')
pl15=plot(xtcga_somatic, option = "screeplot", threshold = 0.1)+ggtitle('P_tcga somatic 96284 SNPs')
pl16=plot(xtst, option = "screeplot", threshold = 0.1)+ggtitle('P_tst 365 SNPs')
grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,pl9,pl10,pl11,pl12,pl13,pl14,pl15,pl16, ncol=4)
dev.off()



