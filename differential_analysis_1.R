library(biomaRt) #Load required libraries
library("BiocParallel")
register(MulticoreParam(40)) #Use all the available cores, its a large genome
library(edgeR)
library(dplyr)
library(plyr)
library(pheatmap)
library(biomaRt)
library(tidyverse)
library(forcats)
library(ggfortify)
library(ClassDiscovery)
library(Pigengene)
setwd('/shire/zillur/rna_seq/again/')
bams
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") #Load and get ids
ids <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','entrezgene', 'refseq_mrna','interpro','interpro_description'),mart = ensembl)
ftc1=featureCounts(bams,nthreads = 40,annot.inbuilt = 'hg38',isPairedEnd = T)
ftc3=featureCounts(bams3,nthreads = 40,annot.ext = '/shire/zillur/rna_seq/grch38/gencode.v31.primary_assembly.annotation.gtf',isGTFAnnotationFile = T,isPairedEnd = T,useMetaFeatures = T,GTF.attrType = 'gene_name')

cnt1=ftc1$counts
colnames(cnt1)=substring(colnames(cnt1),1,5)
colnames(cnt1)[13:51]=substring(colnames(cnt1)[13:51],1,4)
cnt2=cnt1[apply(cnt1!=0,1,all),]
cnt3=cnt2[,2:52] #Low mapping quality
autoplot(prcomp(t(cnt3)),label=T)
cnt4=subset(cnt3,select = -c(38,39)) #Mahalnabish distance
id2=ids
id2$entrezgene_id=as.character(id2$entrezgene_id)
cnt4$hgnc_symbol=id2$hgnc_symbol[match(rownames(cnt4),id2$entrezgene_id)]
cnt5=cnt4[complete.cases(cnt4),]
cnt5$ensembl_gene_id=id2$ensembl_gene_id[match(cnt5$hgnc_symbol,id2$hgnc_symbol)]
pd1=DGEList(counts = cnt5[,1:48],genes = cnt5[,49:51])
pd1$samples$group=c(rep('Control',14),rep('Case',34))
pd2=calcNormFactors(pd1)
pd3=estimateDisp(pd2)
summary(decideTests(etcc1))
pd4=pd3
pd4$samples$drc=drc1$state[match(rownames(pd4$samples),drc1$ID)]
pd4$samples$group=pd4$samples$drc
cc1=model.matrix(~0+pd3$samples$group)
colnames(cc1)=c('Case','Control')
clsg=estimateDisp(pd3,design = cc1)
ftcc=glmQLFit(clsg,design = cc1)
qlcc=glmQLFTest(ftcc,contrast = c(1,-1)) #DE glm approach
summary(decideTests(qlcc))

dc1=model.matrix(~0+pd4$samples$group)
colnames(dc1)=c('High','Low')
dlsg=estimateDisp(pd4,design = dc1)
dtcc=glmQLFit(dlsg,design = dc1)
dlcc=glmQLFTest(dtcc,contrast = c(1,-1)) #DE glm approach
summary(decideTests(dlcc))

kdrc=kegga(dlcc,geneid = dlcc$genes$entrezgene, species='Hs',N=nrow(dlcc))
drcup=topKEGG(kdrc,sort = 'Up',number = 333)
drcup2=drcup[drcup$P.Up<=0.05,]
drcdown=topKEGG(kdrc,sort = 'Down',number = 333)
drcdown2=drcdown[drcdown$P.Down<=0.05,]


#New from ftc3
ntc1=ftc3$counts
write.csv(as.data.frame(ntc1),'/shire/zillur/rna_seq/abstract/count_fin_1.csv')
colnames(ntc1)=gsub('Aligned.sortedByCoord.out.bam','',colnames(ntc1))
colnames(ntc1)[1:14]=gsub('.{4}$','',colnames(ntc1)[1:14])
ntc2$entrezgene=ids$entrezgene_id[match(rownames(ntc2),ids$hgnc_symbol)]
summary(ntc2$entrezgene)
ntc3=ntc2[complete.cases(ntc2),]

autoplot(prcomp(t(ntc3[,1:50])),label=T,shape=F)
sp1=SamplePCA(ntc3[,1:50])
round(cumsum(sp1@variances)/sum(sp1@variances),digits = 2)
m1=mahalanobisQC(sp1,2)
m1=m1[order(m1$p.value),]
write.csv(m1,'/shire/zillur/rna_seq/abstract/mahalnabish_ordered.csv')
dg1=DGEList(counts = ntc3[,1:50],genes = ntc3[,51:52])
dg1$samples$group=c(rep('Control',14),rep('Case',36))
keep <- rowSums(cpm(dg1)>=1) >= 1
dg1=dg1[keep, , keep.lib.size=F]
dg2=calcNormFactors(dg1)
plotMDS(dg2)
des1=model.matrix(~0+dg2$samples$group)
colnames(des1)=c('Case','Control')
dg2=estimateDisp(dg2,design = des1)
plotBCV(dg2)
qlft1=glmQLFit(dg2,design = des1)
plotQLDisp(qlft1)
con1=makeContrasts(Case - Control, levels = des1)
qlf1=glmQLFTest(qlft1,contrast = con1)
t1=topTags(qlf1,n=20000)
summary(decideTests(qlf1))


dg4=dg3
dg4$samples$drc=drc1$state[match(rownames(dg4$samples),drc1$ID)]
dg4$samples$group=dg4$samples$drc
dg4=calcNormFactors(dg4)
plotMDS(dg4)
des2=model.matrix(~0+dg4$samples$group)
colnames(des2)=c('High','Low')
dg4=estimateDisp(dg4,design = des2)
plotBCV(dg4)
qlft2=glmQLFit(dg4,design = des2)
plotQLDisp(qlft2)
con2=makeContrasts(High - Low, levels = des2)
qlf2=glmQLFTest(qlft2,contrast = con2)
t2=topTags(qlf2,n=20000)
summary(decideTests(qlf2))
g1=goana(qlf2,geneid = qlf2$genes$entrezgene)
g1_up=topGO(g1,ontology = 'BP',sort = 'up')
g1_down=topGO(g1,ontology = 'BP',sort = 'down')
k1=kegga(qlf2,geneid = qlf2$genes$entrezgene)
k1_up=topKEGG(k1,sort = 'up')
k1_down=topKEGG(k1,sort = 'down')
write.csv(g1_up,'/shire/zillur/rna_seq/abstract/upregulated_processes_drc_high_1.csv')
write.csv(g1_down,'/shire/zillur/rna_seq/abstract/downregulated_processes_drc_high_1.csv')
write.csv(k1_up,'/shire/zillur/rna_seq/abstract/upregulated_pathways_drc_high_1.csv')
write.csv(k1_down,'/shire/zillur/rna_seq/abstract/downregulated_pathways_drc_high_1.csv')
dr1=t2$table[order(t2$table$PValue),]
dr2=dr1[dr1$PValue<0.05,]
write.csv(dr2,'/shire/zillur/rna_seq/abstract/genes_drc_high_1.csv')

dg5=dg3
dg5$samples$subtype=sb1$Subtype[match(rownames(dg5$samples),sb1$ID_.)]
dg5$samples$subtype=as.character(dg5$samples$subtype)
dg5$samples$subtype[1:14]=rep('Control',14)
sb2=dg5$samples[complete.cases(dg5$samples),]
dg6=DGEList(counts = ntc3[,rownames(sb2)],genes = ntc3[51:52])
dg6$samples$subtype=sb2$subtype
dg6$samples$group=dg6$samples$subtype
keep <- rowSums(cpm(dg6)>=1) >= 1
dg6=dg6[keep, , keep.lib.size=F]
dg6=calcNormFactors(dg6)
plotMDS(dg6)
des3=model.matrix(~0+dg6$samples$group)
colnames(des3)=c('Control','HER2','LumA','LumB','TN')
dg6=estimateDisp(dg6,design = des3)
plotBCV(dg6)
qlft3=glmQLFit(dg6,design = des3)
plotQLDisp(qlft3)
qlf_her2_cont=glmQLFTest(qlft3, contrast = c(-1,1,0,0,0))
t3=topTags(qlf_her2_cont,n=15000)
summary(decideTests(qlf_her2_cont))
qlf_luma_cont=glmQLFTest(qlft3, contrast = c(-1,0,1,0,0))
t4=topTags(qlf_luma_cont,n=15000)
summary(decideTests(qlf_luma_cont))
qlf_lumb_cont=glmQLFTest(qlft3, contrast = c(-1,0,0,1,0))
t5=topTags(qlf_lumb_cont,n=15000)
summary(decideTests(qlf_lumb_cont))
qlf_tn_cont=glmQLFTest(qlft3, contrast = c(-1,0,0,0,1))
t6=topTags(qlf_tn_cont,n=15000)
summary(decideTests(qlf_tn_cont))

anc1=read.csv('/shire/zillur/rna_seq/ancestry_quartile.csv')
dg5$samples$afrq=anc1$afrq[match(rownames(dg5$samples),anc1$Sample_ID)]
dg5$samples$eurq=anc1$eurq[match(rownames(dg5$samples),anc1$Sample_ID)]
dg5$samples$natq=anc1$natq[match(rownames(dg5$samples),anc1$Sample_ID)]
anc2=dg5$samples[complete.cases(dg5$samples),]
dg7=DGEList(counts = ntc3[,rownames(anc2)],genes = ntc3[,51:52])
dg7$samples$afrq=fct_collapse(dg7$samples$afrq, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))
dg7$samples$eurq=fct_collapse(dg7$samples$eurq, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))
dg7$samples$natq=fct_collapse(dg7$samples$natq, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))

des4=model.matrix(~0+dg7$samples$afrq)
keep <- rowSums(cpm(dg7)>=1) >= 1
dg7=dg7[keep, , keep.lib.size=F]
dg7=calcNormFactors(dg7)
plotMDS(dg7)
dg7=estimateDisp(dg7,design = des4)
plotBCV(dg7)
qlft4=glmQLFit(dg7,design = des4)
plotQLDisp(qlft4)
qlf_afrq3_afrq1=glmQLFTest(qlft4, contrast = c(-1,1))
t7=topTags(qlf_afrq3_afrq1,n=15000)
summary(decideTests(qlf_afrq3_afrq1))

des5=model.matrix(~0+dg7$samples$eurq)
qlft5=glmQLFit(dg7,design = des5)
qlf_eurq3_eurq1=glmQLFTest(qlft5, contrast = c(-1,1))
t8=topTags(qlf_eurq3_eurq1,n=15000)
summary(decideTests(qlf_eurq3_eurq1))

des6=model.matrix(~0+dg7$samples$natq)
qlft6=glmQLFit(dg7,design = des6)
qlf_natq3_natq1=glmQLFTest(qlft6, contrast = c(-1,1))
t9=topTags(qlf_natq3_natq1,n=15000)
summary(decideTests(qlf_natq3_natq1))

dn1=read.csv('/shire/zillur/mac_prev/mzillur/rna_analysis/de/initial_results_2/dna_reapir_annotated.csv')
dn2=dn1[,1:3]
dn2$Case_vs_control=t1$table$logFC[match(dn2$gene,rownames(t1$table))]
dn2$DRC_High_vs_Low=t2$table$logFC[match(dn2$gene,rownames(t2$table))]
dn2$HER2_vs_Control=t3$table$logFC[match(dn2$gene,rownames(t3$table))]
dn2$LumA_vs_Control=t4$table$logFC[match(dn2$gene,rownames(t4$table))]
dn2$LumB_vs_Control=t5$table$logFC[match(dn2$gene,rownames(t5$table))]
dn2$TN_vs_Control=t6$table$logFC[match(dn2$gene,rownames(t6$table))]
dn2$AFR_High_vs_Low=t7$table$logFC[match(dn2$gene,rownames(t7$table))]
dn2$EUR_High_vs_Low=t8$table$logFC[match(dn2$gene,rownames(t8$table))]
dn2$NAT_High_vs_Low=t9$table$logFC[match(dn2$gene,rownames(t9$table))]

dn3=dn2[complete.cases(dn2),]
rownames(dn3)=dn3$gene
write.csv(dn3,'/shire/zillur/rna_seq/abstract/dna_repair_logfc_1.csv')
dn4=dn3[!rownames(dn3) %in% c('GTF2H4','NHEJ1'),]
png('/shire/zillur/rna_seq/abstract/dna_repair_heat_4.png', height = 5000,width = 5000)
annotdf1=data.frame(row.names = rownames(dn3),pathway=dn3$pathway2)
newCols1 <- colorRampPalette(grDevices::rainbow(length(unique(annotdf1$pathway))))
mycolors1 <- newCols1(length(unique(annotdf1$pathway)))
names(mycolors1) <- unique(annotdf1$pathway)
mycolors1 <- list(pathway = mycolors1)
rowgap1=(head(as.numeric(cumsum(table(dn3$pathway2))), -1))
mycol2=c('yellow','red','forestgreen','chartreuse','brown1','dodgerblue','darkmagenta')
pheatmap.type(dn3[,(4:12)],color = mycol2, fontsize_row = 28,fontsize_col = 42,annRow = annotdf1,fontsize = 60,annotation_colors = mycolors1,gaps_row=rowgap1)
dev.off()


#pcaplots
dn_cp1=as.data.frame(cpm(dn_cnt1[,1:50],normalized.lib.sizes=T,log=T,gene.length=nrow(dn_cnt1)))
dn_cp1$pathway=dn3$pathway2[match(rownames(dn3),rownames(dn_cp1))]
dn_cp2=as.data.frame(t(dn_cp1))
dn_cp2$condition=dg3$samples$group[match(rownames(dg3$samples),rownames(dn_cp2))]
dn_cp2$drc=dg4$samples$drc[match(rownames(dg4$samples),rownames(dn_cp2))]
dn_cp3=dn_cp2
dn_cp4=dn_cp3[match(rownames(dg6$samples),rownames(dn_cp3)),]
dn_cp4$subtype=dg6$samples$subtype[match(rownames(dg6$samples),rownames(dn_cp4))]
dn_cp5=dn_cp3[match(rownames(dg7$samples),rownames(dn_cp3)),]
dn_cp5$afr=dg7$samples$afrq[match(rownames(dg7$samples),rownames(dn_cp5))]
dn_cp5$eur=dg7$samples$eurq[match(rownames(dg7$samples),rownames(dn_cp5))]
dn_cp5$nat=dg7$samples$natq[match(rownames(dg7$samples),rownames(dn_cp5))]

png('/shire/zillur/rna_seq/abstract/pca_dna_repair_166genes.png')
p1=autoplot(prcomp(dn_cp2[,1:166]),data = dn_cp2,colour='condition')
p2=autoplot(prcomp(dn_cp2[,1:166]),data = dn_cp2,colour='drc')
p3=autoplot(prcomp(dn_cp4[,1:166]),data = dn_cp4,colour='subtype')
p4=autoplot(prcomp(dn_cp5[,1:166]),data = dn_cp5,colour='afr')
p5=autoplot(prcomp(dn_cp5[,1:166]),data = dn_cp5,colour='eur')
p6=autoplot(prcomp(dn_cp5[,1:166]),data = dn_cp5,colour='nat')
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
dev.off()

d_cp1=as.data.frame(cpm(ntc3[,1:50],normalized.lib.sizes=T,log=T,gene.length=nrow(ntc3)))
d_cp2=as.data.frame(t(d_cp1))
d_cp2$condition=dg3$samples$group[match(rownames(dg3$samples),rownames(d_cp2))]
d_cp2$drc=dg4$samples$drc[match(rownames(dg4$samples),rownames(d_cp2))]
d_cp3=d_cp2
d_cp4=d_cp3[match(rownames(dg6$samples),rownames(d_cp3)),]
d_cp4$subtype=dg6$samples$subtype[match(rownames(dg6$samples),rownames(d_cp4))]
d_cp5=d_cp3[match(rownames(dg7$samples),rownames(d_cp3)),]
d_cp5$afr=dg7$samples$afrq[match(rownames(dg7$samples),rownames(d_cp5))]
d_cp5$eur=dg7$samples$eurq[match(rownames(dg7$samples),rownames(d_cp5))]
d_cp5$nat=dg7$samples$natq[match(rownames(dg7$samples),rownames(d_cp5))]

hk1=d_cp2 %>% select(one_of(dput(as.character(rownames(hkg1)))))
hk2=hk1
hk2$condition=dg3$samples$group[match(rownames(dg3$samples),rownames(hk2))]
hk2$drc=dg4$samples$drc[match(rownames(dg4$samples),rownames(hk2))]
hk3=hk2
hk4=hk3[match(rownames(dg6$samples),rownames(hk3)),]
hk4$subtype=dg6$samples$subtype[match(rownames(dg6$samples),rownames(hk4))]
hk5=hk3[match(rownames(dg7$samples),rownames(hk3)),]
hk5$afr=dg7$samples$afrq[match(rownames(dg7$samples),rownames(hk5))]
hk5$eur=dg7$samples$eurq[match(rownames(dg7$samples),rownames(hk5))]
hk5$nat=dg7$samples$natq[match(rownames(dg7$samples),rownames(hk5))]

png('/shire/zillur/rna_seq/abstract/pca_housekeeping_3732genes.png')
p1=autoplot(prcomp(hk2[,1:3732]),data = hk2,colour='condition')
p2=autoplot(prcomp(hk2[,1:3732]),data = hk2,colour='drc')
p3=autoplot(prcomp(hk4[,1:3732]),data = hk4,colour='subtype')
p4=autoplot(prcomp(hk5[,1:3732]),data = hk5,colour='afr')
p5=autoplot(prcomp(hk5[,1:3732]),data = hk5,colour='eur')
p6=autoplot(prcomp(hk5[,1:3732]),data = hk5,colour='nat')
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
dev.off()

tsg1=read.csv('/shire/zillur/rna_seq/Human_TSGs.txt',sep = '\t')
ts1=d_cp2 %>% select(one_of(dput(as.character(tsg1$GeneSymbol))))
ts2=ts1
ts2$condition=dg3$samples$group[match(rownames(dg3$samples),rownames(ts2))]
ts2$drc=dg4$samples$drc[match(rownames(dg4$samples),rownames(ts2))]
ts3=ts2
ts4=ts3[match(rownames(dg6$samples),rownames(ts3)),]
ts4$subtype=dg6$samples$subtype[match(rownames(dg6$samples),rownames(ts4))]
ts5=ts3[match(rownames(dg7$samples),rownames(ts3)),]
ts5$afr=dg7$samples$afrq[match(rownames(dg7$samples),rownames(ts5))]
ts5$eur=dg7$samples$eurq[match(rownames(dg7$samples),rownames(ts5))]
ts5$nat=dg7$samples$natq[match(rownames(dg7$samples),rownames(ts5))]

png('/shire/zillur/rna_seq/abstract/pca_tumor_suppressor_1048genes.png')
p1=autoplot(prcomp(ts2[,1:1048]),data = ts2,colour='condition')
p2=autoplot(prcomp(ts2[,1:1048]),data = ts2,colour='drc')
p3=autoplot(prcomp(ts4[,1:1048]),data = ts4,colour='subtype')
p4=autoplot(prcomp(ts5[,1:1048]),data = ts5,colour='afr')
p5=autoplot(prcomp(ts5[,1:1048]),data = ts5,colour='eur')
p6=autoplot(prcomp(ts5[,1:1048]),data = ts5,colour='nat')
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
dev.off()

og1=readxl::read_xlsx('/shire/zillur/rna_seq/NIHMS540676-supplement-02.xlsx')
og1=og1[-c(1,2),]
colnames(og1)=c(1,2,3,4)

og2=d_cp2 %>% select(one_of(dput(as.character(og1$1))))

og2$condition=dg3$samples$group[match(rownames(dg3$samples),rownames(og2))]
og2$drc=dg4$samples$drc[match(rownames(dg4$samples),rownames(og2))]
og3=og2
og4=og3[match(rownames(dg6$samples),rownames(og3)),]
og4$subtype=dg6$samples$subtype[match(rownames(dg6$samples),rownames(og4))]
og5=og3[match(rownames(dg7$samples),rownames(og3)),]
og5$afr=dg7$samples$afrq[match(rownames(dg7$samples),rownames(og5))]
og5$eur=dg7$samples$eurq[match(rownames(dg7$samples),rownames(og5))]
og5$nat=dg7$samples$natq[match(rownames(dg7$samples),rownames(og5))]


png('/shire/zillur/rna_seq/abstract/pca_onco_49genes.png')
p1=autoplot(prcomp(og2[,1:49]),data = og2,colour='condition')
p2=autoplot(prcomp(og2[,1:49]),data = og2,colour='drc')
p3=autoplot(prcomp(og4[,1:49]),data = og4,colour='subtype')
p4=autoplot(prcomp(og5[,1:49]),data = og5,colour='afr')
p5=autoplot(prcomp(og5[,1:49]),data = og5,colour='eur')
p6=autoplot(prcomp(og5[,1:49]),data = og5,colour='nat')
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
dev.off()


cpn1$group[rownames(cpn1) %in% rownames(ogn1)] <- "Oncogene"
cpn1$group[rownames(cpn1) %in% rownames(dn3)] <- "DNA Repair System"
cpn1$group[rownames(cpn1) %in% rownames(hk3)] <- "Housekeeping genes"
cpn1$group[rownames(cpn1) %in% rownames(ts6)] <- "Tumor Suppressor"

tsn5=Rtsne(cp5[,1:10978],check_duplicates=F,perplexity=1)
tsne_plot3 <- data.frame(x = tsn5$Y[,1], y = tsn5$Y[,2], col = cp5$group)
ggplot(tsne_plot3) + geom_point(aes(x=x, y=y, color=col))
tsne_plot4 <- data.frame(x = tsn5$Y[,1], y = tsn5$Y[,2], col = cp5$drc)
ggplot(tsne_plot4) + geom_point(aes(x=x, y=y, color=col))
tsne_plot5 <- data.frame(x = tsn5$Y[,1], y = tsn5$Y[,2], col = cp5$subtype)
ggplot(tsne_plot5) + geom_point(aes(x=x, y=y, color=col))




