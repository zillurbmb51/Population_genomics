BiocManager::install(c('RnBeads','RnBeads.hg38','DelayedMatrixStats','minfi','methylumi'))
library(RnBeads)
data.dir='/shire/zillur/cell_lines/ega_data/'
idat.dir=file.path(data.dir)
analysis.dir='/shire/zillur/cell_lines/analysis/'
report.dir=file.path(analysis.dir,'reports')
rnb.options(filtering.sex.chromosomes.removal=TRUE,identifiers.column="Sample_ID")
sample.annotation=read.table('/shire/zillur/cell_lines/ega_data/E-MTAB-2706.sdrf.txt',header = T,sep = '\t',fill = T)

rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation,data.dir=idat.dir, data.type="infinium.idat.dir")
idtfiles=read.idat.files(data.dir, barcodes = NULL, sample.sheet = sample.annotation,sep.samples = rnb.getOption("import.table.separator"), 
                         useff = FALSE,
                         verbose = TRUE)
f2=read.csv('/shire/zillur/thesis/verification/msmc/likelihood_all_5.csv',header = F,sep = '\t')
png('/shire/zillur/thesis/verification/msmc/likehoods_and_effective_population_size.png',width = 888)
par(mfrow=c(1,2))
plot(1:150,f2$V1,col='red',type='l',lwd=2,xlab='Iterations',ylab='log likelihood',main='Iterations vs likehoods')
lines(1:60,f2$V2[1:60],col='blue',lwd=2)
lines(1:55,f2$V3[1:55],col='green',lwd=2)
lines(1:50,f2$V4[1:50],col='black',lwd=2)
lines(1:45,f2$V5[1:45],col='yellow',lwd=2)
abline(v=55,lwd=2)
f1=read.table('/shire/zillur/thesis/verification/msmc/res_fin_1.final.txt',header = T)
f1=f1[-1,]
mu=4.6e-09
g=1/6
plot(f1$left_time_boundary/mu*g, (1/f1$lambda)/(2*mu),log='x', type = 'n', xlab='Years ago',ylab = 'Effective population size',main='Years vs effective population size')
lines(f1$left_time_boundary/mu*g,(1/f1$lambda)/(2*mu),type='s',col='red')
dev.off()
library(devtools)
BiocManager::install('crlmm')
library(crlmm)

fls1=readIdatFiles(sampleSheet=sample.annotation, arrayNames=arnam, ids=NULL, path="/shire/zillur/cell_lines/ega_data/",
              arrayInfoColNames=list(barcode="SentrixBarcode_A",
                                     position="SentrixPosition_A"),
              highDensity=FALSE, sep="_",
              fileExt=list(green="Grn.idat", red="Red.idat"),
              saveDate=FALSE, verbose=FALSE)
arnam=sample.annotation$Source.Name
smsht1=read.csv('/shire/zillur/cell_lines/ega_data/E-MTAB-2706.sdrf.txt',sep = '\t')
rdidat1=readIdatFiles(sampleSheet=smsht1, arrayNames=NULL, ids=NULL, path='/shire/zillur/cell_lines/ega_data/',
              arrayInfoColNames=list(barcode="SentrixBarcode_A",
                                     position="SentrixPosition_A"),
              highDensity=FALSE, sep="_",
              fileExt=list(green=c("Grn.idat","grn.idat"), red=c("Red.idat","red.idat")),
              saveDate=FALSE, verbose=FALSE)
