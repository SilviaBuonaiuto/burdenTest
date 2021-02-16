library(tidyverse)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <casesCountsfile> <controlCountsfile> <genesFile> <numberCases> <numberControls>', call.=FALSE)
} 

####### open input files
cases<- read.table(args[1], sep = "\t", header = T)
controls<-read.table(args[2], sep = "\t", header = T )
genes<- read.table(args[3], sep = "\t", header = T) 
numberCases = as.numeric(args[4])
numberControls = as.numeric(args[5])

####Create input files
# INPUT: table og genes and counts
# ca1Het --> count of cases with at least one heterozygous qualifying variant in the gene
# co1Het --> count of controls with at least one heterozygous qualifying variant in the gene
# ca1Hom --> count of cases with at least one homozygous qualifying varinat in the gene
# co1Hom --> count of controls with at least one homozygous qualifying variant in the gene 

caseHet<-genes %>% unite(key, c("position", "allele"), sep = ":/") %>% merge(cases, by = "key") %>% select(SYMBOL, ID) %>% distinct() %>%  group_by(SYMBOL) %>% tally() %>% rename(ca1Het = n)
caseHom<-genes %>% unite(key, c("position", "allele"), sep = ":/") %>% merge(cases, by = "key") %>% filter(ALTcount == 2) %>% select(SYMBOL, ID) %>% distinct() %>%  group_by(SYMBOL) %>% tally() %>% rename(ca1Hom = n)

controlLHet<-genes %>% unite(key, c("position", "allele"), sep = ":/") %>% merge(controls, by = "key") %>% select(SYMBOL, ID) %>% distinct() %>%  group_by(SYMBOL) %>% tally() %>% rename(co1Het = n)
controlLHom<-genes %>% unite(key, c("position", "allele"), sep = ":/") %>% merge(controls, by = "key") %>% filter(ALTcount == 2) %>% select(SYMBOL, ID) %>% distinct() %>%  group_by(SYMBOL) %>% tally() %>% rename(co1Hom = n)

caseAll<-merge(caseHet, caseHom, by = "SYMBOL", all = T) %>% replace(is.na(.), 0)
controlAll<-merge(controlLHet, controlLHom, by = "SYMBOL", all = T) %>% replace(is.na(.), 0)
myd<- merge(controlAll, caseAll, by = "SYMBOL")


######### pvalues recessive and dominant
myd<- myd %>%  rowwise() %>% mutate(cpdom=chisq.test(matrix(c(ca1Het, numberCases-ca1Het, co1Het, numberControls-co1Het), ncol=2, nrow=2))$p.value, cstatdom=chisq.test(matrix(c(ca1Het, numberCases-ca1Het, co1Het, numberControls-co1Het), ncol=2, nrow=2))$statistic, cprec=chisq.test(matrix(c(ca1Hom, numberCases-ca1Hom, co1Hom, numberControls-co1Hom), ncol=2, nrow=2))$p.value, cstatrec=chisq.test(matrix(c(ca1Hom, numberCases-ca1Hom, co1Hom, numberControls-co1Hom), ncol=2, nrow=2))$statistic, log10cpdom=-log10(cpdom),  log10cprec=-log10(cprec) )


######## calulate expected pval under chisquare distribution with 1 d.f.
numberGenes = nrow(myd)
samplingXsquare<-rchisq(numberGenes,df=1)
expPX<-pchisq(samplingXsquare, df=1)
myd<-myd[order(myd$log10cpdom),]
myd$log10cpdomexp<-sort(-log10(expPX))
myd<-myd[order(myd$log10cprec),]
myd$log10cprecexp<-sort(-log10(expPX))

myd %>% write.table("burden_table.tsv", sep = "\t", quote = F, col.names =T, row.names=F)
###plot
pdominant<-ggplot(myd,aes(log10cpdomexp, log10cpdom)) + geom_point() + geom_text_repel(data = subset(myd,log10cpdom >=4), aes(label = SYMBOL)) + geom_abline(intercept = 0, slope = 1)
ggsave("rareVariants_dominant.png", plot = pdominant, width = 6, height = 6 )

precessive<-ggplot(myd,aes(log10cprecexp, log10cprec)) + geom_point() + geom_abline(intercept = 0, slope = 1)
ggsave("rareVariants_recessive.png", plot = precessive, width = 6, height = 6 )
########  calculate lambda
numberGenes = nrow(myd)
samplingXsquare<-rchisq(numberGenes,df=1)
lambdaX=median(myd$cstatdom, na.rm=T)/median(samplingXsquare)

#######correct
expPX<-pchisq(samplingXsquare, df=1)

df <- data.frame(pcorrect = myd$cpdom/lambdaX)
df$logpcorrect <- -log10(df$pcorrect)
df<-df[order(df$logpcorrect),]
df$log10pexp<-sort(-log10(expPX))

####plot after correction
correctDominant<-ggplot(df,aes(log10pexp, logpcorrect)) + geom_point()  + geom_abline(intercept = 0, slope = 1)
ggsave("rareVariants_dominant_correct.png", plot = correctDominant, width = 6, height = 6 )
