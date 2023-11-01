## 孟德尔后：区域关联图、phenoscanner、方向过滤

## 区域关联图----------------------------

load('mr_input_res.Rdata')

data <- vcfR::read.vcfR("./eqtl-a-ENSG00000239713.vcf/eqtl-a-ENSG00000239713.vcf")

## 处理gt
## 处理gt
gt=data@gt
gt=as.data.frame(gt)

colnames(gt)
gt$FORMAT[1]

library(tidyverse)

gt$`eqtl-a-ENSG00000239713`[1]

##!!!
gt=separate(gt,col='eqtl-a-ENSG00000239713',into = c('ES', 'SE',
                                                     'LP','AF','SS',
                                                     'ID'),sep = '\\:')

gc()



gt=na.omit(gt)
colnames(gt)=c('format','beta','se','logpvalue','eaf','samplesize','snp')
gt$beta=as.numeric(gt$beta)
gt$se=as.numeric(gt$se)
gt$logpvalue=as.numeric(gt$logpvalue)
gt$eaf=as.numeric(gt$eaf)
gt$samplesize=as.numeric(gt$samplesize)

gc()
gt$format=NULL


fix=data@fix
fix=as.data.frame(fix)
colnames(fix)
colnames(fix)=c('chr','pos','snp','ref','alt')
fix=fix[,1:5]


## 合并gt fix
eqtl=left_join(fix,gt,by='snp')
eqtl=na.omit(eqtl)
## 查找染色体和位置
#22
#39605007

#MAFs代表次要等位基因的频率，因此其范围在0到0.5之间。值大于0.5的所有EAF通过从1中减去它们而转化为MAF。
eqtl$maf = ifelse(eqtl$eaf < 0.5, 
                  eqtl$eaf,
                  1 - eqtl$eaf)
eqtl$eaf=NULL


#22
#39605007
eqtl=eqtl[eqtl$chr==22,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)
eqtl=eqtl[eqtl$pos > 39605007-1000000 ,]
eqtl=eqtl[eqtl$pos < 39605007+1000000 ,]

my_eqtl=eqtl[,c('snp','p_value','maf')]

colnames(my_eqtl)=c('snp','pvalues','MAF')
my_eqtl=na.omit(my_eqtl)

my_eqtl=my_eqtl[my_eqtl$MAF>0 ,]

##-----接下来gwas疾病
library(TwoSampleMR)
## ！！！！！
coloc_SLE_dat <- extract_outcome_data(snps = c(my_eqtl$snp),
                                      outcomes = "ebi-a-GCST90011866",
                                      proxies = F) %>% 
  mutate(chr.outcome = as.numeric(chr),
         pos.outcome = as.numeric(pos),
         outcome = "Systemic lupus erythematosus", 
         id.outcome = "ebi-a-GCST90011866")

gwas=coloc_SLE_dat
gwas$beta=as.numeric(gwas$beta.outcome)
gwas$se=as.numeric(gwas$se.outcome)
#X=gwas$beta/gwas$se
#P=2*pnorm(q=abs(X), lower.tail=FALSE) #Z=1.96  P=0.05
#gwas$pvalue=P
## check Lp=log10(P)
gwas$varbeta=(gwas$se)^2

gwas=gwas[,c('SNP','pval.outcome',"beta",'varbeta')]

colnames(gwas)=c('snp','pvalues','beta','varbeta')

gwas=na.omit(gwas)

## 作图

library(locuscomparer)


## 以下都整理数据，请关注数据变化
snp_overlap <- intersect(my_eqtl[['snp']],
                         gwas[['snp']])
snp_overlap <- unique(snp_overlap)

exposure_dat <-my_eqtl[my_eqtl[['snp']] %in% snp_overlap,]
outcome_dat <- gwas[gwas[['snp']] %in% snp_overlap,]

exposure_dat <- exposure_dat[order(exposure_dat[['snp']]),]
outcome_dat <- outcome_dat[order(outcome_dat[['snp']]),]

exposure_dat <- data.frame(
  rsid = exposure_dat[['snp']],
  pval = exposure_dat[['pvalues']])

outcome_dat <- data.frame(
  rsid = outcome_dat[['snp']],
  pval = outcome_dat[['pvalues']])



write.table(exposure_dat,"coloc_exposure_test.tsv",sep = "\t",row.names = F,quote = F)
write.table(outcome_dat,"coloc_outcome_test.tsv",sep = "\t",row.names = F,quote = F)

library(locuscomparer)
p <- locuscompare(in_fn1 ="coloc_outcome_test.tsv", 
                  in_fn2 ="coloc_exposure_test.tsv", 
                  title1 = paste0('SLE'," GWAS"),  #！！！ 改改
                  title2 = paste0('APOBEC3G'," eQTL"), # !!!!改改
                  combine =T
)
p

## phenoscanner-------------------------

load('table1.Rdata')
## 固定的函数
mr_phenoscanner <- function(dat, 
                            catalog = c("GWAS", "eQTL", "pQTL", "mQTL", "methQTL"),
                            pvalue = 5e-08,
                            proxies = "EUR", 
                            r2 = 0.8,
                            build = 37)
{
  stopifnot("SNP" %in% names(dat))
  
  snpquery <- unique(dat$SNP)
  query_limit <- ceiling(0.01*length(snpquery))
  
  query_res <- list()
  for (i in 1:query_limit) {
    lo_i <- 100*i - 99
    up_i <- min(100*i,length(snpquery))
    snp_i <- snpquery[lo_i:up_i]
    print(paste0("Searching for ", snp_i))
    
    for (catalogue in catalog) {
      query_i <- try(
        phenoscanner::phenoscanner(
          snpquery = snp_i, catalogue = catalogue,
          proxies = proxies, pvalue = pvalue, r2 = r2, build = build)
      )
      if(!"try-error" %in% class(query_i))
      {
        if("results" %in% names(query_i))
        {
          query_i <- list(query_i$results)
          names(query_i) <- catalogue
          query_res <- append(query_res, query_i)
        }
      }
    }
  }
  
  return(query_res)
}


pleio_pheno <- mr_phenoscanner(table1) # the key information is the SNP column


pleio_pheno_summary_raw <- pleio_pheno$GWAS %>% 
  mutate(catalog = "GWAS") %>% 
  #bind_rows(pleio_pheno$pQTL%>% 
  #            mutate(catalog = "pQTL")) %>% 
  #bind_rows(pleio_pheno$eQTL%>% 
  #            mutate(catalog = "eQTL")) %>% 
  filter(p != "NA",
         ancestry %in% c("Mixed", "European")) %>% 
  mutate(across(p, as.numeric)) %>% 
  group_by(trait, pmid) %>% 
  filter(p == min(p)) %>% 
  ungroup()




# Steiger filtering  方向过滤----------------------------------------------
#比较SNP-exposure和SNP-outcome的R2
#如果exposure与outcome相比，在因果链上的位置离SNP更近
#则一般情况下前者的R2应大于后者

?steiger_filtering
load('mr_input_res.Rdata')
steiger_res= harmonised_dat %>%
  filter(SNP %in% table1$SNP) %>% 
  steiger_filtering()
steiger_res$steiger_pval

steiger_res2= harmonised_dat %>%
  filter(SNP %in% table1$SNP) %>% 
  directionality_test()


## 写出
write.csv(table1,file ='table1.csv',quote = F)

