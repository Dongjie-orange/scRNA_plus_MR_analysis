## 孟德尔后：验证集、双向孟德尔、共定位分析

library(tidyverse)

load('exposure_dat.Rdata')

gene=read.table('./OR.txt',sep = '\t',header = T)
gene=gene$exposure

exposure_dat=exposure_dat[exposure_dat$exposure %in% gene,]

### 验证集----

#提取结局数据,换一个队列
library(TwoSampleMR)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST003156")

# 取交集
exposure_dat=exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]


harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)



## MR

mr_modified <- function(dat = harmonised_dat, prop_var_explained = T)
{
  mr_res <- mr(dat)
  
  pve <- dat %>% 
    dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
    dplyr::group_by(id.exposure) %>% 
    dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
  
  if(prop_var_explained)
  {
    mr_res <- mr_res %>% 
      dplyr::left_join(pve, by = "id.exposure")
  }
  
  return(mr_res)
}




mr_res_vali <- mr_modified(harmonised_dat, prop_var_explained = T)

save(mr_res_vali,harmonised_dat,file ='mr_input_res_vali.Rdata')


result_or=generate_odds_ratios(mr_res_vali)
write.table(result_or[,4:ncol(result_or)],"OR_vali.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图


library(grid)
library(forestploter)


mydata=read.table("OR_vali.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =5 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)


## 其他可视化（一般不太好用，因为一个基因只有很有限的eqtl
mr_res1=mr_res_vali[mr_res_vali$exposure=='APOBEC3G',]
harmonised_dat1=harmonised_dat[harmonised_dat$exposure=='APOBEC3G',]

#绘制散点图
mr_scatter_plot(mr_res1, harmonised_dat1)

#森林图
res_single=mr_singlesnp(harmonised_dat1)
mr_forest_plot(res_single)

#漏斗图
mr_funnel_plot(singlesnp_results = res_single)

#留一法敏感性分析
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))



### 反向孟德尔
# 暴露
bimr_SLE<- extract_instruments('ebi-a-GCST90011866', p1=5e-08, clump=TRUE)

#提取结局数据
outcome_gene<- extract_outcome_data(snps=bimr_SLE$SNP, outcomes="eqtl-a-ENSG00000239713")

# 取交集
bimr_SLE =bimr_SLE [bimr_SLE $SNP %in% outcome_gene$SNP,]


harmonised_SLE_gene <- harmonise_data(bimr_SLE, outcome_gene)

bimr_mr_SLE_gene <- mr(harmonised_SLE_gene)

result_or=generate_odds_ratios(bimr_mr_SLE_gene)
write.table(result_or[,4:ncol(result_or)],"bi_OR.txt",row.names = F,sep = "\t",quote = F)

#将统计结果绘制森林图


library(grid)
library(forestploter)


mydata=read.table("bi_OR.txt",header = T,sep = "\t")
## !!
mydata$outcome='APOBEC3G'
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
forest(mydata[,c(1:3,6,12,13,14)],
       est = mydata$or,
       lower =mydata$or_lci95, 
       upper = mydata$or_uci95,
       sizes =0.3,
       ci_column =6 ,
       ref_line = 1,
       xlim = c(0.05, 3),
)



### 共定位分析----

load('mr_input_res.Rdata')

#如果表型是二分类变量（case和control），输入文件二选一：
#1）rs编号`rs_id`、P值`pval_nominal`、SNP的效应值`beta`、效应值方差`varbeta`；（推荐）
#2）rs编号`rs_id`、P值`pval_nominal`、case在所有样本中的比例`s`，MAF也要，写在list最后

#如果表型是连续型变量，输入文件三选一：
#1）rs编号`rs_id`、P值`pval_nominal`、表型的标准差`sdY`；
#2）rs编号`rs_id`、P值`pval_nominal`、效应值`beta`,效应值方差 `varbeta`, 样本量`N`,次等位基因频率 `MAF`；
#3）rs编号`rs_id`、P值`pval_nominal`、次等位基因频率 `MAF`；(推荐)


## 下载
data <- vcfR::read.vcfR("./eqtl-a-ENSG00000239713.vcf/eqtl-a-ENSG00000239713.vcf")

# 1.SNP ID rs开头
# 2.effect allele，此处相当于ALT
# 3.other allele，此处相当于REF
# 4.beta
# 5.se
# 6.pval


#整理数据
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


#	
#22
#39605007

eqtl=eqtl[eqtl$chr==22,]
eqtl$logpvalue=as.numeric(eqtl$logpvalue)
eqtl$p_value=10^(-eqtl$logpvalue)

eqtl$pos=as.numeric(eqtl$pos)

## 上下1mkb
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


##gwas:疾病
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

##---开始coloc
library(coloc)



input <- merge(my_eqtl, gwas, by="snp", all=FALSE, suffixes=c("_eqtl","_gwas"))


head(input)
library(coloc)

# 如报错，去除重复值
#input=input[!duplicated(input$snp),]

coloc.abf(dataset1=list(pvalues=input$pvalues_gwas, type="cc", s=0.33, N=12653),
                    dataset2=list(pvalues=input$pvalues_eqtl, type="quant", N=30684),
                    MAF=input$MAF)

## 或者用下面的方法

#my_eqtl=as.list(my_eqtl)
#my_eqtl[['type']]='quant'
#my_eqtl[['N']]=31430

#gwas=as.list(gwas)
#gwas[['type']]='cc'
#gwas[['N']]=12653
#gwas$pvalues=NULL
#coloc.abf(dataset1 = gwas,dataset2 = my_eqtl)


## 解读------
#H0：该区域的两个性状都没有遗传关联
#H1 / H2：只有表型1或表型2在该区域具有遗传关联
#H3：两个特征都相关，但因果变量不同
#H4：两个特征都相关并且共享一个因果变量

