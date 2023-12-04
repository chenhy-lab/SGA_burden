# download vcf file and transform the data to the file: demo_burden_data_SGA.txt

source('./Function.R')
panel_gene <- read.delim("data/gene.lst")
SGA627 <- rio::import('data/SGA627.xlsx')
NonSGA1317 <- rio::import('data/NonSGA1317.xlsx')
sample_file <- rio::import('data/demo_burden_data_SGA.txt')

# SGA BurdenTest
case_size <- nrow(SGA627)
ctr_size <- nrow(NonSGA1317)
sample_file <- .myEssentialQc(sample_file)
collapsing_res <- collapsing_burden(sample_file)
collapsing_res <- collapsing_res %>% mutate(risk_score=2*(-log10(Fisher_PTVp))+(-log10(Fisher_MISp)))
collapsing_res_sig <- collapsing_res %>% dplyr::filter(Fisher_PTVp<0.05 | Fisher_MISp<0.05,PTV_OR>1 | MIS_OR>1 | PTV_OR=="Inf" | MIS_OR=="Inf",Fisher_SYNp>0.05,Fisher_NONp>0.05)

# SGA BurdenTest by collapsing analyses for expected p
if(T){
  mutation_dt <- sample_file[,c("Dis_gene","Mutation_type_(VEP)","Mutation_type_(annovar)","uniq_id")]
  group_list <- list()
  all_sample_name <- unique(sample_file$uniq_id)
  
  for (i in 1:10000) {
    set.seed(i)
    tmp_case_sample_name <- sample(all_sample_name,case_size)
    tmp_control_sample_name <- setdiff(all_sample_name,tmp_case_sample_name)
    tmp_sample_file <- within(sample_file,{
      diagnosis <- NA
      diagnosis[sample_file$uniq_id %in% tmp_case_sample_name] <- 'case'
      diagnosis[sample_file$uniq_id %in% tmp_control_sample_name] <- 'control'
    })
    group_list[[i]] <- tmp_sample_file[,c("diagnosis"),drop=F]
  }
}

## expected pvalue
if(T){
  library(snowfall)
  sfInit(parallel = TRUE, cpus = 10, slaveOutfile = "./snowfall_log.txt")  #initialize
  sfLibrary(tidyverse)
  sfExport("mutation_dt","group_list","collapsing_burden")  # load required object.
  collapsing_res_list <- snowfall::sfLapply(group_list,function(x){
    sample_file <- cbind(mutation_dt,x)
    collapsing_res <- collapsing_burden(sample_file=sample_file)
    return(collapsing_res)
  })
  sfStop()
  save(collapsing_res_list,file='data/exp_collapsing_res_list_SGA.rda')
}

if(T){
  load('data/exp_collapsing_res_list_SGA.rda')
  burdenscore_list <- lapply(collapsing_res_list, function(collapsing_res){
    risk_score=2*(-log10(collapsing_res$Fisher_PTVp))+(-log10(collapsing_res$Fisher_MISp))
    return(risk_score)
  })
  burdenscore <- bind_cols(burdenscore_list) %>% as.data.frame()
  colnames(burdenscore) <- 1:10000
  rownames(burdenscore) <- collapsing_res_list[[1]]$gene
  save(burdenscore,file = 'data/burdenscore_SGA_10000.rda')
}

load('data/burdenscore_SGA_10000.rda')
burdenscore <- as.matrix(burdenscore) %>% .[rowSums(.)!=0,]
collapsing_res <- collapsing_res %>% dplyr::filter(gene %in% rownames(burdenscore))
rownames(collapsing_res) <- collapsing_res$gene
collapsing_res$bootstrap_P <- unlist(lapply(collapsing_res$gene,function(x){
  length(which(collapsing_res[x,'risk_score']<burdenscore[x,]))/10000
}))
collapsing_res$fdr <- p.adjust(collapsing_res$bootstrap_P,method ="fdr")
collapsing_res

