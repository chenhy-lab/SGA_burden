library(tidyverse)
library(plyr)
library(MatchIt)
library(glmnet)
library(tidyverse)
library(pROC)
library(caret)
library(DALEX)
library(precrec)

# remove the column with one levels
.myrm_cons_var <- function(df){
  rm_var <- c();two_levels <<- c()
  for (i in colnames(df)) {
    if(length(levels(df[,i]))==1){rm_var <- c(rm_var,i)}
    if(length(table(df[,i]) %>% names()) == 1){rm_var <- c(rm_var,i)}
    if(length(levels(df[,i]))==2){two_levels <<- c(two_levels,i)}
  }
  df <- df[,setdiff(colnames(df),rm_var)]
  return(df)
}

##' Merge all data in sample list
#'
#' \code{merge_file} is a function to merge all dataframe in a list. 
#' 
#'
#' @param sample_list list, produced by function read_file.
#' @return Return a dataframe. 
.mymerge_file <- function(sample_list){
  #sample_file <- data.frame()
  for (i in 1:base::length(sample_list)) {
    if (is.data.frame(sample_list[[i]]) == 1) {
      sample_list[[i]]$uniq_id <- names(sample_list[i])
    }
  }
  sample_list <- purrr::map(sample_list, ~purrr::compact(.)) %>% purrr::keep(~base::length(.) != 0)  #Remove Null value in list
  sample_file <- Reduce("rbind.fill", sample_list)
  return(sample_file)
}

##' Classify the type of variants
#'
#' \code{merge_file} is a function to merge all dataframe in a list. 
#' 
#'
#' @param variant_df data.frame, produced by function sample_file.
#' @return Return a dataframe. 
.myclassify_variant <- function(variant_df){
  if (all(c("Mutation_type_(annovar)", "Mutation_type_(VEP)") %in% colnames(variant_df))) {
    variant_df <- within(variant_df, {
      class1 <- NA
      class2 <- NA
      class3 <- NA
      class4 <- NA
      class1[str_detect(`Mutation_type_(VEP)`, pattern = "intron_variant|splice_region_variant|start_retained_variant|stop_retained_variant|non_coding_transcript_exon_variant|3_prime_UTR_variant|5_prime_UTR_variant") | str_starts(`Mutation_type_(annovar)`, pattern = "intronic|ncRNA|UTR3|UTR5|upstream|downstream|intergenic")] <- 1
      class2[str_detect(`Mutation_type_(VEP)`, pattern = "synonymous_variant") | str_starts(`Mutation_type_(annovar)`, pattern = "synonymous")] <- 2
      class3[str_detect(`Mutation_type_(VEP)`, pattern = "missense_variant|inframe_deletion|inframe_insertion") | str_starts(`Mutation_type_(annovar)`, pattern = "nonsynonymous|nonframeshift")] <- 3
      class4[str_detect(`Mutation_type_(VEP)`, pattern = "splice_acceptor|splice_donor|stop_gained|frameshift_variant") | str_starts(`Mutation_type_(annovar)`, pattern = c("frameshift|splicing|stopgain"))] <- 4
    })
  }
  if ("Mutation_type_(VEP)" %in% colnames(variant_df) & !"Mutation_type_(annovar)" %in% colnames(variant_df)) {
    variant_df <- within(variant_df, {
      class1 <- NA
      class2 <- NA
      class3 <- NA
      class4 <- NA
      class1[str_detect(`Mutation_type_(VEP)`, pattern = "intron_variant|splice_region_variant|start_retained_variant|stop_retained_variant|non_coding_transcript_exon_variant|3_prime_UTR_variant|5_prime_UTR_variant")] <- 1
      class2[str_detect(`Mutation_type_(VEP)`, pattern = "synonymous_variant")] <- 2
      class3[str_detect(`Mutation_type_(VEP)`, pattern = "missense_variant|inframe_deletion|inframe_insertion")] <- 3
      class4[str_detect(`Mutation_type_(VEP)`, pattern = "splice_acceptor|splice_donor|stop_gained|frameshift_variant")] <- 4
    })
  }
  if(!"Mutation_type_(VEP)" %in% colnames(variant_df) & "Mutation_type_(annovar)" %in% colnames(variant_df)){
    variant_df <- within(variant_df,{
      class1 <- NA
      class2 <- NA
      class3 <- NA
      class4 <- NA
      class1[str_starts(`Mutation_type_(annovar)`, pattern = "intronic|ncRNA|UTR3|UTR5|upstream|downstream|intergenic")] <- 1
      class2[str_starts(`Mutation_type_(annovar)`, pattern = "synonymous")] <- 2
      class3[str_starts(`Mutation_type_(annovar)`, pattern = "nonsynonymous|nonframeshift")] <- 3
      class4[str_starts(`Mutation_type_(annovar)`, pattern = c("frameshift|splicing|stopgain"))] <- 4
    })
  }
  
  variant_df <- within(variant_df,{
    class <- NA
    class[pmax(class1,class2,class3,class4,na.rm = T)==4] <- "PTV" 
    class[pmax(class1,class2,class3,class4,na.rm = T)==3] <- "MIS"
    class[pmax(class1,class2,class3,class4,na.rm = T)==2] <- "SYN"
    class[pmax(class1,class2,class3,class4,na.rm = T)==1] <- "NON"
  })
  return(variant_df)
}

## Quanlity control
.myEssentialQc <- function(sample_file,test_type="panel"){
  ## only use panel gene
  if(test_type=='panel'){
    sample_file <- sample_file %>% dplyr::filter(Dis_gene %in% panel_gene$Gene)
  }
  ## remove low quanlity variant
  sample_file <- sample_file %>%
    dplyr::filter(FILTER == "PASS") %>%
    mutate(
      ref_reads = as.numeric(str_split(FORMAT, ",|:", simplify = T)[, 2]),
      allel_reads = as.numeric(str_split(FORMAT, ",|:", simplify = T)[, 3]),
      all_reads = as.numeric(str_split(FORMAT, ",|:", simplify = T)[, 4])
    ) %>%
    mutate(
      allel_freq = allel_reads / (ref_reads + allel_reads),
      true_reads_per = (ref_reads + allel_reads) / all_reads
    ) %>%
    dplyr::filter(allel_freq >= 0.3, true_reads_per > 0.5, Report_Variant != ".")
  
  return(sample_file)
}

get_burdenscore <- function(sample_file,gene){
  dt_score <- sample_file %>% dplyr::filter(Dis_gene %in% gene)
  dt_score <- within(dt_score,{
    class1 <- NA
    class2 <- NA
    class3 <- NA
    class4 <- NA
    class1[str_detect(`Mutation_type_(VEP)`,pattern="intron_variant|splice_region_variant|start_retained_variant|stop_retained_variant|non_coding_transcript_exon_variant|3_prime_UTR_variant|5_prime_UTR_variant") | str_starts(`Mutation_type_(annovar)`,pattern= "intronic|ncRNA|UTR3|UTR5|upstream|downstream|intergenic")] <- 1
    class2[str_detect(`Mutation_type_(VEP)`,pattern="synonymous_variant") | str_starts(`Mutation_type_(annovar)`,pattern= "synonymous")] <- 2
    class3[str_detect(`Mutation_type_(VEP)`,pattern="missense_variant|inframe_deletion|inframe_insertion") | str_starts(`Mutation_type_(annovar)`,pattern = "nonsynonymous|nonframeshift")] <- 3
    class4[str_detect(`Mutation_type_(VEP)`,pattern="splice_acceptor|splice_donor|stop_gained|frameshift_variant") | str_starts(`Mutation_type_(annovar)`,pattern = c("frameshift|splicing|stopgain"))] <- 4
  })
  dt_score <- within(dt_score,{
    class <- NA
    class[pmax(class1,class2,class3,class4,na.rm = T)==4] <- "PTV" 
    class[pmax(class1,class2,class3,class4,na.rm = T)==3] <- "MIS"
    class[pmax(class1,class2,class3,class4,na.rm = T)==2] <- "SYN"
    class[pmax(class1,class2,class3,class4,na.rm = T)==1] <- "NON"
  })
  dt_score <- dt_score[,c("Sample","Dis_gene","Report_Variant","class","revel")]
  dt_score <- dt_score[,c("Sample","Dis_gene","class")] %>% dplyr::filter(class %in% c('PTV','MIS')) %>%
    mutate(class_score=ifelse(class=='MIS',1,2)) %>% dplyr::group_by(Sample) %>% dplyr::summarize(score=sum(class_score))
  return(dt_score)
}


##' Collapsing Burden test.
#'
#' \code{collapsing_burden} is a function to make collapsing burden test.
#' 
#'
#' @param sample_file dataframe, a big dataframe contains all essential reports.
#' 
#' @return Return a dataframe. It's the final results of burden test.
collapsing_burden <- function(sample_file){
  case_file <- sample_file %>% dplyr::filter(diagnosis=="case")
  ctr_file <- sample_file %>% dplyr::filter(diagnosis=="control")
  case_size <- base::length(unique(case_file$uniq_id))
  ctr_size <- base::length(unique(ctr_file$uniq_id))
  variant_n_table <- rbind(case_file[,c("Dis_gene","Mutation_type_(VEP)","Mutation_type_(annovar)","diagnosis","uniq_id")],
                           ctr_file[,c("Dis_gene","Mutation_type_(VEP)","Mutation_type_(annovar)","diagnosis","uniq_id")]) %>% dplyr::rename(gene=Dis_gene) 
  variant_n_table <- within(variant_n_table,{
    class1 <- NA
    class2 <- NA
    class3 <- NA
    class4 <- NA
    class1[str_detect(`Mutation_type_(VEP)`,pattern="intron_variant|splice_region_variant|start_retained_variant|stop_retained_variant|non_coding_transcript_exon_variant|3_prime_UTR_variant|5_prime_UTR_variant") | str_starts(`Mutation_type_(annovar)`,pattern= "intronic|ncRNA|UTR3|UTR5|upstream|downstream|intergenic")] <- 1
    class2[str_detect(`Mutation_type_(VEP)`,pattern="synonymous_variant") | str_starts(`Mutation_type_(annovar)`,pattern= "synonymous")] <- 2
    class3[str_detect(`Mutation_type_(VEP)`,pattern="missense_variant|inframe_deletion|inframe_insertion") | str_starts(`Mutation_type_(annovar)`,pattern = "nonsynonymous|nonframeshift")] <- 3
    class4[str_detect(`Mutation_type_(VEP)`,pattern="splice_acceptor|splice_donor|stop_gained|frameshift_variant") | str_starts(`Mutation_type_(annovar)`,pattern = c("frameshift|splicing|stopgain"))] <- 4
  })
  variant_n_table <- within(variant_n_table,{
    class <- NA
    class[pmax(class1,class2,class3,class4,na.rm = T)==4] <- "PTV" 
    class[pmax(class1,class2,class3,class4,na.rm = T)==3] <- "MIS"
    class[pmax(class1,class2,class3,class4,na.rm = T)==2] <- "SYN"
    class[pmax(class1,class2,class3,class4,na.rm = T)==1] <- "NON"
  })
  
  # table(variant_n_table$class)
  variant_n_table <- variant_n_table[,c("gene","class","diagnosis","uniq_id")] %>% 
    dplyr::filter(!is.na(class)) %>% na.omit()
  
  variant_n_table <- variant_n_table %>% 
    dplyr::group_by(gene,diagnosis,class) %>% unique() %>% 
    dplyr::summarize(n=n()) %>% 
    mutate(variant_name = str_c(diagnosis,"_",class)) %>% 
    ungroup() %>% 
    dplyr::select(gene,variant_name,n) %>% 
    reshape2::dcast(.,gene~variant_name)
  variant_n_table[is.na(variant_n_table)] <- 0
  variant_n_table$case_n <- case_size
  variant_n_table$control_n <- ctr_size
  
  # one-sided Fisher's exact test
  res_list <- purrr::pmap(variant_n_table, function(gene,case_MIS,case_NON,case_PTV,case_SYN,control_MIS,
                                                    control_NON,control_PTV,control_SYN,case_n,control_n){
    tryCatch({
      PTV_Matrix <- matrix(c(case_PTV,control_PTV,case_n-case_PTV,control_n-control_PTV),ncol = 2)
      Fisher_PTVp <- fisher.test(PTV_Matrix,alternative = "greater")$p.value
      PTV_OR <- fisher.test(PTV_Matrix)[["estimate"]][["odds ratio"]]
      
      MIS_Matrix <- matrix(c(case_MIS,control_MIS,case_n-case_MIS,control_n-control_MIS),ncol = 2)    
      Fisher_MISp <- fisher.test(MIS_Matrix,alternative = "greater")$p.value
      MIS_OR <- fisher.test(MIS_Matrix)[["estimate"]][["odds ratio"]]
      
      SYN_Matrix <- matrix(c(case_SYN,control_SYN,case_n-case_SYN,control_n-control_SYN),ncol = 2)
      Fisher_SYNp <- fisher.test(SYN_Matrix,alternative = "greater")$p.value
      
      NON_Matrix <- matrix(c(case_NON,case_NON,case_n-case_NON,control_n-control_NON),ncol = 2)
      Fisher_NONp <- fisher.test(NON_Matrix,alternative = "greater")$p.value
    },error = function(e) {
      cat("Error")
    })
    res_df <- data.frame(Fisher_PTVp=Fisher_PTVp,Fisher_MISp=Fisher_MISp,
                         Fisher_SYNp=Fisher_SYNp,Fisher_NONp=Fisher_NONp,
                         PTV_OR=PTV_OR,MIS_OR=MIS_OR)
    return(res_df)
  })
  fisher_df <- do.call('rbind',res_list)
  variant_n_table <- cbind(variant_n_table,fisher_df)
  
  variant_n_table$Fisher_PTVq <- rep("NA",nrow(variant_n_table))
  variant_n_table$Fisher_MISq <- rep("NA",nrow(variant_n_table))
  variant_n_table$Fisher_SYNq <- rep("NA",nrow(variant_n_table))
  variant_n_table$Fisher_NONq <- rep("NA",nrow(variant_n_table))
  # multiple tests correction for the p-value by the number of genes in test.
  variant_n_table <- variant_n_table %>% mutate_at(c("Fisher_PTVp","Fisher_MISp","Fisher_SYNp","Fisher_NONp"),~as.numeric(.))
  
  variant_n_table$Fisher_PTVq <- p.adjust(variant_n_table$Fisher_PTVp, method = "BH", n = nrow(variant_n_table))
  variant_n_table$Fisher_MISq <- p.adjust(variant_n_table$Fisher_MISp, method = "BH", n = nrow(variant_n_table))
  variant_n_table$Fisher_SYNq <- p.adjust(variant_n_table$Fisher_SYNp, method = "BH", n = nrow(variant_n_table))
  variant_n_table$Fisher_NONq <- p.adjust(variant_n_table$Fisher_NONp, method = "BH", n = nrow(variant_n_table))
  # write output to file
  variant_n_table <- variant_n_table[,c("gene","case_PTV","case_MIS","case_SYN","case_NON",
                                        "control_PTV","control_MIS","control_SYN","control_NON",
                                        "case_n","control_n",
                                        "Fisher_PTVp","Fisher_MISp","Fisher_SYNp","Fisher_NONp",
                                        "Fisher_PTVq","Fisher_MISq","Fisher_SYNq","Fisher_NONq",
                                        "PTV_OR","MIS_OR")]
  return(variant_n_table)
}



