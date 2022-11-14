#library
##################################################
library(tidyverse)
library(ggpubr)
library(patchwork)
##################################################
.DFs_yshift  = function(x) {
    DF <- DFs[[x]]
    DF %>% 
        mutate(animal=x) %>% 
        dplyr::select(animal, cell_cell, member, yshift) -> return_object
    return(return_object)
}

.DFs_test  = function(x) {
    DF <- DFs[[x]]
    #### t-Test####
    T_Result <- t.test(yshift ~ member, data = DF)
    if(T_Result$p.value==0){
        T_pvalue <- 2.2e-16
    }else{
        T_pvalue <- T_Result$p.value
    }
    #### u-Test####
    U_Result <- wilcox.test(yshift ~ member, data = DF)
    if(U_Result$p.value==0){
        U_pvalue <- 2.2e-16
    }else{
        U_pvalue <- U_Result$p.value
    }
    #### F-Test####
    F_Result <- var.test(yshift ~ member, data = DF)
    if(F_Result$p.value==0){
        F_pvalue <- 2.2e-16
    }else{
        F_pvalue <- F_Result$p.value
    }
    p_value_res <- c(T_pvalue,U_pvalue,F_pvalue)
    q_value_res <- p.adjust(p_value_res,"BH")
    #### Ave.####
    ave_yshift <- mean(DF$yshift)
    #### 2SD####
    sd_2_yshift <- sd(DF$yshift)
    #### data.frame####
    data.frame(Animal = x,
               T_pvalue = p_value_res[1],
               T_qvalue = q_value_res[1],
               U_pvalue = p_value_res[2],
               U_qvalue = q_value_res[2],
               F_pvalue = p_value_res[3],
               F_qvalue = q_value_res[3],
               Ave_yshift = ave_yshift,
               SD_2_yshift = sd_2_yshift,
               stringsAsFactors = FALSE,
               row.names = NULL) -> return_object
    return(return_object)
}