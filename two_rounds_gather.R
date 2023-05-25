## calculate the other SPM correlation with the SPM in the 5 list
plot_data = df_non_case_s %>% 
  select( subjectId, EPA_0, AA_0.1, DHA_0, LTB4_0, TXB2_24.9, everything()) %>% 
  tidyr::gather( EPA_0: TXB2_24.9, key="SPM_5", value = "SPM_5_value" ) %>% 
  tidyr::gather( `15-HEPE_0`:`MCTR2_8.8`, key = "other_SPM", value = "other_SPM_value") %>% 
  tidyr::nest(-c(SPM_5, other_SPM)) %>% 
  mutate( data2 = purrr::map(.$data, function(x){
     
    res = broom::tidy(cor.test(x$SPM_5_value, x$other_SPM_value, method=c( "spearman")))
    res = res %>% 
      select( estimate, p.value)
    return(res)
  })) %>% 
  select(SPM_5, other_SPM, data2 ) %>% 
  tidyr::unnest()
