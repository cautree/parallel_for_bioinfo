library(dplyr)
library(readxl)

path_list1 = list.files("data", pattern = "Set 4 i5 R")
fold_path = paste( "data",path_list1, sep="/" )

path_list2 = purrr::map( fold_path, list.files, pattern = "[a-z].xlsx" )
complete_xlsx_path = paste( fold_path, path_list2, sep="/")

sheet_list = purrr::map( complete_xlsx_path, readxl::excel_sheets )

## customized grep to be used in purrr::map
grep_function = function(x){
  res = grep( "_FASTQ",x, value =T)
  return(res)
}

fastq_sheet_list = unlist(purrr::map( sheet_list,grep_function  ))

info_df = data.frame( excel_path = rep(complete_xlsx_path,each=4),
                      sheet_name = fastq_sheet_list,
                      stringsAsFactors = F
                      )


read_data_function = function(x, y){
  
  df = readxl::read_excel( x, sheet = y)
  
  df_s = df %>% 
    select(  well, `PF Clusters` , `% of the plate`) %>% 
    rename( PF_Clusters = `PF Clusters`,
            pct_of_plate = `% of the plate` ) %>% 
    mutate( run_name = stringr::str_extract(x, "\\d{8}_[A-Za-z]{5,}-[A-Za-z]{4}")) %>%  ##regex to get special string
    mutate( i5 = stringr::str_extract(y, "[A-Z]\\d{2}")) %>% 
    select( run_name, i5, well, PF_Clusters, pct_of_plate )
  
  return(df_s)
  
}



df = purrr::map2_dfr( info_df$excel_path, info_df$sheet_name, read_data_function )

df_pp = df %>% 
  tidyr::nest( -c(run_name, i5 )  ) %>% 
  mutate( data2 = purrr::map(.$data, function(x){
    mean_plate = mean(x$pct_of_plate, na.rm =T)  
    median_plate = median(x$pct_of_plate, na.rm = T)  
    sd_plate = sd(x$pct_of_plate, na.rm =T) 
    x = x %>% 
      mutate( super_high_low = ifelse(   (pct_of_plate >  1.3333*mean_plate) | (pct_of_plate < 0.6667*mean_plate)  , 1, 0 ) ) 
    return(x)
    
  })) %>% 
  select( run_name, i5, data2) %>% 
  tidyr::unnest() 


df_table = df_pp %>% 
  tidyr::nest( -run_name) %>% 
  mutate( data2 = purrr::map(.$data, function(x){
    
    y = x %>% 
      filter( super_high_low ==1)
    ## table to data frame
    z = table(y$well, y$i5)
    z = as.data.frame(z)
    
    z = z %>% 
      tidyr::spread( Var2, Freq)
    names(z)[1] = "well"
    return(z)
    
  })) %>% 
  select(run_name, data2 ) %>% 
  tidyr::unnest()
  

df_table_s = df_table %>% 
  tidyr::gather( E08: F09, key="i5", value = "outlier_flag") %>% 
  tidyr::nest( -c(well, i5)) %>% 
  mutate( outlier_counts = purrr::map_dbl(.$data, function(x){
    
    sum = sum(x$outlier_flag, na.rm =T)
    x = x %>% 
      mutate( outlier_n = sum) %>% 
      select( -outlier_flag, -run_name) %>% 
      distinct() %>% 
      pull( outlier_n)
    return(x)
  })) %>% 
  dplyr::select( well, i5, outlier_counts) %>% 
  tidyr::unnest()

df_table_report = df_table_s %>% 
  tidyr::spread( i5, outlier_counts )

readr::write_csv(df_table_report, "results/report_of_outliers_in_each_position.csv" )
