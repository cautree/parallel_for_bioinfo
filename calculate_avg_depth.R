library(dplyr)
library(purrr)
library(tidyr)


## read all the coverage output from samtools depth
file_path= paste("bam_coverage/", list.files("bam_coverage/"), sep="")
samples = stringr::str_replace_all (list.files("bam_coverage/"), ".coverage" ,""  )

read_function = function(x,y){
  
  df =readr::read_table(x, col_names = F)
  names(df) = c("amplicon","position", "depth")
  df =df %>% 
    mutate( sample = y) %>% 
    select( sample, everything())
  return(df)
}

df = purrr::map2_dfr( file_path, samples, read_function)

# all sequence depth
df_group_all = df %>% 
  group_by( sample, amplicon)%>% 
  summarise( avg_depth = mean(depth),
          len = n())

# first 25bp depth
df_first_25 = df %>% 
  tidyr::nest(data=-c(sample, amplicon) ) %>% 
  mutate( data2 = purrr::map( .$data, function(x){
    
    x = x %>% 
      arrange(position) %>% 
      mutate( order = 1:nrow(.))
    x = x %>% 
      filter( order <26)
    return(x)
    
  })) %>% 
  dplyr::select( sample, amplicon, data2) %>% 
  tidyr::unnest()


df_group_first25 = df_first_25 %>% 
  group_by( sample, amplicon)%>% 
  summarise( first25_avg_depth = mean(depth),
             first_len = n())

##last 25bp depth
df_last_25 = df %>% 
  tidyr::nest(data=-c(sample, amplicon) ) %>% 
  mutate( data2 = purrr::map( .$data, function(x){
    
    x = x %>% 
      arrange(position) %>% 
      mutate( order = 1:nrow(.))
    
    y = x[ (nrow(.)-24): nrow(.), ]

    return(y)
    
  })) %>% 
  dplyr::select( sample, amplicon, data2) %>% 
  tidyr::unnest()


df_group_last25 = df_last_25 %>% 
  group_by( sample, amplicon)%>% 
  summarise( last25_avg_depth = mean(depth),
            last_len = n())