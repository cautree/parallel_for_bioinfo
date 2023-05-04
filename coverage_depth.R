library(dplyr)
library(purrr)
library(tidyr)

function_check = function(x){
  sum_thresh = sum(x >200, na.rm = T)
  return(sum_thresh)
}

path = "bam_coverage/low_conc_10ng/"

## read all the coverage output from samtools depth
file_path= paste(path, list.files(path), sep="")
samples = stringr::str_replace_all (list.files(path), ".coverage" ,""  )

read_function = function(x,y){
  
  df =readr::read_table(x, col_names = F)
  names(df) = c("amplicon","position", "depth")
  df_wide = df %>% 
    tidyr::spread( amplicon, depth)
  
  thresh = as.data.frame(sapply( df_wide[-1], function_check ))
  names(thresh) = "counts"
  thresh = thresh %>% 
    tibble::rownames_to_column( var = "amplicon") %>% 
    filter( counts >50)
  
  df_s = df_wide %>% 
    select( position, thresh$amplicon)
  
  names(df_s)[-1] = paste(y, names(df_s)[-1], sep="_")
  
  return(df_s)
  
}

