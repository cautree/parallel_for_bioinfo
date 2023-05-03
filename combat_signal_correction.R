library(dplyr)
library(sva)

df_orig = data.table::fread("/an/vital/metabolomics/QA/yliu/jupiter_meta/check_data/jupyter_experimental_data_20200226.csv")
df_orig = as.data.frame(df_orig)

#get baseline data for different batch
df_batch1 = df_orig %>% 
  dplyr::filter( batch==1) %>% 
  dplyr::filter( year==0)

df_batch2 = df_orig %>% 
  dplyr::filter( batch==2) %>% 
  dplyr::filter( year==0)

df = dplyr::bind_rows(df_batch1, df_batch2 )

get_imputated_sample= function(df_meta){
  
  ##get meta data, which has number in the name
  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  ##get other info data, which do not have number in the name
  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]
  
  #remove more than 20% missing
  df_data <- df_data[,colSums(is.na(df_data))< nrow(df_data)*.2]
  df_data[] <- lapply(df_data, NA2mean)
  
  res = dplyr::bind_cols( df_info, df_data)
  return(res)
  
}


#log transform, centered to median, sd =1
get_normalized_data = function(df_meta){
  
  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]
  
  #remove more than 20% missing
  df_data <- df_data[,colSums(is.na(df_data))< nrow(df_data)*.2]
  
  df_data = as.data.frame(apply( df_data,
                                 2,
                                 function(y) ( log(y) - median(log(y), na.rm = T))/ sd( log(y), na.rm = T)))
  
  res = dplyr::bind_cols( df_info, df_data)
  return(res)
  
}

NA2mean <- function(x) replace(x, is.na(x), min(x, na.rm = TRUE)*0.25)

get_combat_corrected_data= function(df_meta){
  
  df_meta = get_imputated_sample(df_meta)
  df_meta = get_normalized_data(df_meta)
  
  df_info = df_meta[ names(df_meta)[!grepl("[0-9]", names(df_meta))] ]
  df_data = df_meta[ names(df_meta)[grepl("[0-9]", names(df_meta))] ]
  
  df_info = df_info %>%
    tidyr::separate( plate_well, c("plate", "well"), sep="_", remove = FALSE) %>%
    dplyr::mutate(plate = as.numeric(plate))
  
  batch = df_info$plate
  
  edata = df_meta[names(df_meta) %in% c('subjectId', names(df_data))]
  
  rownames(edata) = edata$subjectId
  edata$subjectId = NULL
  
  edata = t(edata)
  
  combat_edata = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  combat_edata = as.data.frame(t(combat_edata))
  
  res = dplyr::bind_cols(df_info, combat_edata)
  
  res =res %>%
    dplyr::select(-plate, -well) %>%
    dplyr::select(subjectId, plate_well, everything())
  
  return(res)
  
}



df_corrected = get_combat_corrected_data(df)
