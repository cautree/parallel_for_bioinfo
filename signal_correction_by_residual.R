library(dplyr)
library(broom)
library(purrr)
library(tidyr)



lm_model = function(x){
  
  lm(meta_reading ~ plate, data =x)
}

df = df_meta %>%
  tidyr::gather( key="mzid", value= "meta_reading", -plate_well) %>%
  tidyr::separate(plate_well, c( "plate", "well"), sep="_", remove = FALSE) %>% 
  dplyr::select(-well) %>%
  dplyr::select(plate_well, plate, everything()) %>%
  tidyr::nest(-mzid) %>%
  dplyr::mutate( mod1 = purrr::map(.$data, ~lm_model(.x))) %>%
  dplyr::mutate(tidy_mod1 = purrr::map(.$mod1, ~broom::augment(.x))) %>%
  dplyr::mutate(res = purrr::map2(.$data, .$tidy_mod1, function(x,y){
    ## original data for late compare
    x_s = x %>%
      dplyr::rename(orig_plate = plate,
                    orig_meta_reading = meta_reading)
    ## residualized data, and original data
    y_s = y %>%
      dplyr::select(meta_reading, plate, .resid)
    
    z = dplyr::bind_cols(x_s, y_s)
    return(z)
    
  }) )%>% 
  dplyr::select(mzid, res) %>%
  tidyr::unnest() # can check for the original data and corrected ones matched well or not


df_corrected = df %>%
  dplyr::select(mzid, plate_well, .resid) %>%
  tidyr::spread( key= mzid, value=.resid)

