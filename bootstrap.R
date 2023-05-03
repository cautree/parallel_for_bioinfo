

#delete observations from the data with HII_bl>1.0.

# use thresholds between 0 and 1.0. Use lets say 0, 0.05, 0.1, 0.15, 0.2, 0.25 etc.
# 
# bootstrap the data.
# 
# use one of the thresholds to dichotomize HII above(=1) or below it(=0).
# 
# predict the outcome~HII_2cat+other_predictors. You do not use clogit but logistic
# regression. Then you run it on bootstrapped data, save the coefficients
# 
# calculate predicted probability of the outcome for the out of bag bootstrap data [OOB data] (subjects that were not selected by bootstrap).
# 
# use that predicted probabilities to calculate AUC for the OOB data.
# 
# Repeat 4.-7. for the next threshold, etc until you run out of thresholds.
# 
# repeat the procedure 3.-8. for all thresholds. You do this for B=500 times and will end up with 500 AUCs for each threshold.
# 
# calc the median AUC for each threshold and select the threshold that produces the highest median AUC.


library(dplyr)
library(rsample)
library(pROC)
library(AUC)
library(broom)
library(purrr)


##====================setting up
keep_var =c("hdl_caco","HIIb","age","drug","race3cat","sbpb","smokecig","bmib","glucb",
            "ldlcb","fxchd","lntrigb")

df = readr::read_csv("june5.csv")

df_small = df %>% 
  dplyr::filter(HIIb <=1) %>%  ## select only hiib less or each to 1
  dplyr::select(keep_var)

df_small = df_small %>% 
  dplyr::mutate(hiib_cut_a = cut(HIIb,breaks = seq(0.2, 1, by = 0.05)))

table(df_small$hiib_cut_a)
levels(df_small$hiib_cut_a)

df_small$hiib_cut =df_small$hiib_cut_a
df_small$hiib_cut
levels(df_small$hiib_cut) = 1:16  # change the cut to numeric, for easy handling later
df_small$hiib_cut = as.numeric(df_small$hiib_cut)

max(df_small$hiib_cut)
min(df_small$hiib_cut)

## for later merge purpose
thresh_infor = df_small %>% 
  dplyr::select(hiib_cut, hiib_cut_a) %>% 
  dplyr::distinct(hiib_cut, hiib_cut_a) %>% 
  dplyr::arrange(hiib_cut) %>% 
  dplyr::rename(threshold = hiib_cut)

thresh_infor 

##==================test using 12:(0.75,0.8] as cut
# thresh is numeric 1 to 16, 1 is for (0.75,0.8], etc
boots = rsample::bootstraps(df_small, times =5)
names(boots)

fit_glm_on_bootstrap <- function(split, thresh) {
  #get the bootstrapping sample
  df = analysis(split)
  
  #use one of the thresholds to dichotomize HII above(=1) or below it(=0).
  df = df %>%
    dplyr::mutate( HIIb_cat = ifelse(hiib_cut<=thresh, 0,1)) %>% 
    dplyr::mutate(HIIb_cat = as.factor(HIIb_cat))
  
  ## the glm model for classification  
  glm(hdl_caco~ HIIb_cat +age+ drug+ as.factor(race3cat)+ sbpb+  smokecig+ bmib + glucb +ldlcb+ fxchd + lntrigb, df, family = binomial)
}

## get the model from the bootstrapping samples
boot_models <- boots %>% 
  dplyr::mutate(model = purrr::map(splits, fit_glm_on_bootstrap,12), #12 means cut of (0.75,0.8]
                coef_info = purrr::map(model, broom::tidy))

## get the model coefs(not used in this code, but can be useful to check)
boot_coefs <- boot_models %>% 
  tidyr::unnest(coef_info)
#### get the out of bag sample

get_assessment <- function(split, thresh) {
  
  
  df = assessment(split)
  
  #use one of the thresholds to dichotomize HII above(=1) or below it(=0).
  df= df%>%
    dplyr::mutate( HIIb_cat = ifelse(hiib_cut<=thresh, 0,1)) %>% 
    dplyr::mutate(HIIb_cat = as.factor(HIIb_cat))
}

## get the oob for each bootstrap
assessment_df = boots %>%
  dplyr::mutate(oob = purrr::map(splits, get_assessment,12)) %>% #12 means cut of (0.75,0.8]
  dplyr::select(oob, id) %>% 
  dplyr::rename(oob_id = id)


## merge model and oob
auc_df =dplyr::bind_cols(boot_models, assessment_df)


auc_df_report = auc_df %>%
  dplyr::mutate( fit.predict = purrr::map2(.$model, .$oob, function(x,y) predict(x, newdata=y, type="response") ))%>%
  dplyr::mutate(auc = purrr::map2(.$fit.predict, .$oob, function(x,y) broom::tidy(AUC::roc(x, as.factor(y$hdl_caco) )) )) %>% 
  dplyr::mutate(auc_value_pre = purrr::map2(.$fit.predict, .$oob, function(x,y) pROC::roc(y$hdl_caco, x) )) %>% 
  dplyr::mutate(auc_value = purrr::map_dbl(.$auc_value_pre, function(x) pROC::auc(x))) %>% 
  dplyr::select(id, auc_value)
###=================end the test=========


##function for all threshold=========================
# thresh is numeric 1 to 16, 1 is for (0.2,0.25], etc
get_auc = function( thresh, n_boot=500) {
  
  ## bootstrapping n_boot times
  set.seed(43)
  boots = rsample::bootstraps(df_small, times =n_boot)
  
  ## get the model in the bootstrapping samples
  boot_models <- boots %>% 
    dplyr::mutate(model = purrr::map(splits, fit_glm_on_bootstrap,thresh),
                  coef_info = purrr::map(model, broom::tidy))
  
  ## get the oob for each bootstrap
  assessment_df = boots %>%
    dplyr::mutate(oob = purrr::map(splits, get_assessment,thresh)) %>% 
    dplyr::select(oob)
  
  ## merge model and oob
  auc_df =dplyr::bind_cols(boot_models, assessment_df)
  
  ## using the model to get the prediction for oob, and calculate auc
  auc_df_report = auc_df %>%
    dplyr::mutate( fit.predict = purrr::map2(.$model, .$oob, function(x,y) predict(x, newdata=y, type="response") ))%>%
    dplyr::mutate(auc = purrr::map2(.$fit.predict, .$oob, function(x,y) broom::tidy(AUC::roc(x, as.factor(y$hdl_caco) )) )) %>% 
    dplyr::mutate(auc_value_pre = purrr::map2(.$fit.predict, .$oob, function(x,y) pROC::roc(y$hdl_caco, x) )) %>% 
    dplyr::mutate(auc_value = purrr::map_dbl(.$auc_value_pre, function(x) pROC::auc(x))) %>% 
    dplyr::select(id, auc_value) %>% 
    dplyr::mutate(threshold = thresh)
  
}

# thresh from 3((0.3,0.35]) to 14((0.85,0.9])
thresh_vec = 3:14

res =thresh_vec %>%
  purrr::map_dfr(., get_auc) %>% 
  dplyr::group_by (threshold) %>% 
  dplyr::summarise(median_auc = median(auc_value, na.rm = T)) %>%   # get the auc median
  dplyr::left_join(thresh_infor, by = "threshold")

