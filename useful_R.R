## create many plots at the same time and save them in a folder
plot_df = df %>% 
  tidyr::pivot_longer( `TT105-2_A01`: `TT105-2_H12`, names_to = "well", values_to = "depth", values_drop_na = TRUE) %>% 
  tidyr::nest( data = -well) %>% 
  mutate(plot = purrr::map2(
    data, well, 
    ~ ggplot(data = .x, aes(x = position, y = depth)) +
      ggtitle(glue("Depth for {.y}")) +
      geom_line())) 
  
purrr::walk2(plot_df$plot, plot_df$well, ~{
     pdf(paste0("coverage_plots/TT105-2_md/",  .y, "_depth.pdf") )
     print(.x)
     dev.off()
})


## turn many columns of factor vars in to numeric at the same time
n_covars <- map_dfc(covars[,-c(1:2)],as.numeric)


## turn everything into numberic var first, then do the correlation and plot
corrplot(cor(n_covars[,-1],n_phenos),
         method = "circle",
         cl.cex=0.6,cl.ratio=0.5,cl.align.text = "l", tl.col = "black")
         
         
##  complex heatmap
n_exposome <- map_dfc(expsms[,-1],as.numeric)
cor_exposome <- cor(n_exposome)
htc <- ComplexHeatmap::Heatmap(cor_exposome,
                               name = "Cor",
                show_column_names = FALSE,
              show_row_names = FALSE)
ComplexHeatmap::draw(htc)

## prop.table
msitab <- with(coad$covars,table(MSI_status,sex))
prop.table(msitab,2)


## correlation plots
cormat <- cor(coad$covars[,1:5],use="pairwise.complete.obs")
corrplot::corrplot(cormat, method="square")


## is.element to check if the first arguemnt in the set (second argument)
is.element(coad$fdnam$Chr,c("X","Y"))

## get the rank of each element in a vector, then get the pisition that pass the threashhold
rfmad <- rank(-fmad)
fidx <- which(rfmad <= 500)


## get rid the first row
tail -n +2 nums.txt

## see both head and tail of the files
(head -n 2; tail -n 2) < Mus_musculus.GRCm38.75_chr1.bed


## read files, and as.data.frame at the same time
as.data.frame(fread(reportPath,stringsAsFactors=FALSE))


## use regex
lookaround <- "A penguin costs 2.99, a whale costs 5.99, I only have 3.50 left."
stringr::str_extract_all(lookaround, "\\d\\.\\d{2}")
#[1] "2.99" "5.99" "3.50"

## Look behind "(?<=if preceded by this)match_this"
prices <- str_extract_all(lookaround, "(?<=costs)\\s\\d\\.\\d{2}")
#[1] " 2.99" " 5.99"


## Look Aheads: "match this(?=if followed by this)"
animals <- str_extract_all(lookaround, "\\w+\\s(?=costs)")
#[1] "penguin " "whale "


##in R studio, if you highlight some text and do CMD-Option-Return on a Mac, the text gets sent to the RStudio Unix Terminal



# run model
race_model = function(df){
  aov(formula = meta_reading ~ race_category, data = df)
}

gender_model = function(df){
  aov(formula = meta_reading ~ gender, data = df)
}
df_mod_for_race = df %>% 
  tidyr::gather( `EPA_0`: `MCTR2_8.8`, key = "meta", value = "meta_reading")  %>% 
  tidyr::nest(-meta) %>% 
  mutate( mod = purrr::map(.$data, race_model)) %>% 
  mutate( mod_tidy= purrr::map(.$mod, broom::tidy)) %>% 
  select( meta, mod_tidy) %>% 
  tidyr::unnest() %>% 
  filter( term =="race_category") %>% 
  select(meta, p.value) %>% 
  mutate( response = "race") %>% 
  arrange(p.value)



  
df_mod_for_gender = df %>% 
  tidyr::gather( `EPA_0`: `MCTR2_8.8`, key = "meta", value = "meta_reading")  %>% 
  tidyr::nest(-meta) %>% 
  mutate( mod = purrr::map(.$data, gender_model)) %>% 
  mutate( mod_tidy= purrr::map(.$mod, broom::tidy)) %>% 
  select( meta, mod_tidy) %>% 
  tidyr::unnest() %>% 
  filter( term =="gender") %>% 
  select(meta, p.value) %>% 
  mutate( response = "gender") %>% 
  arrange(p.value)


df_gender_race = dplyr::bind_rows(df_mod_for_race, df_mod_for_gender)
df_gender_race = df_gender_race %>% 
  rename( term = meta) %>% 
  mutate(p.value = -1 * log10(p.value)) %>% 
  mutate( estimate = as.numeric(rescale(p.value))) %>% 
  dplyr::mutate(estimate = round(estimate,2))





## get p and cor from spearmn cor
correlation = psych::corr.test( df_non_case_s[-1], use = "complete.obs", method = "spearman")

cor_df = as.data.frame(correlation$r)
p_df = as.data.frame(correlation$p)

cor_df_long = cor_df %>% 
  tibble::rownames_to_column( var = "SPM1") %>% 
  gather(EPA_0:MCTR2_8.8, key="SPM2", value = "SPM_cor" ) %>% 
  filter( SPM1 != SPM2) %>% 
  filter( SPM1 %in% c("EPA_0","AA_0.1", "DHA_0", "LTB4_0", "TXB2_24.9"))

p_df_long = p_df %>% 
  tibble::rownames_to_column( var = "SPM1") %>% 
  gather(EPA_0:MCTR2_8.8, key="SPM2", value = "SPM_p" ) %>% 
  filter( SPM1 != SPM2) %>% 
  filter( SPM1 %in% c("EPA_0","AA_0.1", "DHA_0", "LTB4_0", "TXB2_24.9"))


plot_data = cor_df_long %>% 
  left_join( p_df_long, by = c("SPM1", "SPM2")) %>% 
  dplyr::mutate(SPM_cor = round(SPM_cor,2))
names(plot_data) = c("response", "term",     "estimate", "p.value" )
