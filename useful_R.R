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