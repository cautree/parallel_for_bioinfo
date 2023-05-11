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
