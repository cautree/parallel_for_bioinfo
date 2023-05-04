library(boot)
library(table1)
library(flextable)

melanoma2 <- melanoma

# Factor the basic variables that
# we're interested in
melanoma2$status <- 
  factor(melanoma2$status, 
         levels=c(2,1,3),
         labels=c("Alive", # Reference
                  "Melanoma death", 
                  "Non-melanoma death"))

res = table1(~ factor(sex) + age + factor(ulcer) + thickness | status, data=melanoma2)
write.table (res , "my_table_1_file.csv", col.names = T, row.names=F, append= T, sep=',')
t1flex(res) %>% 
  save_as_docx(path="table1.docx")