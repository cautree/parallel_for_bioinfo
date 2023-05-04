library(tidymodels)
library(tidyverse)
library(broom)

data(Orange)

Orange <- as_tibble(Orange)
Orange

table(Orange$Tree)


cor(Orange$age, Orange$circumference)
#> [1] 0.914


ggplot(Orange, aes(age, circumference, color = Tree)) +
  geom_line()


## get the correlation of age and circumference for differenc trees
cor_table = Orange %>% 
  group_by(Tree) %>%
  summarize(correlation = cor(age, circumference))

ct <- cor.test(Orange$age, Orange$circumference)
broom::tidy(ct)

## below two are the same
nested <- 
  Orange %>% 
  nest( data = -Tree)
nested <- 
  Orange %>% 
  nest(data = c(age, circumference))


## below two are the same
nested %>% 
  mutate(test = map(data, ~ cor.test(.x$age, .x$circumference)))

nested %>% 
  mutate(test = map(.$data, ~ cor.test(.x$age, .x$circumference)))


nested %>% 
  mutate(
    test = map(data, ~ cor.test(.x$age, .x$circumference)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  Orange %>% 
  nest(data = c(age, circumference)) %>% 
  mutate(
    test = map(data, ~ cor.test(.x$age, .x$circumference)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test)


######## regression model
lm_fit <- lm(age ~ circumference, data = Orange)
summary(lm_fit)

broom::tidy(lm_fit)



Orange %>%
  nest(data = c(-Tree)) %>% 
  mutate(
    fit = map(data, ~ lm(age ~ circumference, data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) %>% 
  select(-data, -fit)



data(mtcars)
mtcars <- as_tibble(mtcars)  # to play nicely with list-cols
mtcars

mtcars %>%
  nest(data = c(-am)) %>% 
  mutate(
    fit = map(data, ~ lm(wt ~ mpg + qsec + gear, data = .x)),  # S3 list-col
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) %>% 
  select(-data, -fit)



regressions <- 
  mtcars %>%
  nest(data = c(-am)) %>% 
  mutate(
    fit = map(data, ~ lm(wt ~ mpg + qsec + gear, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )


regressions %>% 
  select(tidied) %>% 
  unnest(tidied)


regressions %>% 
  select(glanced) %>% 
  unnest(glanced)
##output like this
# r.squared adj.r.squared sigma statistic  p.value    df    logLik   AIC   BIC deviance df.residual  nobs
# <dbl>         <dbl> <dbl>     <dbl>    <dbl> <dbl>     <dbl> <dbl> <dbl>    <dbl>       <int> <int>
# 1     0.833         0.778 0.291     15.0  0.000759     3  -0.00580  10.0  12.8    0.762           9    13
# 2     0.625         0.550 0.522      8.32 0.00170      3 -12.4      34.7  39.4    4.08           15    19

regressions %>% 
  select(augmented) %>% 
  unnest(augmented)