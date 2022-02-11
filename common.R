# Import data -------------------------------------------------------------
#install.packages("arrow")
library(arrow)
# load("~/Desktop/bookdown-demo-main/data/WV6_Data_R_v20201117.rdata")
# write_parquet(WV6_Data_R_v20201117, here::here("data", "WVS6"))
data_parquet <- read_parquet(here::here("data", "WVS6"))

# Global options ----------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(reactable)

# Data - non transpose  ---------------------------------------------------
data_raw <- data_parquet

data = data_raw %>%
  filter(C_COW_ALPHA %in% c("USA","IND","NIG")) %>%
  rename("country" = C_COW_ALPHA) %>%
  dplyr::select(country, V4:V9, V12:V22, V45:V56, V70:V79, V95:V146, -V144, V153:V160J, V192:V216, V228A:V228K) %>%
  #remove cols which have all NA's
  janitor::remove_empty(which = "cols") %>%
  #remove cols which have 20% or more missing data (NA)
  purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=20) %>%
  #impute the missing values with
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>%
  remove_rownames()
