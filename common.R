
# Import data -------------------------------------------------------------

library(osfr)
osf_auth("pJq7p3qcpFA2noR8Sp7dDOUG3w1MmEyDpzsHGljrhBMoMuHASVLhwpLTWUDOJlw5zTeqk5")

# Create a temporary directory
tmp_dir <- tempdir()

osf_retrieve_file("https://osf.io/ybwrc/") %>%
  osf_download(
    path = tmp_dir)
file <- list.files(tmp_dir, pattern = "\\.rdata", full.names = TRUE)
# read in the file
load(file)
#WVS_data <- WV6_Data_R_v20201117


# Global options ----------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(reactable)

# Data - non transpose  ---------------------------------------------------
data_raw <- WV6_Data_R_v20201117

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
