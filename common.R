# Import data -------------------------------------------------------------
#install.packages("arrow")
library(arrow)
# load("~/Desktop/bookdown-demo-main/data/WV6_Data_R_v20201117.rdata")
# write_parquet(WV6_Data_R_v20201117, here::here("data", "WVS6"))
data_raw <- read_parquet(here::here("data", "WVS6"))

# Global options ----------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(reactable)
library(psych)
library(hrbrthemes)

# Data - non transpose  ---------------------------------------------------
data_val <- data_raw %>%
  filter(C_COW_ALPHA %in% c("USA","IND","NIG")) %>%
  rename("country" = C_COW_ALPHA) %>%
  dplyr::select(country, V4:V9, V12:V22, V45:V56, V70:V79, V95:V146, -V144, V153:V160J, V192:V227, V228A:V228K) %>%
  #remove cols which have all NA's
  janitor::remove_empty(which = "cols") %>%
  #remove cols which have 20% or more missing data (NA)
  purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=20) %>%
  #impute the missing values with
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>%
  remove_rownames()

# Data Preparation
  # remove columns with near zero variance
#data_val <- data_val[,-caret::nearZeroVar(data_val)]


# Indices  ----------------------------------------------------------------

# 1. Emancipative values

     #>Choice
  # Toleration of Abortion: V204,
  # Toleration of Divorce: V205,
  # Toleration of Homosexuality: V203,
     #>Equality
  # Women’s Equality (Politics): -V51,
       #On the whole, men make better political leaders than women do.
  # Women’s Equality (Education): -V50,
       # A university education is more important for a boy than for a girl.
  # Women’s Equality (Jobs): -V45,
      # V45. When jobs are scarce, men should have more right to a job than women.
     #>Voice
  # Priority More Say (Local): - V226,
  # Priority More Say (National): - V227,
  # Freedom of Speech:
    #> Autonomy
  # Independence a Desired Quality: V12,
  # Obedience NOT a Desired Quality: -V21, [set 2 = 0]
  # Imagination a Desired Quality: V15


# emp_val = data %>%
#   mutate(V12 = replace(V21, V21 == 2, 0),
#          V51 = max(V51) - V51  +1,
#          V50 =  max(V50) - V50 +1,
#          V45 = max(V45) - V45 +1,
#          V226 = max(V226) - V226 +1,
#          V227 = max(V227) - V227 +1) %>%
#   rowwise %>%
#   mutate(choice = sum(V204, V205, V203),
#          equality= sum(V51, V50, V45),
#          voice = sum(V226, V227),
#          autonomy = sum(V12, V21, V15)) %>%
#   dplyr::select(choice, equality, voice, autonomy, country)
# #emp_val
#
# emp_val_noise = emp_val %>%
#   map_if(.,is.numeric,jitter) %>%
#   as.data.frame()
