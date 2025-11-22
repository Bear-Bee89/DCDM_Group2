library(dplyr)
library(readr)
library(stringr)

df_raw <- read_csv(
  "IMPC_parameter_description.txt",
  col_names = TRUE,
  trim_ws = FALSE
)

df_clean <- df_raw %>%
  # Trim whitespace for all columns
  mutate(across(everything(), ~str_trim(.x))) %>% 
  
  # Convert empty strings to NA **only for character columns**
  mutate(across(where(is.character), ~na_if(.x, ""))) %>% 
  
  # Remove duplicate rows
  distinct()

###
df_clean <- df_clean %>%
  mutate(across(
    where(is.character),
    ~ {
      x <- str_trim(.)                         # remove whitespace
      x <- ifelse(tolower(x) %in% c("na", "n/a", "null", ""), NA, x)
      x
    }
  ))

#####
df_clean <- df_clean %>%
mutate(across(where(is.character), ~ gsub("\\^", "", .)))
#####
write.csv(df_clean,
          file = "IMPC_parameter_description_cleaned.csv",
          row.names = FALSE)
######
any(duplicated(df_clean))
######