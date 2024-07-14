library(vegan)
library(tidyverse)

######## Cellular components ##########

#load data
data <- read.csv('../data/cellular_component_investment_data.csv')
colnames(data)[1] <- "sample.id"

metadata <- read.csv('../data/mtstp_analysis_metadata.tsv', sep='\t')
metadata <- subset(metadata, infection.status == "infected")

#find failed sample and remove from metadata
metadata <- metadata[metadata$sample.id != "mtstpEiu113", ]

test <- metadata %>%
  left_join(data, by = c("sample.id"))

trt <- test %>%
  select(c("sample.id","developmental.stage", "plant"))

data2 <- test %>%
  select(-c("sample.id", "plant", "infection.status", "developmental.stage", "lineage"))

euclidian.matrix <- dist(data2)

ad.test <- adonis2(euclidian.matrix ~ developmental.stage, data = trt)
ad.test

######## Biological processes ##########

#load data
data <- read.csv('../data/biological_processes_investment_data.csv')
colnames(data)[1] <- "sample.id"

metadata <- read.csv('../data/mtstp_analysis_metadata.tsv', sep='\t')
metadata <- subset(metadata, infection.status == "infected")

#find failed sample and remove from metadata
metadata <- metadata[metadata$sample.id != "mtstpEiu113", ]

test <- metadata %>%
  left_join(data, by = c("sample.id"))

trt <- test %>%
  select(c("sample.id","developmental.stage", "plant"))

data2 <- test %>%
  select(-c("sample.id", "plant", "infection.status", "developmental.stage", "lineage"))

euclidian.matrix <- dist(data2)

ad.test <- adonis2(euclidian.matrix ~ developmental.stage, data = trt)
ad.test
