#!/usr/bin/env Rscript

library(tidyverse)

#df0 <- read_tsv("results.tsv", na = "na")
# Load necessary library
library(dplyr)

# Get a list of files with the pattern "exp_01_build_index_results.*.tsv"
file_list <- list.files(pattern = "^exp_03_lookup_results..*\\.tsv$")

# Initialize an empty list to store dataframes
df_list <- list()

# Iterate over the files and load them into dataframes
for (file in file_list) {
  # Read the TSV file into a dataframe
  df <- read.delim(file, header = TRUE, sep = "\t")

  # Append the dataframe to the list
  df_list[[file]] <- df
}

df0 <- bind_rows(df_list)

df1 <- df0 %>%
    mutate(genome = str_replace(genome, "minikraken4GB_k31", "minikr4gib")) %>%
    mutate(genome = str_replace(genome, "minikraken8GB_k31", "minikr8gib")) %>%
    mutate(SIalg = paste(I_alg, "-", S_alg)) %>%
    mutate(bits_per_kmer=8000*Q_mem_kb/kmer_count) %>%
    mutate(bits_per_ref_kmer=bits_per_kmer*rate) %>%
    mutate(Q_time_per_query_ms = 1000*1000*Q_time_s / num_queries) %>%
    mutate(plot=ifelse(rate=="0.1", ifelse(qType=="Str", "subsampl-stream", "subsampl-isol"), ifelse(qType=="Str", "stream", "isol")))

df2 <- df1 %>%
    mutate (alg = case_when(
        SIalg == "fmsi - global" ~ as.character(2),
        SIalg == "SSHash - eulertigs" ~ as.character(1),
        SIalg == "SSHash - prophasm" ~ as.character(1),
    )) %>%
    mutate (algorithm = case_when(
        SIalg == "fmsi - global" ~ "FMSI(superstr)",
        SIalg == "SSHash - eulertigs" ~ "SSHash(spss)",
        SIalg == "SSHash - prophasm" ~ "SSHash(spss)", # in a few cases, eulertigs were infeasible to compute using ggcat (HG T2T for k=17 and minikraken8GB)
    )) %>%
    mutate (qTypeLabel = case_when(
        qType == "Pos" ~ "+",
        qType == "Neg" ~ "-",
        qType == "Str" ~ "", # not used now
    ))

df2 %>%
    write_tsv("lookup_final_results.tsv")
