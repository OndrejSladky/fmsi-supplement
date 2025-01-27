#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

options(tibble.width = Inf, width = 300) # for printing

h <- 7
w <- 16 # 7.5 without the legend, 15 with legend
u <- "cm"

FMSIfasterMS <- "FMSI(faster MS)"

fn <- paste0("fig_constr_time.pdf")

ps <- 1.8 # point size
fs <- 8 # font size

MAX_Y_VAL <- 52
MAX_X_VAL <- 45

df0 <- read_tsv("final_results.tsv", na = "na")

#genomes <- c("hg-t2t", "ec-pg-hq")
ks <- c(31)
algorithms <- c("CBL", "SSHash", "SBWT-small", "SBWT-fast", "FMSI", FMSIfasterMS)

df1 <- df0 %>%
    filter((rate == 1.0 & genome == "hg-t2t") | (genome == "ec-pg-hq") | (genome == "mtg-ilm"))  %>%
    filter(qType == "Pos") %>%
    filter(k %in% ks) 
    #filter(genome %in% genomes)

df2 <- df1 %>%
    filter(algorithm != "FMSI(spss)") %>%
    mutate(
        algorithm = recode_factor(
            algorithm,
            "FMSI(superstr)"   = "FMSI",
            "SBWT-small"       = "SBWT-small",
            "SBWT-fast"        = "SBWT-fast",
            "SSHash(spss)"     = "SSHash",
            "CBL"              = "CBL"
        )
    ) %>%
    mutate(
        S_alg = recode_factor(
            S_alg,
            "global"    = "KmerCamel (global MS)",
            "eulertigs" = "GGCAT (eulertigs)",
        )
    ) %>%
    filter(algorithm %in% algorithms)

df3 <- df2 %>%
    bind_rows(
        df2 %>% filter(genome == "hg-t2t") %>% filter(algorithm == "FMSI") %>%
            mutate(
                algorithm = FMSIfasterMS,
                S_time_s = as.double(2800) ## GGCAT comp. of eulertigs (1706 s) + estimate of ~1000 s to comp. global MS from eulertigs (source: Karel :-)
            )
    ) %>%
    bind_rows(
        df2 %>% filter(genome == "ec-pg-hq") %>% filter(algorithm == "FMSI") %>%
            mutate(
                algorithm = FMSIfasterMS, # no constr. time update for subsampled
            )
    )

df4 <- df3 %>%
    #mutate(S_time_s = case_when(is.na(S_time_s) ~ 0, !is.na(S_time_s) ~ S_time_s)) %>%
    mutate(S_time_h = S_time_s / 3600) %>%
    mutate(I_time_h = I_time_s / 3600) %>%
    mutate(S_mem_GB = S_mem_kb * 1024 / 10^9) %>%
    mutate(I_mem_GB = I_mem_kb * 1024 / 10^9) %>%
    mutate(rate = ifelse(rate == 1.0, "100%", "10%")) %>%
    mutate_if(is.numeric, round, digits = 3)


df5 <- df4 %>%
    select(genome, rate, k, algorithm, S_alg, S_time_h, S_mem_GB, I_time_h, I_mem_GB)


df6 <- df5 %>%
    mutate(algorithm = fct_relevel(algorithm, algorithms))
#     pivot_longer(
#         cols = c("S_time_s", "I_time_s"),
#         names_to = "part",
#         values_to = "seconds"
#     ) %>%
#     mutate(part = factor(part) %>% fct_rev()) %>%

df <- df6

show(df)


# Install required packages if not already installed
if (!requireNamespace("kableExtra", quietly = TRUE)) install.packages("kableExtra")
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")

# Load necessary libraries
library(kableExtra)
library(knitr)


# Create a LaTeX table
create_latex_table <- function(df, output_file) {
  # Convert dataframe to a LaTeX table using knitr and kableExtra
  latex_table <- df %>%
    kable("latex", booktabs = TRUE, align = "c") %>% # Use booktabs for professional table style
    kable_styling(latex_options = c("hold_position")) # Keeps table at current position

  # Save the LaTeX table to a file
  cat(latex_table, file = output_file)

  message("LaTeX table saved to ", output_file)
}

# Call the function to create and save the LaTeX table
output_file <- "table_constr_memtime.tex"  # Specify the output file name
create_latex_table(df, output_file)
