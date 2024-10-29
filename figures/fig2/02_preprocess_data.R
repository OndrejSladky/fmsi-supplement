#!/usr/bin/env Rscript

library(tidyverse)

df0 <- read_tsv("results.tsv", na = "na")

df1 <- df0 %>%
    mutate(genome = str_replace(genome, "sars-cov-2_pangenome_k32", "SC2-pang")) %>%
    mutate(genome = str_replace(genome, "spneumo_pangenome_k32", "Spneumo-pang")) %>%
    mutate(genome = str_replace(genome, "spneumoniae", "Spneumo-gen")) %>%
    mutate(genome = str_replace(genome, "escherichia_coli.k63", "Ecoli-pang")) %>%
    mutate(SIalg = paste(I_alg, "-", S_alg)) %>%
    mutate(bits_per_kmer=8000*Q_mem_kb/kmer_count) %>%
    mutate(bits_per_ref_kmer=bits_per_kmer*rate) %>%
    mutate(Q_time_per_query_ms = 1000*1000*Q_time_s / num_queries) %>%
    mutate(plot=ifelse(rate=="0.1", "subsampl", ifelse(qType=="Str", "stream", "isol")))

df2 <- df1 %>%
    mutate (alg = case_when(
        SIalg == "sbwt - none" ~ as.character(d+3), # (index with RCs)", # d \in [0,2] for SBWT
        SIalg == "cbl - none" ~ as.character(6),
        SIalg == "SSHash - prophasm" ~ as.character(2),
        SIalg == "bwa - prophasm" ~ as.character(1),
        SIalg == "prophex - prophasm" ~ as.character(0),
        SIalg == "fmsi - local" ~ as.character(d+7),
        SIalg == "fmsi - global" ~ as.character(7),
    )) %>%
    # mutate (algorithm = case_when(
    #     SIalg == "sbwt - none" ~ paste("SBWT ", ifelse(d==0, "(time efficient)", "(memory efficient)")), # (index with RCs)", the default varinat
    #     SIalg == "cbl - none" ~ "CBL",
    #     SIalg == "SSHash - prophasm" ~ "SSHash (on ProphAsm output)",
    #     SIalg == "bwa - prophasm" ~ "BWA (on ProphAsm output)",
    #     SIalg == "prophex - prophasm" ~ "ProPhex (on ProphAsm output)",
    #     SIalg == "fmsi - local" ~ paste("FMSI (on kmercamel's local with d=", d, ")", sep=""),
    #     SIalg == "fmsi - global" ~ "FMSI (on kmercamel's global)",
    # )) %>%
    mutate (algorithm = case_when(
        SIalg == "sbwt - none" ~ paste0("SBWT", ifelse(d==0, "-fast", "-small")), # (index with RCs)", the default varinat
        SIalg == "cbl - none" ~ "CBL",
        SIalg == "SSHash - prophasm" ~ "SSHash(spss)",
        SIalg == "bwa - prophasm" ~ "BWA(spss)",
        SIalg == "prophex - prophasm" ~ "ProPhex(spss)",
        SIalg == "fmsi - local" ~ "FMSI(spss)",
        SIalg == "fmsi - global" ~ "FMSI(superstr)",
    )) %>%
    mutate (qTypeLabel = case_when(
        qType == "Pos" ~ "+",
        qType == "Neg" ~ "-",
        qType == "Str" ~ "", # not used now
    ))

df2 %>%
    write_tsv("final_results.tsv")
