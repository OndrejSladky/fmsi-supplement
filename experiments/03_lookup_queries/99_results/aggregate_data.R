#!/usr/bin/env Rscript

library(tidyverse)

options(tibble.width = Inf, width = 300) # for printing

# merge size and memtime stats ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (dataset name).", call.=FALSE)
}
suffix=""
if (length(args)==2) {
    suffix=args[2]
}
dataset=args[1]

get_hostname <- function(){ 
    return(as.character(Sys.info()["nodename"])) 
} 
filesuffix=paste(".", dataset, ".", get_hostname(), suffix, ".tsv", sep="")


df.size_stats <- read_tsv(paste("size_stats", filesuffix, sep="")) %>%
    mutate(index_bytes_nostreaming = index_bytes)
    
# first for FMSI -----------------------------------------------------------

write("-- FMSI + CAMEL --", stdout())

df.camel_memtime <- read_tsv(paste("camel_memtime", filesuffix, sep="")) %>%
    mutate(S_time_s = `user(s)`+`sys(s)`) %>%
    mutate(S_mem_kb = `max_RAM(kb)`) %>%
    select(genome, rate, S_alg, k, d, S_time_s, S_mem_kb)
#print(df.camel_memtime)
df.fmsi_memtime <- read_tsv(paste("fmsi_memtime", filesuffix, sep="")) %>%
    mutate(I_time_s = `user(s)`+`sys(s)`) %>%
    mutate(I_mem_kb = `max_RAM(kb)`) %>%
    select(genome, rate, S_alg, k, d, I_time_s, I_mem_kb)
df.fmsi_query_memtime <- read_tsv(paste("fmsi_query_memtime", filesuffix, sep="")) %>%
    mutate(Q_time_s = `user(s)`+`sys(s)`) %>%
    mutate(Q_mem_kb = `max_RAM(kb)`) %>%
    select(genome, rate, S_alg, k, d, Q_time_s, Q_mem_kb, qType, num_queries)

df.fmsi_stats0 <- df.size_stats %>%
    filter(I_alg == "fmsi") %>%
    full_join(df.camel_memtime) %>%
    full_join(df.fmsi_memtime) %>%
    full_join(df.fmsi_query_memtime)
    
df.fmsi_stats <- df.fmsi_stats0 %>%
    mutate(SI_time_s = S_time_s + I_time_s) %>%
    mutate(SI_mem_kb = apply( df.fmsi_stats0[c('S_mem_kb', 'I_mem_kb')], 1, max )) %>%
    arrange(genome, rate, k, S_alg, d, qType) 
#show(df.fmsi_stats)
    

# rSPSS -----------------------------------------------------------

write("-- rSPSS --", stdout())

df.rspss_comp_memtime <- read_tsv(paste("rspss_comp_memtime", filesuffix, sep="")) %>%
    mutate(S_time_s = `user(s)`+`sys(s)`) %>%
    mutate(S_mem_kb = `max_RAM(kb)`) %>%
    select(genome, rate, S_alg, k, d, S_time_s, S_mem_kb)




# SSHash -----------------------------------------------------------

write("-- SSHash --", stdout())

df.sshash_memtime <- read_tsv(paste("sshash_memtime", filesuffix, sep="")) %>%
    mutate(I_time_s = `user(s)`+`sys(s)`) %>%
    mutate(I_mem_kb = `max_RAM(kb)`) %>%
    select(genome, rate, S_alg, k, d, I_time_s, I_mem_kb)
df.sshash_query_memtime <- read_tsv(paste("sshash_query_memtime", filesuffix, sep="")) %>%
    mutate(Q_time_s = `user(s)`+`sys(s)`) %>%
    mutate(Q_mem_kb = `max_RAM(kb)`) %>%
    select(genome, rate, S_alg, k, d, Q_time_s, Q_mem_kb, qType, num_queries)

df.sshash_stats0 <- df.size_stats %>%
    filter(I_alg == "SSHash") %>%
    full_join(df.rspss_comp_memtime) %>%
    full_join(df.sshash_memtime) %>%
    full_join(df.sshash_query_memtime)
    
df.sshash_stats <- df.sshash_stats0 %>%
    mutate(SI_time_s = S_time_s + I_time_s) %>%
    mutate(SI_mem_kb = apply( df.sshash_stats0[c('S_mem_kb', 'I_mem_kb')], 1, max )) %>%
    arrange(genome, rate, k, S_alg, d, qType) 

#show(df.sshash_stats)


write("-- final binding --", stdout())

df.stats <- df.fmsi_stats %>%
    bind_rows(df.sshash_stats)%>%
    mutate(index_bytes = case_when(qType == "Str" ~ index_bytes_streaming, qType != "Str" ~ index_bytes_nostreaming)) %>%
    select(!c(index_bytes_streaming, index_bytes_nostreaming)) %>%
    mutate(Q_time_mus_per_query = Q_time_s / num_queries * 10^6) %>%
    mutate(Q_mem_bits_per_kmer = Q_mem_kb * 1000 * 8 / kmer_count) %>%
    arrange(genome, rate, k, I_alg, S_alg, d, qType)
    
df.stats %>% 
    write_tsv(paste("exp_03_lookup_results", filesuffix, sep=""),  na = "na")
