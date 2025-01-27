#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 4
w <- 7 # 7.5 without the legend, 15 with legend
u <- "cm"


fn <- paste0("fig_minikraken.pdf")

# $ wc -c minikraken_20171*/database*
#  536870928 minikraken_20171019_8GB/database.idx
# 8053064728 minikraken_20171019_8GB/database.kdb
#  536870928 minikraken_20171101_4GB_dustmasked/database.idx
# 3758097436 minikraken_20171101_4GB_dustmasked/database.kdb
#

kraken4gb <- 536870928 + 3758097436
kraken8gb <- 536870928 + 8053064728

ps <- 1.8 # point size
fs <- 8 # font size

MAX_Y_VAL <- 52
MAX_X_VAL <- 45

df0 <- read_tsv("lookup_final_results.tsv", na = "na")

genomes <- c("minikr4gib", "minikr8gib")
algorithms <- c("SSHash", "FMSI")

df1 <- df0 %>%
    filter(genome %in% genomes)


df2 <- df1 %>%
    filter(qType == "Pos") %>%
    filter(plot == "isol") %>%
    filter(algorithm != "FMSI(spss)") %>%
    mutate(algorithm = str_replace(algorithm, "FMSI.*", "FMSI")) %>%
    mutate(algorithm = str_replace(algorithm, "^BWA.*", "BWA")) %>%
    mutate(algorithm = str_replace(algorithm, "SSHash.*", "SSHash")) %>%
    filter(algorithm %in% algorithms)

df3 <- df2 %>%
    select(genome, k, kmer_count, algorithm, index_bytes) %>%
    rename(core = index_bytes) %>%
    mutate(dict2taxid = kmer_count * 2)

df4 <- df3 %>%
    bind_rows(
        df3 %>% filter(genome == "minikr4gib", algorithm == "FMSI") %>%
            mutate(
                algorithm = "Kraken",
                core = kraken4gb,
                dict2taxid = 0
            )
    ) %>%
    bind_rows(
        df3 %>% filter(genome == "minikr8gib", algorithm == "FMSI") %>%
            mutate(
                algorithm = "Kraken",
                core = kraken8gb,
                dict2taxid = 0
            )
    )

df5 <- df4 %>%
    pivot_longer(
        cols = c("core", "dict2taxid"),
        names_to = "part",
        values_to = "bytes"
    ) %>%
    mutate(part = factor(part) %>% fct_rev()) %>%
    mutate(algorithm = factor(algorithm, levels = c("Kraken", "SSHash", "FMSI"))) %>%
    mutate(genome = ifelse(genome == "minikr8gib", "MiniKraken8GiB db", "MiniKraken4GiB db"))


df <- df5


theme_set(
    theme_bw(base_size = fs) +
        theme(
            #axis.text.x = element_text(angle = 45),
            panel.grid.major = element_line(color = "grey90", size = 0.25),
            #panel.grid.minor = element_line(color = "grey90", size = 0.25),
            panel.grid.minor = element_blank(),
            legend.key.size = unit(0.2, "cm"),
            legend.position = "top",
            legend.spacing = unit(0.0, "cm"),
            # Overall spacing between items
            legend.text = element_text(margin = margin(
                t = 0,
                r = 0,
                b = 0,
                l = 0
            )),
            legend.margin = margin(
                t = 0,
                r = 0,
                b = 0,
                l = 0
            ),
            legend.title = element_blank(),
            strip.text.x = element_text(hjust = 0),
            plot.margin = margin(
                t = 2,
                r = 2,
                b = -7,
                l = 1,
                unit = "pt"
            )
        )
)


df_labels <- df %>%
    group_by(genome, algorithm) %>%
    summarize(total_bytes = sum(bytes),
              bytes_per_kmer = sum(bytes / kmer_count)) %>%
    ungroup()

ggplot() +
    geom_bar(
        data = df,
        aes(
            x = algorithm,
            y = bytes / kmer_count,
            fill = algorithm,
            alpha = part
        ),
        stat = "identity"
    ) +
    geom_text(
        data = df_labels,
        aes(
            x = algorithm,
            y = bytes_per_kmer,
            label = paste(round(total_bytes / 1e9, 2), "GB")
        ),
        position = position_identity(),
        vjust = -0.5,
        size = 0.3 * fs,
        color = "black"
    ) +
    scale_fill_nejm() + facet_grid(. ~ genome) +
    guides(algorithm = guide_legend(nrow = 1), color = guide_legend(nrow = 1)) +
    scale_alpha_manual(values = c("dict2taxid" = 0.6, "core" = 1.0)) +
    scale_y_continuous(limits = c(0, 17.5), expand = c(0, 0)) +
    labs(y = "bytes per k-mer", x = "")



## == == == == == == == == == == == == == ==
## Step 5: Save plot --------------------------------------------------------------------
## == == == == == == == == == == == == == ==

ggsave(fn,
       height = h,
       width = w,
       unit = u)
