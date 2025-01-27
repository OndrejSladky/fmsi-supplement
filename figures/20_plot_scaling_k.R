#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 7.5
#w <- 10 # 7.5 without the legend, 15 with legend
#w <- 30 # 7.5 without the legend, 15 with legend
w <- 14 # 7.5 without the legend, 15 with legend
u <- "cm"

#dataset <- "Ecoli-pang-HQ"
#dataset <- "Ecoli-pang-all"
#dataset <- "HG-T2T"
#dataset <- "Spneumo-pang"

fn <- paste0("fig_k_scaling.pdf")

ps <- 2.0 # point size
fs <- 11 # font size

MAX_Y_VAL <- 42
MAX_X_VAL <- 45

MARGIN <- 4


genomes <- c("ec-pg-hq", "mtg-ilm", "hg-t2t")

my_shapes <- c(15, 25, 2, 3, 10, 9, 19, 2)
my_colors <- c(1, 4, 4, 2, 3)
custom_palette <- pal_nejm("default")(max(my_colors))[my_colors]



df0 <- read_tsv("final_results.tsv", na = "na")

df1 <- df0 %>%
    mutate(bits_per_kmer_ind = 8 * index_bytes / kmer_count) %>%
    pivot_longer(
        cols = c("bits_per_kmer", "bits_per_kmer_ind"),
        names_to = "measurement",
        values_to = "bits_per_kmer"
    ) %>%
    mutate(measurement = ifelse(measurement == "bits_per_kmer_ind", "disk", "memory"))

df2 <- df1 %>%
    filter(genome %in% genomes) %>%
    mutate(genome = factor(genome, levels = genomes))


df3 <- df2 %>%
    filter(qType == "Pos") %>%
    filter(plot == "isol") %>%
    filter(algorithm != "FMSI(spss)") %>%
    filter(k>12) %>%
    mutate(algorithm = str_replace(algorithm, "FMSI.*", "FMSI")) %>%
    mutate(algorithm = str_replace(algorithm, "^BWA.*", "BWA")) %>%
    mutate(algorithm = str_replace(algorithm, "SSHash.*", "SSHash")) %>%
    filter(algorithm != "BWA") %>%
    mutate(
        algorithm = recode_factor(
            algorithm,
            "FMSI(superstr)"   = "FMSI",
            "SBWT-small"       = "SBWT-small",
            "SBWT-fast"        = "SBWT-fast",
            "SSHash(spss)"     = "SSHash",
            "CBL"              = "CBL"
        )
    )

df <- df3 %>%
    mutate(bits_per_kmer_adj = ifelse(bits_per_kmer <= MAX_Y_VAL, bits_per_kmer, MAX_Y_VAL +
                                          MARGIN / 2))


# Define interpolation sequence
k_min <- min(df$k, na.rm = TRUE)
k_max <- max(df$k, na.rm = TRUE)
k_interp <- seq(k_min, k_max, by = 0.01) # Adjust step size as needed


theme_set(
    theme_bw(base_size = fs) +
        theme(
            #axis.text.x = element_text(angle = 45),
            panel.grid.major = element_line(color = "grey90", size = 0.25),
            #panel.grid.minor = element_line(color = "grey90", size = 0.25),
            panel.grid.minor = element_blank(),
            legend.key.size = unit(0.4, "cm"),
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
            strip.text.x = element_text(hjust = 0)
        )
)

rect_df <- df %>%
    distinct(genome) %>%
    mutate(
        xmin = -Inf,
        xmax = Inf,
        ymin = MAX_Y_VAL,
        ymax = Inf # Set to the maximum y-axis limit
    )

ggplot() +
    geom_line(
        data = df,
        aes(
            x = k,
            y = bits_per_kmer,
            color = algorithm,
            alpha = measurement
        ),
        size = ps * 0.3
    ) +
    scale_x_continuous(
        "k",
        breaks = seq(-1, 400, 4),
        expand = c(0, 0)
    ) +
    geom_rect(
        data = rect_df,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax
        ),
        fill = "white",
        inherit.aes = FALSE
    ) +
    # Plot original data points
    geom_point(
        data = df,
        aes(
            x = k,
            y = bits_per_kmer_adj,
            color = algorithm,
            shape = algorithm,
            alpha = measurement
        ),
        size = ps
    ) +
    scale_y_continuous(
        "bits per k-mer",
        breaks = c(seq(0, MAX_Y_VAL, 8), MAX_Y_VAL + MARGIN / 2),
        labels = c(seq(0, MAX_Y_VAL, 8), paste0(">", MAX_Y_VAL)),
        limits = c(0, MAX_Y_VAL + MARGIN),
        expand = c(0, 0)
    ) +
    theme(legend.position = "top") +
    scale_shape_manual(values = my_shapes) +
    scale_color_manual(values = custom_palette) +
    facet_grid(. ~ genome) +
    geom_hline(yintercept = MAX_Y_VAL, linewidth = 0.2 * ps, ) +
    geom_hline(
        yintercept = 2,
        linewidth = 0.3 * ps,
        color = "gray",
        linetype = "dashed"
    ) +
    guides(
        algorithm = guide_legend(nrow = 1),
        color = guide_legend(nrow = 1),
        alpha = guide_legend(nrow = 1)
    )



## == == == == == == == == == == == == == ==
## Step 5: Save plot --------------------------------------------------------------------
## == == == == == == == == == == == == == ==

ggsave(fn,
       height = h,
       width = w,
       unit = u)
