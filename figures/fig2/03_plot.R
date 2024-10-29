#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 6
w <- 15 # 7.5 without the legend, 15 with legend
u <- "cm"

ps <- 2.0 # point size
fs <- 10 # font size


# == == == == == == == == == == == == == == == == ==
# Options (for copy-paste below) ----
# == == == == == == == == == == == == == == == == ==

variables <- "ktext ~ plot"
variables <- "genome ~ plot"

datasets <- c("Spneumo-pang")
datasets <- c("SC2-pang")
datasets <- c("Ecoli-pang", "Spneumo-pang", "SC2-pang")
datasets <- c("Ecoli-pang")

ks <- c(15)
ks <- c(23)
ks <- c(31)
ks <- c(15, 23, 31)


# == == == == == == == == == == == == == == == == ==
# Current plot combo configuration ----
#    ..adjust for current plot here...
# == == == == == == == == == == == == == == == == ==

# Small panel (not used)
#variables <- "genome ~ plot"; datasets <- c("Ecoli-pang"); ks <- c(23); h <- 6; w <- 15

# Fig 2, species x modes for k=23
variables <- "genome ~ plot"; datasets <- c("Ecoli-pang", "Spneumo-pang", "SC2-pang"); ks <- c(23); h <- 12; w <- 12

# Sup. Fig 1; Spneumo: modes x k
variables <- "ktext ~ plot"; datasets <- c("Spneumo-pang"); ks <- c(15, 23, 31); h <- 12; w <- 12
variables <- "ktext ~ plot"; datasets <- c("Ecoli-pang"); ks <- c(15, 23, 31); h <- 12; w <- 12
variables <- "ktext ~ plot"; datasets <- c("SC2-pang"); ks <- c(15, 23, 31); h <- 12; w <- 12










# ==
# Plotting -----
# ==

## == == == == == == == == == == == == == ==
## Step 1: Prepare variables, functions, etc. --------------------------------------------------------------------
## == == == == == == == == == == == == == ==

comboname <- gsub(" ", "", variables)  # Replace spaces with underscores
comboname <- gsub("~", "-vs-", comboname)

basename <- paste("plot",
                  comboname,
                  paste(datasets, collapse = "-"),
                  paste(ks, collapse = "-"),
                  sep = "__")

fn <- paste0(basename, ".pdf")

number_formatter <- function(y) {
    ifelse(y %% 1 == 0, as.character(as.integer(y)), as.character(y))
}

my_shapes <- c(9, 12, 19, 15, 2, 25, 10)



## == == == == == == == == == == == == == ==
## Step 2: Load data ----
## == == == == == == == == == == == == == ==

df0 <- read_tsv("final_results.tsv", na = "na")

df <- df0 %>%
    mutate(genome = fct_rev(genome)) %>%
    #filter(k == 23) %>%
    filter(genome %in% datasets) %>%
    filter(k %in% ks) %>%
    #filter(bits_per_kmer < 512) %>%
    filter(!(plot == "subsampl" & I_alg == "bwa")) %>%
    #mutate(plot2 = factor(plot, levels = c("stream", "isol", "subsampl"))) %>%
    mutate (
        plot = case_when(
            plot == "stream" ~ "a) streamed k-mer qs",
            plot == "isol" ~ "b) isolated k-mer qs",
            plot == "subsampl" ~ "c) subsampl. ref (10%)"
        )
    ) %>%
    mutate(ktext=paste0("k=",k))



## == == == == == == == == == == == == == ==
## Step 3: Set up theme --------------------------------------------------------------------
## == == == == == == == == == == == == == ==

theme_set(
    theme_bw(base_size=fs) +
        theme(
            axis.text.x = element_text(
                angle = 45
            ),
            #panel.grid.major = element_line(color = "grey70", size = 0.25),
            #panel.grid.minor = element_line(color = "grey90", size = 0.25),
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
            strip.text.x = element_text(hjust = 0)
        )
)



## == == == == == == == == == == == == == ==
## Step 4: Plot --------------------------------------------------------------------
## == == == == == == == == == == == == == ==

# Get unique combinations of 'algorithm' and 'I_alg'
alg_prog_df <- df %>% select(algorithm, I_alg) %>% unique()

# Get unique programs
prog_levels <- unique(df$I_alg)

# Assign colors to programs using the 'nejm' palette
nejm_colors <- pal_nejm("default")(length(prog_levels))
prog_colors <- setNames(nejm_colors, prog_levels)

# Map 'algorithm' to colors via 'I_alg'
alg_colors <- setNames(prog_colors[match(alg_prog_df$I_alg, names(prog_colors))], alg_prog_df$algorithm)


minor_breaks <- 2 ^ seq(-10, 10, by = 2)  # Every second exponent
major_breaks <- 2 ^ seq(-9, 9, by = 2)   # Exponents in between

ggplot(
    df,
    aes(
        bits_per_kmer,
        #bits_per_ref_kmer,
        Q_time_per_query_ms,
        color = algorithm,
        shape = algorithm,
        alpha = ifelse(qTypeLabel == "-", 0, 1)
    )
) +
    geom_point(size = ps) +
    facet_grid(variables) +
    scale_x_continuous(
        trans = 'log2',
        name = 'memory [bits / k-mer]',
        breaks = major_breaks,
        minor_breaks = minor_breaks,
        lim = c(2, NA)
    ) +
    scale_y_continuous(
        name = expression(paste('time per query [', mu, 's]')),
        trans = 'log2',
        breaks = major_breaks,
        minor_breaks = minor_breaks,
        labels = number_formatter
    ) +
    scale_shape_manual(values = my_shapes) +
    scale_color_manual(values = alg_colors) +
    guides(
        color = guide_legend("Group"),
        shape = guide_legend("Group"),
        alpha = 'none'
    )



## == == == == == == == == == == == == == ==
## Step 5: Save plot --------------------------------------------------------------------
## == == == == == == == == == == == == == ==

ggsave(fn,
       height = h,
       width = w,
       unit = u)
