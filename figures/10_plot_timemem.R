#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

h <- 6
w <- 15 # 7.5 without the legend, 15 with legend
u <- "cm"

ps <- 2.0 # point size
fs <- 10 # font size



# == == == == == == == == == == == == == == == == ==
# Current plot combo configuration ----
#    ..adjust for current plot here...
# == == == == == == == == == == == == == == == == ==

# Small panel (not used)
#variables <- "genome ~ plot"; datasets <- c("Ecoli-pang"); ks <- c(23); h <- 6; w <- 15

# Fig 2, species x modes for k=23
#variables <- "genome ~ plot"; datasets <- c("ec-pg-all","ec-pg-hq","sc2-pg","ng-pg","sp-pg","mtg-ilm","rna-ilm","hg-t2t","hg-ilm","minikr4gib","minikr8gib"); ks <- c(23); h <- 12; w <- 12
#variables <- "genome ~ plot"; datasets <- c("ec-pg-all","mtg-ilm","hg-t2t"); ks <- c(23); h <- 12; w <- 12
variables <- "plot ~ genome"; datasets <- c("ec-pg-hq","mtg-ilm","hg-t2t"); ks <- c(23); h <- 12; w <- 12

# Sup. Fig 1; Spneumo: modes x k
#variables <- "ktext ~ plot"; datasets <- c("Spneumo-pang"); ks <- seq(from = 11, to = 67, by = 4); h <- 12; w <- 12
# variables <- "ktext ~ plot"; datasets <- c("Ngono-pang"); ks <- c(15, 23, 31); h <- 12; w <- 12
# variables <- "ktext ~ plot"; datasets <- c("Ecoli-pang-all"); ks <- c(15, 23, 31); h <- 12; w <- 12
# variables <- "ktext ~ plot"; datasets <- c("Ecoli-pang-HQ"); ks <- c(15, 23, 31); h <- 12; w <- 12
# variables <- "ktext ~ plot"; datasets <- c("SC2-pang"); ks <- c(15, 23, 31); h <- 12; w <- 12
# variables <- "ktext ~ plot"; datasets <- c("HG-T2T"); ks <- c(15, 23, 31); h <- 12; w <- 12
# variables <- "ktext ~ plot"; datasets <- c("minikraken4GB"); ks <- c(15, 23, 31); h <- 12; w <- 12


# new mode with many ks
#variables <- "ktext ~ plot"; datasets <- c("Spneumo-pang"); ks <- seq(from = 11, to = 67, by = 4); h <- 48; w <- 12
#variables <- "ktext ~ plot"; datasets <- c("HG-T2T"); ks <- seq(from = 11, to = 43, by = 4); h <- 36; w <- 12



#datasets <- c("Ecoli-pang-all", "Ecoli-pang-HQ", "Ngono-pang", "Spneumo-pang", "SC2-pang", "HG-T2T", "HG-illumina", "RNAseq", "microbiome", "minikraken4GB");
## for loop over datasets
#for (g in datasets)
#{

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

#fn <- paste0(basename, ".pdf")
fn <- "fig_memtime.pdf"

number_formatter <- function(y) {
    ifelse(y %% 1 == 0, as.character(as.integer(y)), as.character(y))
}

my_shapes <- c(15,  25, 2, 3, 10, 9, 19,  2)
my_colors <- c(1,  4, 4, 2, 3)
custom_palette <- pal_nejm("default")(max(my_colors))[my_colors]



## == == == == == == == == == == == == == ==
## Step 2: Load data ----
## == == == == == == == == == == == == == ==

df0 <- read_tsv("final_results.tsv", na = "na")

df1 <- df0 %>%
    filter(genome %in% datasets) %>%
    mutate(genome = factor(genome, levels = datasets)) %>%
    filter(k %in% ks) %>%
    filter(I_alg != "bwa") %>%
    filter(algorithm != "FMSI(spss)") %>%
    #mutate(plot2 = factor(plot, levels = c("stream", "isol", "subsampl"))) %>%
    mutate(
        plot = recode_factor(
            plot,
            "isol"             = "isol",
            "stream"           = "stream",
            "subsampl-isol"    = "subsampl-isol",
            "subsampl-stream"  = "subsampl-stream"
        )
    ) %>%
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
    mutate(ktext = paste0("k=", k))

df <- df1


## == == == == == == == == == == == == == ==
## Step 3: Set up theme --------------------------------------------------------------------
## == == == == == == == == == == == == == ==

theme_set(
    theme_bw(base_size = fs) +
        theme(
            axis.text.x = element_text(angle = 45),
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
    scale_color_manual(values = custom_palette) +
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

#}
