#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

library(xtable)

desired_genome_order <- c("sc2-pg",
                          "sp-pg",
                          "ec-pg-hq",
                          "ec-pg-all",
                          "mtg-ilm",
                          "rna-ilm",
                          "hg-ilm",
                          "hg-t2t")


df0 <- read_tsv("final_results.tsv", na = "na")

# Correction for kilobytes in GNU time
df1 <- df0 %>%
    mutate(bits_per_kmer = bits_per_kmer * (1024 / 1000))

df2 <- df1 %>%
    mutate(bits_per_kmer_ind = 8 * index_bytes / kmer_count) %>%
    pivot_longer(
        cols = c("bits_per_kmer", "bits_per_kmer_ind"),
        names_to = "measurement",
        values_to = "bits_per_kmer"
    ) %>%
    mutate(measurement = ifelse(measurement == "bits_per_kmer_ind", "disk", "memory"))

df3 <- df2 %>%
    filter(qType == "Pos") %>%
    filter(plot == "isol") %>%
    filter(algorithm != "FMSI(spss)") %>%
    filter(15 <= k) %>%
    filter(k <= 43) %>%
    filter(genome %in% desired_genome_order) %>%
    mutate(genome = fct_relevel(genome, desired_genome_order)) %>%
    #filter(genome != "hg-t2t" | k <= 39) %>%
    #filter(genome != "sp-pg" | k <= 56) %>%
    mutate(algorithm = str_replace(algorithm, "FMSI.*", "FMSI")) %>%
    mutate(algorithm = str_replace(algorithm, "^BWA.*", "BWA")) %>%
    mutate(algorithm = str_replace(algorithm, "SSHash.*", "SSHash")) %>%
    filter(algorithm != "BWA") %>%
    mutate(genome = factor(genome, levels = desired_genome_order))


df4 <- df3 %>%
    select(c(
        "genome",
        "algorithm",
        "k",
        "measurement",
        "bits_per_kmer",
        "kmer_count"
    ))

df5 <- df4 %>%
    pivot_wider(names_from = "measurement", values_from = "bits_per_kmer")

df <- df5



################





convert_to_si <- Vectorize(function(x) {
    si_units <- c(
        "p" = -12,
        "n" = -9,
        "µ" = -6,
        "m" = -3,
        " " = 0,
        "k" = 3,
        "M" = 6,
        "G" = 9,
        "T" = 12
    )
    if (x == 0) {
        "0"
    } else{
        exponent <- floor(log10(abs(x))) # Find the order of magnitude
        closest_unit <- max(which(exponent >= unlist(si_units))) # Find the closest SI unit
        scaled_value <- x / (10 ^ (si_units[closest_unit]))
        formatted_value <- formatC(scaled_value, digits = 3, format = "g")
        result <- paste0(formatted_value, " ", names(si_units)[closest_unit])
        result
    }
})
sum0 <- df %>%
    group_by(genome, algorithm) %>%
    summarise(
        #kmer_count31 = format(kmer_count[k == 31], big.mark = ","),
        #kmer_count31M = format(round(kmer_count[k == 31] / 1e6), big.mark = ","),
        kmer_count31 = convert_to_si(kmer_count[k == 31]),
        ks = list(unique(k)),
        min_kmer_count = min(kmer_count),
        max_kmer_count = max(kmer_count),
        disk_min_bits_per_kmer = min(disk),
        disk_max_bits_per_kmer = max(disk),
        mem_min_bits_per_kmer = min(memory),
        mem_max_bits_per_kmer = max(memory),
        .groups = 'drop'  # This option removes the grouping structure afterwards
    )


sum1 <- sum0 %>%
    write_tsv("bits_per_kmer_range_summary.tsv")

rround <- function(x) {
    y <- round(10 * x) / 10
    sprintf("%.1f", y)
}



sum2 <- sum1 %>%
    bind_rows(
        sum1 %>%
            group_by(algorithm) %>%
            summarise(
                genome = "Overall",
                ks = NA,
                min_kmer_count = min(min_kmer_count),
                max_kmer_count = max(max_kmer_count),
                disk_min_bits_per_kmer = min(disk_min_bits_per_kmer),
                disk_max_bits_per_kmer = max(disk_max_bits_per_kmer),
                mem_min_bits_per_kmer = min(mem_min_bits_per_kmer),
                mem_max_bits_per_kmer = max(mem_max_bits_per_kmer),
                .groups = 'drop'
            )
    ) %>%
    mutate(across(ends_with("per_kmer"), rround))


getKDesc <- Vectorize(function(values) {
    if (is.null(values) || all(is.na(values))) {
        return("NA")
    }
    unique_values <- sort(unique(as.numeric(na.omit(values))))
    n <- length(unique_values)

    if (n < 3) {
        return("NA")
    } else {
        return(
            paste0(
                "\\KS{",
                unique_values[1],
                "}{",
                unique_values[2],
                "}{",
                unique_values[n],
                "}"
            )
        )
    }
}, vectorize.args = "values")

sum3 <- sum2 %>%
    transmute(
        dataset = genome,
        desc = getKDesc(ks),
        kmer_count31 = kmer_count31,
        algorithm,
        disk_ram = paste0(
            "\\RESULT{",
            mem_min_bits_per_kmer,
            "}{",
            mem_max_bits_per_kmer,
            "}{",
            disk_min_bits_per_kmer,
            "}{",
            disk_max_bits_per_kmer,
            "}"
        )
    ) %>%
    write_tsv("bits_per_kmer_range_summary_v2.tsv")


# correction for CBL
sum4 <- sum3 %>%
    mutate(desc = ifelse(
        dataset == "hg-t2t" & algorithm == "CBL",
        sum3 %>%
            filter(dataset == "hg-t2t" & algorithm=="FMSI") %>%
            select("desc"),
        desc
    ))

sum5 <- sum4 %>%
    mutate(desc = ifelse(dataset == "Overall", "", desc)) %>%
    pivot_wider(names_from = "algorithm", values_from = "disk_ram") %>%
    select("dataset",
           "desc",
           "kmer_count31",
           "CBL",
           "SSHash",
           "SBWT-fast",
           "SBWT-small",
           "FMSI") %>%
    rename(`k` = desc, `31-mers` = kmer_count31)

sum_final <- sum5


########

# Generate LaTeX table with \hline before 'Overall' rows

# Create xtable object
latex_table <- xtable(sum_final)

# Adjust column alignment as needed
align(latex_table) <- c("|l|", "|lr|r|", rep("r|", ncol(sum_final) - 1))

# Identify the positions where 'Overall' rows are located
# Assuming 'dataset' column starts with "Overall" for summary rows
overall_rows <- which(startsWith(as.character(sum_final$dataset), "Overall"))

# Create add.to.row list
if (length(overall_rows) > 0) {
    # Positions in add.to.row are 0-based and refer to rows after which to insert
    # To insert before 'Overall', use row index -1
    # Ensure that positions are valid (>=0)
    insert_positions <- overall_rows - 1
    insert_positions <- insert_positions[insert_positions >= 0]

    addtorow <- list()
    addtorow$pos <- as.list(insert_positions)
    addtorow$command <- rep("\\hline \n", length(insert_positions))
} else {
    # If no 'Overall' rows, no hline to add
    addtorow <- NULL
}

# Print the table without sink (optional, for console preview)
print(
    latex_table,
    include.rownames = FALSE,
    include.colnames = TRUE,
    floating = FALSE,
    sanitize.text.function = identity,
    add.to.row = addtorow
)

# Now, write to the .tex file
sink("bits_per_kmer_range_summary_v2.tex.tmp")

# Print the table with sink
print(
    latex_table,
    include.rownames = FALSE,
    include.colnames = TRUE,
    floating = FALSE,
    sanitize.text.function = identity,
    add.to.row = addtorow
)

sink()
