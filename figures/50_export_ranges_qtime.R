#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)
library(xtable)

# Desired genome order
desired_genome_order <- c(
    "sc2-pg",
    "sp-pg",
    "ec-pg-hq",
    "ec-pg-all",
    "mtg-ilm",
    "rna-ilm",
    "hg-ilm",
    "hg-t2t"
)

# 1) Read the input data
df0 <- read_tsv("final_results.tsv", na = "na")

# 2) Lightly clean/standardize the 'algorithm' field, factorize 'genome', etc.
df1 <- df0 %>%
    filter(15 <= k, k <= 43) %>%
    filter(genome %in% desired_genome_order) %>%
    mutate(
        genome    = fct_relevel(genome, desired_genome_order),
        algorithm = str_replace(algorithm, "FMSI.*", "FMSI"),
        algorithm = str_replace(algorithm, "^BWA.*", "BWA"),
        algorithm = str_replace(algorithm, "SSHash.*", "SSHash")
    ) %>%
    # Drop FMSI(spss) if desired:
    filter(algorithm != "FMSI(spss)")

# 3) Identify which rows are non‐streamed vs. streamed queries:
#    - Non‐streamed: typically (plot == "isol" AND qType == "Pos")
#    - Streamed:     typically (plot == "stream" AND qType == "Str")
#    Adjust if your data definitions differ.
df2 <- df1 %>%
    mutate(
        # Define a 'measurement' column to label which type of query time:
        measurement = case_when(
            plot == "isol"   & qType == "Pos" ~ "qtime_ns",  # non‐streamed
            plot == "stream" & qType == "Str" ~ "qtime_s",   # streamed
            TRUE                         ~ NA_character_
        )
    ) %>%
    # Keep only rows that correspond to either non‐streamed or streamed queries.
    filter(!is.na(measurement))

# 4) Pivot so that for each (genome, algorithm, k) we get columns qtime_ns and qtime_s.
#    Here, we use Q_time_per_query_ms (adjust if you prefer a different column).
df3 <- df2 %>%
    select(genome, algorithm, k, measurement, Q_time_per_query_ms, kmer_count) %>%
    pivot_wider(
        names_from  = measurement,
        values_from = Q_time_per_query_ms
    )

# 5) Summarize min/max query times across valid k for each genome/algorithm.
#    Also collect k‐mer counts, etc., as in your original script.
#    If desired, we can store "kmer_count" from the row where k==31,
#    or compute min/max over them, etc.
df4 <- df3 %>%
    group_by(genome, algorithm) %>%
    summarise(
        kmer_count31  = kmer_count[k == 31],  # pick out 31-mers if they exist
        ks            = list(unique(k)),
        min_kmer_count = min(kmer_count, na.rm = TRUE),
        max_kmer_count = max(kmer_count, na.rm = TRUE),
        # Non‐streamed times
        ns_min_qtime  = min(qtime_ns, na.rm = TRUE),
        ns_max_qtime  = max(qtime_ns, na.rm = TRUE),
        # Streamed times
        s_min_qtime   = min(qtime_s, na.rm = TRUE),
        s_max_qtime   = max(qtime_s, na.rm = TRUE),
        .groups = "drop"
    )

# Write out a quick TSV version
write_tsv(df4, "qtime_range_summary.tsv")

# 6) Add an overall row grouping by algorithm
df5 <- df4 %>%
    bind_rows(
        df4 %>%
            group_by(algorithm) %>%
            summarise(
                genome         = "Overall",
                kmer_count31   = NA,
                ks             = NA,
                min_kmer_count = min(min_kmer_count, na.rm = TRUE),
                max_kmer_count = max(max_kmer_count, na.rm = TRUE),
                ns_min_qtime   = min(ns_min_qtime, na.rm = TRUE),
                ns_max_qtime   = max(ns_max_qtime, na.rm = TRUE),
                s_min_qtime    = min(s_min_qtime, na.rm = TRUE),
                s_max_qtime    = max(s_max_qtime, na.rm = TRUE),
                .groups = "drop"
            )
    )

# If you like, round nicely:
rround <- function(x) {
    y <- round(10 * x) / 10
    sprintf("%.1f", y)
}

df6 <- df5 %>%
    mutate(
        ns_min_qtime = rround(ns_min_qtime),
        ns_max_qtime = rround(ns_max_qtime),
        s_min_qtime  = rround(s_min_qtime),
        s_max_qtime  = rround(s_max_qtime)
    )

# 7) Create a LaTeX-friendly version
#    (Mirroring your original approach, but adapting column names.)
#    If you want to replicate your \RESULT{} style macros, rename accordingly.
getKDesc <- Vectorize(function(values) {
    if (is.null(values) || all(is.na(values))) {
        return("NA")
    }
    unique_values <- sort(unique(as.numeric(na.omit(values))))
    n <- length(unique_values)
    if (n < 3) {
        return("NA")
    } else {
        return(paste0("\\KS{", unique_values[1], "}{", unique_values[2], "}{", unique_values[n], "}"))
    }
}, vectorize.args = "values")

df7 <- df6 %>%
    transmute(
        dataset      = genome,
        k           = getKDesc(ks),
        `31-mers`    = format(kmer_count31, big.mark = ",", scientific = FALSE),
        algorithm,
        # Build a column with the min/max for non-streamed vs. streamed:
        qtime_ranges = paste0(
            "\\RESULT{",
            ns_min_qtime, "}{", ns_max_qtime, "}{",
            s_min_qtime,  "}{", s_max_qtime,
            "}"
        )
    )

# Write out the second TSV
write_tsv(df7, "qtime_range_summary_v2.tsv")

# 8) Possibly do final pivot for table layout and correction for "CBL" if desired
#    (As in your original code).  Example:
df8 <- df7 %>%
    # Example correction for "hg-t2t & algorithm == 'CBL'", if needed:
    # mutate(k = ifelse(
    #   dataset == "hg-t2t" & algorithm == "CBL",
    #   <some other k-value logic>,
    #   k
    # )) %>%
    mutate(k = ifelse(dataset == "Overall", "", k)) %>%
    pivot_wider(
        names_from  = "algorithm",
        values_from = "qtime_ranges"
    ) %>%
    select(
        "dataset",
        "k",
        "31-mers",
        # reorder or rename columns as you prefer
        "CBL",
        "SSHash",
        "SBWT-fast",
        "SBWT-small",
        "FMSI"
    )

# 9) Make a LaTeX table
latex_table <- xtable(df8)
align(latex_table) <- c("|l|", "|lr|r|", rep("r|", ncol(df8) - 1))

# Insert a horizontal rule before rows whose 'dataset' is "Overall"
overall_rows <- which(startsWith(as.character(df8$dataset), "Overall"))
if (length(overall_rows) > 0) {
    insert_positions <- overall_rows - 1
    insert_positions <- insert_positions[insert_positions >= 0]

    addtorow <- list()
    addtorow$pos <- as.list(insert_positions)
    addtorow$command <- rep("\\hline \n", length(insert_positions))
} else {
    addtorow <- NULL
}

# Print to console (optional)
print(
    latex_table,
    include.rownames = FALSE,
    include.colnames = TRUE,
    floating = FALSE,
    sanitize.text.function = identity,
    add.to.row = addtorow
)

# 10) Write out the final .tex file
sink("qtime_range_summary_v2.tex.tmp")
print(
    latex_table,
    include.rownames = FALSE,
    include.colnames = TRUE,
    floating = FALSE,
    sanitize.text.function = identity,
    add.to.row = addtorow
)
sink()
