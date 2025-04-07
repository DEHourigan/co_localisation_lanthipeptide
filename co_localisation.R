########################################
## Packages - don't package shame people
########################################
.libPaths(c("/data/san/data0/users/david/rstudio/packages", .libPaths()))
newlib <- "/data/san/data0/users/david/rstudio/packages"

packages <- c("cowplot","data.table", "multidplyr", "readr", "formattable", "fs", "dplyr", "ggplot2", "purrr", "ggthemes", "BiocManager", "gplots", "gridExtra", "grid", "forcats", "tidyr", "dtplyr", "topGO", "SparseM", "Biostrings", "patchwork","GenomicRanges", "seqinr", "stringr", "readxl", "thacklr", "gggenomes")

load_packages <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org", 
      lib = newlib)
      library(package, character.only = TRUE)
    }
  }
}
load_packages(packages)

setwd("/data/san/data2/users/david/co_localisation")

########################################
## Functions used
########################################
writeFasta <- function(data, filename) {
  fastaLines <- c()
  for (rowNum in 1:nrow(data)) {
    fastaLines <- c(fastaLines, as.character(paste(">", data[rowNum, "name"], sep = "")))
    fastaLines <- c(fastaLines, as.character(data[rowNum, "seq"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


########################################
# function to append taxonomy given a nucleotide accession
########################################
library(taxonomizr)
# wget was used in terminal ## getAccession2taxid(baseUrl='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/')
# setwd("/data/san/data0/users/david/taxonomizer") #SET UP DATABASE ONCE #
sqlFile <- "/data/san/data0/users/david/taxonomizer/accessionTaxa.sql"
  
append_taxonomy <- function(df, sqlFile) {
  df$taxid <- taxonomizr::accessionToTaxa(df$subject.x, sqlFile)
  taxa <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  for(taxon in taxa) {
    df[[taxon]] <- getTaxonomy(df$taxid, sqlFile, desiredTaxa = taxon)
  }
  return(df)
}


setwd("/data/san/data2/users/david/co_localisation")
pfam_descs <- fread("tables/pfam_desc.tsv", col.names = c("pfam", "clan","clan_name","desc_pfam","name"))

########################################
# Load in rodeo table for class II lanthipeptide
########################################
LanM_40k_window = fread("data/processed/LanM_paper_40k_window.tsv")
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
LanM_phyla <- LanM_40k_window %>% 
  dplyr::select(phylum)
LanM_contigs <- LanM_40k_window %>%
  dplyr::select(Nucleotide_acc) %>%
  distinct()
fwrite(LanM_contigs, "data/processed/LanM_contigs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

########################################
# Load in PF00365_df
########################################
PF00365_df <- fread("data/negative_dataset/PF00365_df.tsv", sep = "\t") 

print_num_distinct <- function(df, filter_df, filter_column) {
  num_distinct_original <- df %>% dplyr::select(Nucleotide_acc) %>% n_distinct()
  num_distinct_filtered <- df %>% 
    filter(phylum %in% filter_df[[filter_column]]) %>%
    dplyr::select(Nucleotide_acc) %>%
    n_distinct()
  
  cat("\nNumber of distinct Nucleotide_acc in the original dataframe:", num_distinct_original)
  cat("\nNumber of distinct Nucleotide_acc in the filtered dataframe:", num_distinct_filtered, "\n")
}

print_num_distinct(PF00365_df, LanM_phyla, "phylum")

########################################
# create tables
########################################
# of note: Abi from paper Diep PF02517
create_table <- function(df, filter_df, filter_column) { # Define the function
  result_table <- df %>% 
    filter(phylum %in% filter_df[[filter_column]]) %>%
    dplyr::select(Nucleotide_acc, subject.y) %>% 
    distinct() 

  colnames(result_table) <- c("Nucleotide_acc", "pfam")
  result_table$group <- "control"
  result_table$pfam = sub("\\..*", "", result_table$pfam)

  return(result_table)
}

PF00365_table <- create_table(PF00365_df, LanM_phyla, "phylum")
PF00365_table$id <- "PF00365"

### only take controls with Pfams as annotated by PfamA db
PF00365_table <- PF00365_table %>% 
  group_by(Nucleotide_acc) %>%
  filter(any(pfam == "PF00365")) 

n_distinct(PF00365_table$Nucleotide_acc)
n_distinct(LanM_40k_window$Nucleotide_acc)

n_distinct(PF00365_table$Nucleotide_acc) + n_distinct(LanM_40k_window$Nucleotide_acc)

n_distinct(PF00365_table$pfam)
n_distinct(LanM_40k_window$pfam)
n_distinct(PF00365_table$pfam) + n_distinct(LanM_40k_window$pfam)

LanM_table <- LanM_40k_window %>% 
  dplyr::select(Nucleotide_acc, pfam) %>% # this is to ensure presence of pfam is counted once
  distinct()
LanM_table$group <- "lanII bacteriocin"
LanM_table$id <- "lan"

# find the pfams that are shared between two
shared_pfams <- unique(intersect(PF00365_table$pfam, LanM_table$pfam))


ctrl_v_lanm = rbind(PF00365_table, LanM_table) 
ctrl_df = rbind(PF00365_table)

write_tsv(ctrl_v_lanm, 
  "data/processed/ctrl_v_LanM_vs_Pks.tsv")
write_tsv(ctrl_df, 
  "data/processed/LanM_ctrl_df_pfams_nucs.tsv")

PF00365_table$Nucleotide_acc %>% unique()

########################################
# Plot what the control df looks like
########################################
library(ggsci)
library(ggplot2)
set.seed(1)
fig_1a <- PF00365_df %>%
filter(Nucleotide_acc %in% PF00365_table$Nucleotide_acc) %>%
  group_by(Nucleotide_acc) %>%
  dplyr::summarise(unique = n_distinct(locus_tag)) %>%
  ggplot(aes(x = "All Genomes", y = unique)) +
  geom_jitter(width = 0.3, size = 0.3, alpha = 0.3, color="black") +
  geom_violin(alpha = 0.7, fill="goldenrod", color = "black") +
  theme_bw() +
  labs(x = "", y = "Protein Count", title = "PF00365") +
   theme(
    plot.title = element_text(size = 12),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) 
# ggsave(fig_1a, filename = "figures/PF00365_control_window_protein_count_LanM.png", 
#     width = 6, height = 10, units = "cm")


########################################
# Plot what the LanM df looks like
########################################
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
fig_1b <- LanM_40k_window %>%
  group_by(unique_window) %>%
  summarise(unique = n_distinct(Protein_acc)) %>%
  ggplot(aes(x = "All Genomes", y = unique)) +
  geom_jitter(width = 0.3, size = 0.3, alpha = 0.5, color="black") +
  geom_violin(alpha = 0.7, fill="red", color="black") +
  theme_bw() +
  scale_color_npg() +
  scale_fill_npg() +
  labs(x = "", y = "Protein Count", title = "LanM") +
   theme(
    plot.title = element_text(size = 12),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)
  )
# ggsave(fig_1b, file="figures/LanM_bacteriocin_window_protein_count.png",
#   width = 6, height = 10, units = "cm")

fig_1d <- LanM_40k_window %>%
  dplyr::select(genus, Nucleotide_acc,start_window) %>% 
  mutate(genus = ifelse(is.na(genus), "Other", genus)) %>%
  distinct() %>%
  group_by(genus) %>%
  dplyr::summarise(count = n()) %>% 
  filter(count > 10) %>%
  arrange(desc(count)) %>%
  ggplot(., aes(x = reorder(genus, count), y = count)) +
  geom_bar(stat = "identity", color="black", fill="red", alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")

# ggsave(fig_1c, filename = "figures/LanM_paper_40k_window_genus_count.png",
#     height = 4, width = 4)

########################################
# Table Supplementary (plasmids)
########################################
# read in all plasmid files into a df from "data/platon_out"
plasmid_files <- list.files(path = "data/platon_out", pattern = "*.tsv", full.names = TRUE)
plasmid_df <- do.call(rbind, lapply(plasmid_files, fread, sep = "\t", header = TRUE)) %>%
  janitor::clean_names() %>%
  dplyr::rename(Nucleotide_acc = "id") 

plasmid_df <- plasmid_df %>%
  left_join(LanM_40k_window %>% dplyr::select(Nucleotide_acc, 'Genus/Species') %>% distinct())

plasmid_chrom_df <- LanM_40k_window %>%
  mutate(Type = ifelse(Nucleotide_acc %in% plasmid_df$Nucleotide_acc, "Plasmid", "")) %>%
  dplyr::select(Nucleotide_acc, 'Genus/Species', Type) %>%
  distinct()

# write the table
fwrite(plasmid_df, "tables/Supplementary_table_S3.csv")

fig_1c <- plasmid_chrom_df %>%
  group_by(Type) %>%
  dplyr::summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = "Total", y = percentage, fill = Type)) +
  geom_bar(stat = "identity", color = "black", aes(alpha=0.7)) +  
  theme_bw() +
  labs(x = "", y = "Percentage", title = "Location") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"  # Remove legend
  ) +
  geom_text(aes(label = Type), position = position_stack(vjust = 0.5), size = 3, color = "black") +
  scale_fill_manual(values = c("goldenrod", "red"))

########################################
#  Table Supplementary (GIs)
########################################
gi_files <- list.files(path = "data/gbff", pattern = "*.txt", full.names = TRUE)
gi_df <- data.table::rbindlist(lapply(gi_files, function(file) {
  df <- fread(file, sep = "\t", header = FALSE, fill = TRUE)
  df$Nucleotide_acc <- sub("\\.txt$", "", basename(file))
  df
}), fill = TRUE) %>%
dplyr::rename(genomic_island = V1,
  gi_start = V2,
  gi_end = V3) 

# find LanM that are present in GIs. i.e. start or end within gi_start and gi_end from gi_df
LanM_gi_df <- LanM_40k_window %>%
  left_join(gi_df, by = "Nucleotide_acc") %>%
  group_by(Nucleotide_acc,genomic_island ) %>%
  filter(start < gi_end & end > gi_start) %>%
  mutate(Type = ifelse(!is.na(genomic_island), "GI", "")) %>%
  dplyr::select(Nucleotide_acc, 'Genus/Species', Type, gi_start, gi_end) %>%
  distinct() %>%
  arrange(`Genus/Species`)

# write the table
fwrite(LanM_gi_df, "tables/Supplementary_table_GIs_S4.csv")

# plot the kernal densirty of GI lengths
fig_1e <- LanM_gi_df %>%
  filter(Type == "GI") %>%
  ggplot(aes(x = gi_length)) +
  geom_density(fill = "red", alpha = 0.7) +
  theme_bw() +
  labs(x = "Genomic Island Length", y = "Density") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )


########################################
# # Figure 1
########################################
f1_layout <- "ABDDD"

cowplot_figure_1 <- fig_1a + fig_1b  + fig_1d +
  plot_layout(design = f1_layout) +
  plot_annotation(tag_levels = 'A')


# send final figure to github repo
ggsave2("figures/figure_2.png", 
  cowplot_figure_1, width = 20, height = 10, units = "cm", dpi=600)
ggsave2("figures/figure_2.tiff", 
  cowplot_figure_1, width = 20, height = 10 units = "cm", dpi=600)


########################################
# Load in LanM vs control pfam dataset as a checkpoint
########################################
ctrl_v_lanm = fread("data/processed/ctrl_v_LanM_vs_Pks.tsv")


########################################
# Create contingency table for lanM
########################################
create_contingency_table <- function(x, y, control_name) {
  LanM_table <- y %>% 
    dplyr::select(Nucleotide_acc, pfam) %>% 
    distinct()
  LanM_table$group <- "lanthipeptide II"

  x <- x %>%
    dplyr::select(Nucleotide_acc, pfam) %>%
    distinct()
  x$group <- control_name

  # make contingency table for pfams
  k <- rbind(LanM_table, x) %>% 
    distinct() 
  # Filter out rows with complete cases
  k <- k[complete.cases(k), ]  
  k[is.na(k)] <- "Unknown"

  contingency_table <- table(k$group, k$pfam)
  total_df <- sum(contingency_table[1, ])
  total_control <- sum(contingency_table[2, ])
  
  cat("\nTotal for lanII bacteriocin:", total_df)
  cat("\nTotal for", control_name, ":", total_control, "\n")
  
  return(contingency_table)
}


cont_control_PF00365 <- create_contingency_table(PF00365_table, LanM_40k_window, "PF00365")


########################################
# Stats on the LanM vs control pfam dataset
########################################

LanMin = LanM_40k_window %>%
	dplyr::select(Nucleotide_acc, pfam) %>%
	mutate(group = "lanthipeptide",
		id = "lanthipeptide II") %>%
	filter(grepl("^PF", pfam))

tablex = rbind(LanMin, PF00365_table) %>% dplyr::select(-id)
fwrite(tablex, "data/processed/lanthipeptide_pfam_table.tsv")

# check point
tablex = fread("data/processed/lanthipeptide_pfam_table.tsv")


# Summarize counts per PFAM per group
pfam_counts <- tablex %>%
  group_by(pfam, group) %>%
  summarize(count = n(), .groups = "drop")

# Calculate mean and 2*SD for each group
group_stats <- pfam_counts %>%
  group_by(group) %>%
  summarize(mean_count = mean(count), sd_count = sd(count), .groups = "drop") %>%
  mutate(upper_limit = mean_count + 2 * sd_count) %>%
  ungroup()

# Plot histogram of PFAM counts, faceted by group, with 2*SD lines
ggplot(pfam_counts, aes(x = count, fill = group)) +
  geom_histogram(bins = 16, alpha = 0.6, position = "identity", color = "black") +
  facet_wrap(~group, scales = "free_x") +
  labs(x = "PFAM Count", y = "Frequency", title = "Histogram of PFAM Counts by Group") +
  theme_minimal() +
  scale_x_log10() +  # Log scale if needed for better visualization
  geom_vline(data = group_stats, aes(xintercept = upper_limit, color = group), linetype = "dashed", size = 0.8) +
  scale_color_manual(values = c("red", "blue"))  # Adjust colors as needed




plot_size_of_pfam_pools_df = tablex %>%
  group_by(group) %>%
  summarize(total = n())

n_distinct(tablex$Nucleotide_acc)
n_distinct(tablex$pfam)

# Step 1: Calculate total contigs per group
total_group_contigs <- tablex %>%
  distinct(Nucleotide_acc, group) %>%
  group_by(group) %>%
  summarise(total_group_contigs = n_distinct(Nucleotide_acc))


tablex %>% 
  filter(group == "lanthipeptide") %>%
  group_by(pfam) %>%
  summarise(total_lanthipeptide_contigs = n_distinct(Nucleotide_acc))

tablex %>%
  filter(group == "control") %>%
  distinct(Nucleotide_acc)

result <- tablex %>%
  distinct(Nucleotide_acc, pfam, group) %>%
  group_by(group, pfam) %>%
  summarise(contigs_with_pfam = n_distinct(Nucleotide_acc)) %>%
  left_join(total_group_contigs, by = "group") %>%  # Join the total contig counts for each group
  mutate(percent_contigs = (contigs_with_pfam / total_group_contigs) * 100) %>% # Calculate percentage
  arrange(group, desc(percent_contigs))

percent_test = result %>% 
  arrange(desc(percent_contigs)) %>% 
  filter(group == "lanthipeptide") %>%
  left_join(., pfam_descs) %>%
  dplyr::select(pfam, clan, name, percent_contigs, contigs_with_pfam, total_group_contigs)

percent_control = result %>% 
  arrange(desc(percent_contigs)) %>%
  filter(group == "control") %>%
  left_join(., pfam_descs) %>%
  dplyr::select(pfam, clan, name, percent_contigs, contigs_with_pfam, total_group_contigs)

# Merge the dataframes and count the fold difference between pfams per group
merged_df <- merge(percent_test, percent_control, by = c("pfam", "clan", "name"), suffixes = c("_lanthipeptide", "_control"))

# Calculate the fold difference
merged_df2 <- merged_df %>%
  mutate(fold_difference = percent_contigs_lanthipeptide / percent_contigs_control) %>%
  dplyr::select(-group_lanthipeptide, -group_control) %>% 
  arrange(desc(fold_difference))

fwrite(merged_df2, "fold_difference_table.tsv")



# Calculate total counts and counts per group
counts <- tablex %>%
  group_by(pfam) %>%
  summarize(
    total_genes = n(),
    genes_in_lanthipeptide = sum(group == "lanthipeptide"),
    genes_in_control = sum(group == "control")
  ) %>%
  filter(pfam != "") %>%
  arrange(desc(total_genes))

# Calculate fractions for each group
N <- sum(counts$total_genes)
n_lanthipeptide <- sum(counts$genes_in_lanthipeptide)
n_control <- sum(counts$genes_in_control)
f_lanthipeptide <- n_lanthipeptide / N
f_control <- n_control / N

# Calculate expected numbers and standard deviations for each group
counts_2 <- counts %>%
  mutate(
    m_prime_lanthipeptide = total_genes * f_lanthipeptide,
    s_prime_lanthipeptide = sqrt(total_genes * f_lanthipeptide * (1 - f_lanthipeptide)),
    m_prime_control = total_genes * f_control,
    s_prime_control = sqrt(total_genes * f_control * (1 - f_control))
  )


significant_threshold <- 0.05

counts_3 <- counts_2 %>%
  mutate(
    Z_score_lanthipeptide = ifelse(s_prime_lanthipeptide > 0, 
                                   (genes_in_lanthipeptide - m_prime_lanthipeptide) / s_prime_lanthipeptide, 
                                   NA),
    Z_score_control = ifelse(s_prime_control > 0, 
                             (genes_in_control - m_prime_control) / s_prime_control, 
                             NA),
    p_value = pnorm(-abs(Z_score_lanthipeptide)),
    bonferroni_p_value = pmin(1, p_value * n()),
    p_value_control = pnorm(-abs(Z_score_control)),
    bonferroni_p_value_control = pmin(1, p_value_control * n()),
    
    # Check significance after Bonferroni correction
    significant = bonferroni_p_value < significant_threshold,

    log_Z_score_lanthipeptide = ifelse(Z_score_lanthipeptide >= 0, log1p(Z_score_lanthipeptide), NA),
    log_Z_score_control = ifelse(Z_score_control >= 0, log1p(Z_score_control), NA),
  ) %>%
  filter(!is.na(pfam)) %>%
  left_join(pfam_descs, by="pfam")
  

# View(counts_3)
fwrite(counts_3, "tables/Supplementary_table_S1.csv")


counts_3$Z_score_lanthipeptide <- as.numeric(counts_3$Z_score_lanthipeptide)


library(ggrepel)
counts_3_color <- counts_3 %>%
    mutate(color_group = ifelse(log_Z_score_lanthipeptide > 1.64, ">1.64 z-score in Lanthipeptide", "<1.64 z-score in Lanthipeptide")) 

# Define the vector of pfams to be labeled
labelled_points <- c("PF00005", "PF05147", "PF13575", "PF03412","PF00365","PF16934",
  "PF07733", "PF00224")


# Scatter plot Lanthipeptide
scatter_plot <- ggplot(counts_3_color %>% filter(!is.na(log_Z_score_lanthipeptide)), 
             aes(x = total_genes, y = log_Z_score_lanthipeptide, fill = color_group)) +
  geom_point(
    alpha = 0.3,
    shape=21,
    color = "black") +
  labs(title = "Scatter Plot: Log Z-scores vs. Total Genes", x = "Total Genes", y = "Log Z-Score (Lanthipeptide)") +
  theme_classic(base_size = 8) +
  scale_fill_manual(values = c(">1.64 z-score in Lanthipeptide" = "red", "<1.64 z-score in Lanthipeptide" = "grey")) +
  theme(legend.title = element_blank(),
    legend.position = c(0.65, 0.20)) + # Move legend to bottom right corner
  geom_text_repel(data = counts_3_color %>% filter(pfam %in% labelled_points & !is.na(log_Z_score_lanthipeptide)),
          aes(label = pfam),
          vjust = -0.5,
          hjust = 1,
          size = 2.5,
          color = "black")

  # ggsave(filename = "figures/scatter_plot_colored.png", 
  #        plot = scatter_plot, width = 10, height = 7, units = "cm",
  #        dpi=600)

labelled_points_control = c(
  "PF07733", "PF00224","PF00365","PF03328",
  # "PF07733", "PF14579", "PF01336",
  # "PF02780","PF02129","PF00391","PF10776"
  "PF01039","PF00056"
  )

# Scatter plot Control
scatter_plot_control <- ggplot(counts_3_color %>% filter(!is.na(log_Z_score_control)), 
  aes(
    x = total_genes, 
    y = log_Z_score_control, 
    fill = factor(log_Z_score_control > 1.64, levels = c(FALSE, TRUE), labels = c("<1.64 z-score in Control", ">1.64 z-score in Control")))) +
    geom_point(
      alpha = 0.3,
      shape=21,
      color = "black") +
    labs(title = "Scatter Plot: Log Z-scores vs. Total Genes", x = "Total Genes", y = "Log Z-Score (Control)") +
    theme_classic(base_size = 8) +
    scale_fill_manual(values = c(">1.64 z-score in Control" = "goldenrod", "<1.64 z-score in Control" = "grey")) +
    theme(legend.title = element_blank(),
      legend.position = c(0.65, 0.20))  +
    geom_text_repel(data = counts_3_color %>% filter(pfam %in% labelled_points_control) %>% filter(!is.na(log_Z_score_control)),
                    aes(label = pfam),
                    vjust = -0.5,
                    hjust = 1,
                    size = 2.5,
                    color = "black")

counts_3_color %>%
  filter(pfam %in% c("PF00391", "PF10776","PF01420"))

########################################
# Plot expected vs observed counts
########################################
expected_vs_observed_plot <- ggplot(counts_3) +
  geom_point(aes(x = m_prime_control, y = m_prime_control, fill = "Expected"), 
             alpha = 0.3, shape = 21, color="black") +
  geom_point(aes(x = m_prime_control, y = genes_in_control, fill = "Control"), 
             shape = 21, alpha = 0.3, color="black") +
  geom_point(aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, fill = "Lanthipeptide"), 
             shape = 21, color="black") +
  scale_fill_manual("", values = c("Expected" = "darkgrey", "Control" = "goldenrod", "Lanthipeptide" = "red")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  coord_cartesian(clip = "off") +
  labs(x = "Expected Number of Pfams", y = "Observed Number of Pfams", 
    title = "Expected vs. Observed") +
  theme_classic(base_size = 8) +
  geom_label_repel(data = subset(counts_3, (name == "ABC transporter" & significant == TRUE) | name == "Lanthionine synthetase C-like protein"),
                   aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = name),
                   color = "black", size = 2, segment.color = "grey88") +
  geom_label_repel(data = subset(counts_3, (name == "ABC transporter" & significant == TRUE) | name == "Lanthionine synthetase C-like protein"), 
                   aes(x = m_prime_control, y = genes_in_control, label = name),
                   color = "black", size = 2, segment.color = "grey88") +
  theme(legend.position = c(0.18, 0.94), 
        legend.background = element_blank(),
        legend.spacing = unit(0.1, 'cm'),
        legend.key.size = unit(0.3, 'cm'))

# ggsave(filename = "figures/expected_vs_observed_plot.png", 
#   plot = expected_vs_observed_plot, width = 10, height = 7, units = "cm", dpi=600)


subset_count = counts_3 %>%
  filter(m_prime_lanthipeptide < 500 &  m_prime_control < 500) %>%
  filter(genes_in_lanthipeptide > m_prime_lanthipeptide) %>%
  filter(bonferroni_p_value < 0.05)

# Vector of pfam labels to highlight
pfam_labels_2 <- c(
  "Competence-damaged protein", 
  "Methylase_S",
  "Methyltransf_25", 
  "Mersacidin", 
  "Lantibiotic_a",
  "Peptidase_S8",
  "AbiJ_NTD3", 
  "MqsA_antitoxin", 
  "Phage_integrase",
  "EcoR124_C",
  "rve",
  "DEAD",
  "CinA",
  "RsgA_GTPase",
  "Acetyltransf_1",
  "Acetyltransf_10",
  "Acetyltransf_3",
  "Acetyltransf_6",
  "Acetyltransf_7",
  "Eco57I")

expected_vs_observed_plot_subset <- ggplot(subset_count) +
  geom_point(aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide),
    shape = 21,
    color = "black",
    fill = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Expected Number of Pfams", y = "Observed Number of Pfams", 
       title = "Expected vs. Observed (< 500 Expected)") +
  theme_classic(base_size = 8) +
  geom_label_repel(data = subset(subset_count, desc_pfam %in% pfam_labels_2), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black",
             size = 2,
             segment.color = "grey88",
             box.padding = 0.35,
             point.padding = 0.5,
             max.overlaps = Inf,
             force = 100)

cowplot_figure_3 = cowplot::plot_grid(
  scatter_plot_control, scatter_plot,
  expected_vs_observed_plot, expected_vs_observed_plot_subset,
  labels = c("A", "B", "C", "D"),
  ncol = 2,
  label_size = 10)

ggsave(filename = "figures/figure_3.png", 
  plot = cowplot_figure_3, width = 16, height = 12, units = "cm", dpi=600)
ggsave(filename = "figures/figure_3.tiff", 
  plot = cowplot_figure_3, width = 16, height = 12, units = "cm", dpi=600)

########################################
# chi-squared 
########################################
# null hypothesis is that there is no correlation between the two groups
contingency_table <- table(tablex$group, tablex$pfam)
chi_results_mc <- chisq.test(contingency_table, 
  simulate.p.value = TRUE, 
  B = 2000)
print(chi_results_mc) # reject the null hypothesis meaning there is a correlation


########################################
# Print Cramér's V
########################################
cramers_v <- sqrt(chi_results_mc$statistic / (sum(contingency_table) * (min(dim(contingency_table)) - 1)))
print(cramers_v)
# this value of >0.5 means there is a association between group and Pfam)
# confirm method with package
# install.packages("lsr")
library(lsr)
cramersV(contingency_table)


########################################
# Exploration of dataset
########################################
def = fread("tables/phage_defense.txt", header=FALSE, 
  col.names = c("pfam"))
pfams_proteins = LanM_40k_window %>% 
  dplyr::select(-Query) %>%
  distinct()

# what are the highest pfams that are part of defence systems?
pfams_highest = pfams_proteins %>%
  filter(pfam %in% def$pfam) %>%
  dplyr::select(pfam,  desc,Protein_acc) %>%
  group_by(pfam, desc) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
fwrite(pfams_highest, "tables/pfams_highest.tsv")

def_red = tail(pfams_highest, 104)

pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(Nucleotide_acc, genus,species,"Genus/Species", unique_window) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>% n_distinct()


pfams_def_genus_counts = pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(pfam, desc, genus) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
fwrite(pfams_def_genus_counts, "tables/pfams_def_genus_counts.tsv")

#### what have more 2 or more pfams from defence
pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(Nucleotide_acc, genus,species,"Genus/Species") %>%
  summarize(n = n()) %>%
  filter(n > 4) %>%
  arrange(desc(n))  %>%
  distinct(Nucleotide_acc) 

pfams_proteins %>%
  filter(pfam %in% def_red$pfam) %>%
  distinct() %>%
  group_by(genus) %>%
  summarize(n = n()) %>%
  filter(n > 4) %>%
  distinct(genus) 

coloc_df = pfams_proteins %>% 
  dplyr::select(Nucleotide_acc, genus, species, Protein_acc, pfam) %>% 
  distinct()

coloc_df %>%
  filter(pfam %in% def_red$pfam) %>%
  # merge pfams if they on the same protein
  group_by(Nucleotide_acc, Protein_acc) %>%
  summarise(pfams = paste(pfam, collapse = ";")) %>%
  group_by(Nucleotide_acc) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


#### Bacteriocin
bacteriocins = pfams_proteins %>% 
  filter(pfam == "PF10439") %>%
  dplyr::select(Protein_acc) %>%
  distinct()

fwrite(bacteriocins, "tables/bacteriocins.tsv")

#### COMPETENCE
competence = c("PF07508","PF07508","PF05952","PF00154")
pfams_proteins %>%
  filter(pfam %in% competence) %>%
  distinct() %>%
  group_by(Nucleotide_acc, genus,species,"Genus/Species", unique_window) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>% n_distinct()

pfams_competence_genus_counts = pfams_proteins %>%
  filter(pfam %in% competence) %>%
  distinct() %>%
  group_by(pfam, desc, genus) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
fwrite(pfams_competence_genus_counts, "tables/pfams_competence_genus_counts.tsv")

sum(is.na(tablex))
# Check unique values
length(unique(tablex$Nucleotide_acc))
length(unique(tablex$pfam))
length(unique(tablex$group))

# View structure of the dataframe
str(tablex)

# Frequency counts
table(tablex$pfam)
table(tablex$group)

# Cross-tabulation
table(tablex$group, tablex$pfam)

# Cross-tabulation
tbl <- table(tablex$group, tablex$pfam)

# Chi-Square Test of Independence with Monte Carlo simulation
chi_results_mc <- chisq.test(tbl, simulate.p.value = TRUE, B = 2000)
print(chi_results_mc)
print("These counts are not independent and reject the null hypotheis")

# Using expected frequencies from the Monte Carlo Chi-Square test
expected_mc <- chi_results_mc$expected

# Printing expected frequencies
print(expected_mc)

# Over-representation check (observed > expected)
over_represented <- tbl > expected_mc
print(over_represented)

# Reshape the data for visualization and further analysis
melted_data <- reshape2::melt(tbl)
colnames(melted_data) <- c("group", "pfam", "count")
melted_expected <- reshape2::melt(expected_mc)
colnames(melted_expected) <- c("group", "pfam", "expected")

# Merging datasets
comparison <- merge(melted_data, melted_expected, by=c("group", "pfam"))
comparison <- comparison %>%
  mutate(ratio = count / expected)

# remove pfams with NA
comparison <- comparison %>%
	left_join(., pfam_descs) %>%
	filter(!is.na(desc)) 

# Plotting the ratio of observed to expected counts, with colors by group
ratio_plot_observed_vs_expected <- ggplot(comparison, aes(x = reorder(pfam, -ratio), y = ratio, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("goldenrod", "red")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x.top = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.title.x = element_text()
  ) +
  xlab("Pfam") +
  labs(title = "Ratio of Observed to Expected Counts")

ggsave(
  filename = "figures/ratio_plot_observed_vs_expected.png",
  plot = ratio_plot_observed_vs_expected,
  width = 15, height = 10, units = "cm"
)

comparison_filter = comparison %>% filter(pfam %in% PF00365_fishers_LanM_rodeo$pfam)

comparison %>% 
	filter(ratio > 1) %>%
	filter(group == "lanthipeptide") %>%
	arrange(desc(ratio)) %>% 
	filter(pfam %in% c("PF05016", "PF01420","PF05147"))


########################################
# Fisher's exact test
########################################
perform_test <- function(df, control_group) {
  
  # group counts
  total_circular_bacteriocin <- sum(df["lanthipeptide II",])
  total_control <- sum(df[control_group,])
  
  test_func <- function(column) {
    present_counts <- as.numeric(column)
    absent_counts <- c(total_circular_bacteriocin, total_control) - present_counts
    variable_table <- rbind(present_counts, absent_counts)
  
    # Fisher's test
    test_result <- fisher.test(variable_table)
  
    # Proportion difference
    proportion_diff <- present_counts[1]/total_circular_bacteriocin - present_counts[2]/total_control
  
    return(data.frame(variable = names(column), p_value = test_result$p.value, proportion_diff = proportion_diff))
  }
  results <- apply(df, 2, test_func)
  return(results)
}

# Running the function on each dataframe
results_PF00365_fishers <- perform_test(cont_control_PF00365, "PF00365")
results_PF00365_fishers_df <- do.call(rbind, lapply(results_PF00365_fishers, function(x) {
  x$control <- "PF00365"
  x
}))



all_fishers_df <- results_PF00365_fishers_df %>%
  tibble::rownames_to_column(var = "pfam")
 
all_fishers_df_clean = all_fishers_df %>% 
  dplyr::select(pfam, p_value, proportion_diff, control) %>% 
  distinct()

all_fishers_df_clean$pfam <- gsub("\\.\\d+$", "", all_fishers_df_clean$pfam)
all_fishers_df_clean = all_fishers_df_clean %>%
  distinct()

all_fishers_df_clean$adjusted_p_value <- p.adjust(all_fishers_df_clean$p_value, 
  method = "bonferroni") # Apply Bonferroni correction

pfam_descs_2 = LanM_40k_window %>% dplyr::select(pfam, desc) %>% distinct() 

PF00365_fishers_LanM_rodeo = left_join(all_fishers_df_clean, pfam_descs_2) %>%
  distinct() %>% 
  arrange(adjusted_p_value) %>%
  filter(grepl("^PF", pfam)) %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(proportion_diff > 0) %>%
  filter(control == "PF00365")

# fwrite(PF00365_fishers_LanM_rodeo, "tables/PF00365_fishers_LanM_rodeo.tsv", sep = "\t")
fwrite(PF00365_fishers_LanM_rodeo, "tables/Supplementary_table_S2.tsv", sep = "\t")



########################################
# Volcano Plot
########################################
library(ggrepel)
vplot1 = PF00365_fishers_LanM_rodeo %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(proportion_diff > 0) %>%
  left_join(dplyr::select(LanM_40k_window, desc, pfam) %>% distinct()) %>%
  arrange(desc(adjusted_p_value)) %>%
  filter(control == "PF00365") %>%
  ggplot(aes(x = log10(proportion_diff), y = -log10(adjusted_p_value))) +
  geom_jitter(
    alpha = 0.9,
    size = 1,
    color = "black",
    fill = "red",
    shape = 21
  ) +
  geom_label_repel(
    data = . %>% filter(desc %in% c(
      "Lanthionine synthetase C-like protein",
      "Type-A lantibiotic",
      "Acetyltransferase (GNAT) family",
      "Acetyltransferase (GNAT) domain",
      "Phage_integrase",
      "RsgA GTPase",
      "Recombinase",
      "Bacterial mobilisation protein (MobC)",
      "RelB antitoxin",
      "Resolvase, N terminal domain",
      "Integrase core domain",
      "Competence-damaged protein",
      "ParB/Sulfiredoxin domain",
      "Bacteriophage CI repressor helix-turn-helix domain",
      "Recombinase zinc beta ribbon domain",
      "Bacillus competence pheromone ComX",
      "MobC-like protein",
      "ABC transporter"
    )),
    aes(label = desc),
    vjust = 1,
    hjust = 1,
    size = 2
  ) +
  xlab("log10(Proportion difference)") +
  ylab("-log10(Adjusted P-Value)") +
  ggtitle("Volcano Plot: LanM-associated Pfams") +
  scale_color_manual(values = c("black", "#a22d2d")) +
  theme_classic(base_size = 8, base_family = "arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  coord_cartesian(clip = "off")

vplot2 = PF00365_fishers_LanM_rodeo %>% 
  filter(proportion_diff < 0.0015) %>%
  filter(adjusted_p_value < 0.05) %>% 
  filter(proportion_diff > 0) %>%
  ggplot(aes(x = log10(proportion_diff), y = -log10(adjusted_p_value))) +
  geom_jitter(alpha = 0.9, 
    shape=21,
    size = 1, 
    aes(color="black",
      fill= "red")) +
  geom_label_repel(
    data = . %>% filter(desc %in% c(
      "Lanthionine synthetase C-like protein",
      "Type-A lantibiotic",
      "Acetyltransferase (GNAT) family",
      "Acetyltransferase (GNAT) domain",
      "Phage_integrase",
      "RsgA GTPase",
      "Recombinase",
      "Bacterial mobilisation protein (MobC)",
      "RelB antitoxin",
      "Resolvase, N terminal domain",
      "Integrase core domain",
      "Competence-damaged protein",
      "ParB/Sulfiredoxin domain",
      "Bacteriophage CI repressor helix-turn-helix domain",
      "Recombinase zinc beta ribbon domain",
      "Bacillus competence pheromone ComX",
      "MobC-like protein",
      "ABC transporter"
    )),
      max.overlaps = Inf,
      min.segment.length = 0,
      arrow = arrow(length = unit(0.015, "npc")),
      aes(label = desc), 
        color="black",
        fill="white",
        size = 2, 
           ) + 
  xlab("log10(Proportion difference)") +
  ylab("-log10(Adjusted P-Value)") +
  ggtitle("Volcano Plot: LanM-associated Pfams") +
  scale_color_manual(values = c("black", "#a22d2d")) +
  theme_classic(base_size = 8, base_family = "arial") +
  scale_fill_identity() +
  theme(legend.position = "none",
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) +
  coord_cartesian(clip="off")



########################################
# PFAM TO GO for enriched
########################################
library(ragp)
library(topGO)
library(GO.db)
PF00365_fishers_LanM_rodeo = fread("tables/Supplementary_table_S2.tsv")
top_n = 10
DF_GO <- pfam2go(data_pfam = PF00365_fishers_LanM_rodeo, pfam = "pfam")
fwrite(DF_GO, "tables/Supplementary_table_GO.tsv", sep = "\t")

unique(DF_GO$GO_name)

hgt_GO = c("GO:DNA integration",
 "GO:DNA recombination",
 "GO:DNA replication initiation",
 "GO:DNA strand exchange activity",
 "GO:DNA transposition",
 "GO:transposase activity",
 "GO:acyltransferase activity, transferring groups other than amino-acyl groups",
 "GO:double-strand break repair",
 "GO:endonuclease activity",
 "GO:DNA binding",
 "GO:5-methyltetrahydropteroyltriglutamate-homocysteine S-methyltransferase activity",
 "GO:methyltransferase activity",
 "GO:methyltransferase activity",
 "GO:transposase activity",
 "GO:viral process",
 "GO:DNA integration",
 "GO:double-strand break repair",
 "GO:negative regulation of DNA-templated transcription"
 )

top_GO_terms_hgt <- DF_GO %>%
  filter(!is.na(GO_name)) %>%
  filter(GO_name != "hypothetical") %>%
  filter(GO_name != "") %>%
  group_by(GO_name) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(desc(total)) %>% 
  filter(GO_name %in% c(hgt_GO)) %>%
  top_n(top_n, total)

top_GO_terms <- DF_GO %>%
  filter(!is.na(GO_name)) %>%
  filter(GO_name != "hypothetical") %>%
  filter(GO_name != "") %>%
  group_by(GO_name) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(desc(total)) %>% 
  top_n(top_n, total)

# merge the two dataframes
top_GO_terms <- rbind(top_GO_terms, top_GO_terms_hgt)

top_GO_terms_df <- DF_GO %>%
  filter(!is.na(GO_name)) %>%
  filter(GO_name != "hypothetical") %>%
  filter(GO_name != "") %>%
  group_by(GO_name) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(desc(total)) 
fwrite(top_GO_terms_df, "tables/Supplementary_table_GO_terms.tsv", sep = "\t")

########################################
# GO (1) Biological Process
########################################

# Plot
top_GO_terms$font_face <- ifelse(top_GO_terms$GO_name %in% hgt_GO, "bold", "plain")

GO_plot <- ggplot(top_GO_terms_hgt, 
  aes(y = reorder(GO_name, total), x = total)) +
  geom_bar(stat = "identity", position = "dodge",
    color = "black", fill="red") +
  labs(y = "", x = "GO terms of LanM-associated Pfams") +
  theme_classic(base_size = 8, base_family = "arial") +
  scale_fill_manual(values = c("lan" = "#d02727", "PF00365" = "goldenrod", "PF0090" = "#686a68", "PF01225" = "#4c9f7b")) +
  ggtitle(paste("Top", top_n, "GO: Biological Process")) +
  theme(axis.text.y = element_text(face = "bold"))

ggsave("figures/LanM_GO_BP_v_pks.png", 
    GO_plot, width = 14, height = 7, units = "cm")

ggsave("figures/LanM_GO_BP_v_pks.png", 
  GO_plot, width = 14, height = 7, units = "cm")


cowplot_figure_4 = cowplot::plot_grid(
  cowplot::plot_grid(vplot1, vplot2, labels = c("A", "B"), ncol = 2, label_size = 10),
  GO_plot,
  labels = c("", "C"),
  ncol = 1,
  rel_heights = c(1,0.8),
  label_size = 10
)

# save figure 4
ggsave2("figures/figure_4.png", 
  cowplot_figure_4, width = 20, height = 16, units = "cm", dpi=600)
ggsave2("figures/figure_4.tiff", 
  cowplot_figure_4, width = 20, height = 16, units = "cm", dpi=600)



########################################
# Look for patterns of systems    1) PHAGE DEFENSE SYSTEMS
########################################
ctrl_v_lanm = fread("data/processed/ctrl_v_LanM_vs_Pks.tsv")

library("readxl")
rsphagepaper1 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 1)
rsphagepaper2 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 2)
rsphagepaper3 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 3)
rsphagepaper4 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 4)
rsphagepaper5 <- read_excel("/data/san/data0/users/david/intelligence/rotem_sorek/aar4120_tabless1-s5.xlsx", sheet = 5)

# Create a vector with the 5 objects
list_of_tibbles  <- c(rsphagepaper1, rsphagepaper2, rsphagepaper3, rsphagepaper4, rsphagepaper5)

phage_defense_pfams1 =  rsphagepaper1 %>% unique() %>% filter(str_starts(Family, "pfam")) %>% 
    dplyr::select(Family) %>% distinct() %>% arrange(Family) %>%
    distinct() %>%
    mutate(Family = str_replace_all(Family, "pfam", "PF"))

phage_defense_pfams2 =  rsphagepaper2 %>% unique() %>% filter(str_starts(Family, "pfam")) %>% 
    dplyr::select(Family) %>% distinct() %>% arrange(Family) %>%
    distinct() %>% 
    mutate(Family = str_replace_all(Family, "pfam", "PF"))

phage_defense_pfams3 =  rsphagepaper3 %>% unique() %>% filter(str_starts(`Anchor protein families`, "pfam")) %>% 
    dplyr::select(`Anchor protein families`) %>% distinct() %>% arrange(`Anchor protein families`) %>%
    distinct() %>% 
    mutate(Family = str_replace_all(`Anchor protein families`, "pfam", "PF")) %>% dplyr::select(Family)

phage_defense_pfams4 =  rsphagepaper4 %>% unique() %>% filter(str_starts(`Anchor protein families`, "pfam")) %>% 
    dplyr::select(`Anchor protein families`) %>% distinct() %>% arrange(`Anchor protein families`) %>%
    distinct() %>% 
    mutate(Family = str_replace_all(`Anchor protein families`, "pfam", "PF")) %>% dplyr::select(Family) %>%
    separate_rows(Family, sep = ";")


phage_defense_pfams5 =  rsphagepaper5 %>% unique() %>% filter(str_starts(`Associated domains`, "pfam")) %>% 
    dplyr::select(`Associated domains`) %>% distinct() %>% arrange(`Associated domains`) %>%
    distinct() %>% 
    mutate(`Associated domains` = str_replace_all(`Associated domains`, "pfam", "PF")) %>% 
    separate_rows(`Associated domains`, sep = ";")


phage_defense_pfams_df = rbind(phage_defense_pfams1, phage_defense_pfams2, phage_defense_pfams3, phage_defense_pfams4) 
# Create a data frame from the list
phage_defense_pfams_df = phage_defense_pfams_df %>% distinct() %>% arrange(Family) %>% 
  filter(str_detect(Family, "^PF")) %>% distinct()



########################################
# How many pfams are assocaited with phage defense systems
########################################
PF00365_fishers_LanM_rodeo %>% 
  filter(pfam %in% phage_defense_pfams_df$Family) %>% 
  filter(!grepl("AAA|Helix-turn-helix domain|Histidine kinase|antibiotic", desc)) %>%
  View()

# do the same for counts
counts_3 %>%
  filter(pfam %in% phage_defense_pfams_df$Family) %>% 
  filter(Z_score_control > 1.64 & significant == TRUE) %>% View
  filter(!grepl("AAA|Helix-turn-helix domain|Histidine kinase|antibiotic", desc)) 


########################################
# Align competence proteins
########################################
library(Biostrings)
library(msa)

input_file = "/data/san/data0/users/david/intelligence/lanthipeptide_rodeo_2020paper/PF05952.faa"
output_dir = "/data/san/data0/users/david/intelligence/figures"
proteins = readAAStringSet(input_file)

aligned_proteins = msa(proteins, type = "protein", method = "ClustalW", 
  verbose= TRUE, order= "aligned" )

output_base_name <- tools::file_path_sans_ext(basename(input_file))
output_tex = file.path(output_dir, paste0(output_base_name, "_alignment_output.tex"))
msaPrettyPrint(aligned_proteins, 
  file=output_tex, output="tex", verbose=FALSE, askForOverwrite=FALSE,
  showNames="left", showNumbering="none", showLogo="top",
  paperWidth = 12, paperHeight = 4,
               showConsensus="bottom", logoColors="hydropathy",
               shadingMode = "identical",shadingModeArg = 90)

tinytex::pdflatex(output_tex)



########################################
# Heatmap for defense systems
########################################
library(dplyr)
library(tidyr)
library(pheatmap)
library(data.table)
library(ggplot2)

########################################
# Load in rodeo table for class II lanthipeptide
########################################
LanM_40k_window = fread("data/processed/LanM_paper_40k_window.tsv")
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
LanM_phyla <- LanM_40k_window %>% 
  dplyr::select(phylum)
LanM_contigs <- LanM_40k_window %>%
  dplyr::select(Nucleotide_acc) %>%
  distinct()
fwrite(LanM_contigs, "data/processed/LanM_contigs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Load genus information
genus_df <- LanM_40k_window %>%
  dplyr::select(Nucleotide_acc, genus, family, order) %>% 
  distinct()
fwrite(genus_df, "tables/genus_df.tsv", sep = "\t")
# Load defense system data
defense_files <- list.files(path = "/data/san/data2/users/david/co_localisation/data/defense", 
                            pattern = "*finder_systems.tsv", full.names = TRUE, recursive = TRUE)

defense_df_files <- do.call(rbind, lapply(defense_files, function(file) {
  df <- fread(file, sep = "\t", header = TRUE)
  df$Nucleotide_acc <- basename(dirname(file))  # Extract accession
  return(df)
}))

# Merge with genus data
defense_df <- left_join(genus_df, defense_df_files, by = "Nucleotide_acc") %>% 
  dplyr::select(Nucleotide_acc, genus, type, sys_id, family, order)

# Identify accessions without a system & assign "no_system"
defense_df <- defense_df %>%
  mutate(type = ifelse(is.na(type), "no_system", type))

# Count unique sys_id per Nucleotide_acc for each (family, type) combination
defense_counts <- defense_df %>%
  group_by(family, type, Nucleotide_acc) %>% 
  summarise(unique_sys_count = n_distinct(sys_id), .groups = "drop")

# Sum unique sys_id counts for each (family, type) while filtering out missing family values
contingency_table <- defense_counts %>%
  filter(!is.na(family)) %>%
  group_by(family, type) %>%
  summarise(total_sys_count = sum(unique_sys_count), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = total_sys_count, values_fill = list(total_sys_count = 0)) %>%
  textshape::column_to_rownames("family")

# Convert to proportions within each family (sum to 1 per row)
prop_table <- prop.table(as.matrix(contingency_table), margin = 1) * 100

# Define top genera and types
top_20_genera <- c("Streptomycetaceae", "Streptococcaceae", "Paenibacillaceae", "Bacillaceae", 
                   "Oscillatoriaceae", "Pseudomonadaceae", "Lachnospiraceae",
                   "Lactobacillaceae", "Bifidobacteriaceae", "Staphylococcaceae", 
                   "Weeksellaceae",
                   "Oscillospiraceae", "Enterococcaceae", "Actinomycetaceae",
                   "Eubacteriaceae", "Desulfallaceae", "Aerococcaceae", "Eggerthellaceae",
                   "Thermomonosporaceae", "Nocardioidaceae")

top_20_types <- c("RM", "Cas", "Wadjet", "Avs", "CBASS", "Ceres", "Lanthiphage",
                  "AbiD", "AbiE", "AbiH", "Abi2", "AbiAlpha", "AbiJ", "AbiN", "AbiU",
                  "Gabija", "SanaTA", "Detocs", "MazEF", "pAgo", "Prometheus", "Zorya",
                  "Shedu", "Septu", "Lamassu-Fam", "SoFIC", "BREX", "DarTG", 
                  "SspBCDE", "VP1851", "Tiamat", "Thoeris", "RosmerTA", "Pycsar",
                  "PsyrTA", "PrrC", "Retron", "Anti_CRISPR", "Anti_RM", "Anti_Pycsar", "no_system")

# Subset of interest
filtered_genera <- intersect(rownames(prop_table), top_20_genera)
filtered_types <- intersect(colnames(prop_table), top_20_types)
filtered_prop_table <- prop_table[top_20_genera, top_20_types, drop = TRUE]
# drop the column name no-system
filtered_prop_table <- filtered_prop_table[, -which(colnames(filtered_prop_table) == "no_system")]

# drop no-system column here
fig5_a <- pheatmap(filtered_prop_table, 
  cluster_rows = TRUE, 
  treeheight_row = 0,
  cluster_cols = TRUE,
  treeheight_col = 0,
  display_numbers = FALSE, 
  color = c("white", colorRampPalette(c("grey88", "red"))(99)),
  angle_col = 45,
  fontsize = 7)

# Save the heatmap
ggsave("figures/defense-genus_heatmap.png", plot = fig5_a$gtable, width = 20, height = 8, units = "cm")



########################################
# subset the data for system that are in LanM_40k_window 
########################################
top_20_genera_40k <- c(
  "Streptomycetaceae", 
  "Actinomycetaceae",
  "Streptococcaceae",
  "Enterococcaceae",
  "Lactobacillaceae",
  "Paenibacillaceae", 
  "Bacillaceae", 
  "Peptoniphilaceae",
  "Lachnospiraceae",
  "Clostridiaceae",
  "Corynebacteriaceae",
  "Staphylococcaceae",
  "Weeksellaceae",
  "Flavobacteriaceae",
  "Micrococcaceae"
  )
top_20_types_40k <- c("RM","Wadjet","CBASS",
  "AbiD", "AbiE", "Abi2", "AbiQ",
  "Gabija", "ShosTA",
  "Hna",
  "VP1853",
  "Nhi","PD-T7-3",
  "Lamassu-Fam",
  "RosmerTA",
  "Anti_CRISPR","Anti_RM")

defense_df_subset_close <- do.call(rbind, lapply(defense_files, function(file) {
  df <- fread(file, sep = "\t", header = TRUE)
  df$Nucleotide_acc <- basename(dirname(file))
  return(df)
})) %>% 
  mutate(sys_begin_number = as.numeric(gsub(".*_(\\d+)$", "\\1", sys_beg))) %>%
  left_join(genus_df) %>% 
  dplyr::select(Nucleotide_acc, genus, type, subtype, activity, sys_id, family, order, sys_beg)

defense_df_subset_close$Protein_acc <- sub(".*_prot_(.*)_\\d+$", "\\1", defense_df_subset_close$sys_beg)

fwrite(defense_df_subset_close,
  "tables/defense_df_subset_close.tsv", sep = "\t")

defense_df_subset_40k = defense_df_subset_close %>%
  filter(Protein_acc %in% LanM_40k_window$Protein_acc)

defense_df_subset_40k = left_join(genus_df, defense_df_subset_40k, by = c("Nucleotide_acc","genus","family","order")) %>% 
  dplyr::select(Nucleotide_acc, genus, type, subtype, activity, sys_id, family, order, sys_beg)

defense_df_subset_40k <- defense_df_subset_40k %>%
  mutate(type = ifelse(is.na(type), "no_system", type))
defense_df_subset_40k <- defense_df_subset_40k[, c("Nucleotide_acc","family","type")] %>%
  distinct()


contingency_table_subset <- table(defense_df_subset_40k$family, defense_df_subset_40k$type)
prop_df_40kb <- prop.table(contingency_table_subset, margin = 1) * 100
prop_df_40kb[is.nan(prop_df_40kb)] <- 0 
df_40kb <- prop_df_40kb[top_20_genera_40k, top_20_types_40k]

fig5_b <- pheatmap(df_40kb, 
  cluster_rows = TRUE, 
  treeheight_row = 0,
  cluster_cols = TRUE,
  treeheight_col = 0,
  display_numbers = FALSE,
  color = c("white", colorRampPalette(c("grey88", "red"))(99)),
  angle_col = 45,
  fontsize = 6)

ggsave("figures/defense-genus_heatmap_40k.png",
  fig5_b, width = 11, height = 6, units = "cm", dpi=600)

########################################
# plot distance from LanM all
########################################
# Load defense system data
# read in contig lengths
contig_length = fread("tables/contig_lengths.tsv", sep=" ") %>%
  dplyr::rename(
    Nucleotide_acc = name,
    total_contig_length = length
  )

all_contigs_files <- list.files(path = "data/tsv", 
                           pattern = "*.tsv", full.names = TRUE, recursive = TRUE)

all_contigs_df <- do.call(rbind, lapply(all_contigs_files, function(file) {
  df <- fread(file, sep = "\t", header = TRUE)
  return(df)
})) %>%
dplyr::rename(
  Nucleotide_acc = nucleotide_acc,
  Protein_acc = protein_id) %>%
  group_by(Nucleotide_acc) %>%
  left_join(contig_length) 

defense_hmm_files <- list.files(path = "data/defense", 
                            pattern = "*hmmer.tsv", full.names = TRUE, recursive = TRUE)

defense_hmm_df <- do.call(rbind, lapply(defense_hmm_files, function(file) {
  df <- fread(file, sep = "\t", header = TRUE)
  return(df)
})) %>%
mutate(Protein_acc = sub(".*_prot_(.*)_\\d+$", "\\1", hit_id)) %>%
dplyr::rename(
  Nucleotide_acc = replicon)

LanM_proteins <- LanM_40k_window %>% filter(pfam == "PF05147") %>%
  pull(Protein_acc)  %>% unique()


# Classify proteins in all_contigs_df
def_lanm_df <- all_contigs_df %>%
  dplyr::select(Nucleotide_acc, Protein_acc, start, end, total_contig_length) %>%
  distinct() %>%
  mutate(type = case_when(
    Protein_acc %in% defense_hmm_df$Protein_acc ~ "defense",
    Protein_acc %in% LanM_proteins ~ "lanM",
  )) %>%
  drop_na() %>%
  mutate(midpoint = (start + end) / 2)

# Compute distances
get_distances <- function(df1, df2, label) {
  inner_join(df1, df2, by = "Nucleotide_acc", suffix = c("_query", "_lanM")) %>%
    mutate(distance = abs(midpoint_query - midpoint_lanM), type = label)
}

merged_df <- bind_rows(
  get_distances(filter(def_lanm_df, type == "defense"), filter(def_lanm_df, type == "lanM"), "defense"),
) %>%
  left_join(genus_df) %>%
  group_by(Nucleotide_acc, Protein_acc_query) %>%
  filter(distance == min(distance)) %>%
  ungroup()

fwrite(merged_df, "tables/defense_distance_all.tsv", sep = "\t")
# Print unique contig counts
cat("Unique defense contigs:", n_distinct(merged_df %>% filter(type == "defense") %>% pull(Nucleotide_acc)), "\n")
cat("Unique anti-defense contigs:", n_distinct(merged_df %>% filter(type == "anti-defense") %>% pull(Nucleotide_acc)), "\n")


# Plot distances
fig_defense_dist <- ggplot(merged_df, aes(x = distance, fill = type)) +
  geom_density(alpha = 0.7) +
  labs(x = "Distance (bp)", y = "Density") +
  theme_bw() +
  scale_fill_manual(values = c("defense" = "red", "anti-defense" = "goldenrod")) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  guides(fill = "none")

ggsave("figures/defense_distance_all.png",
  fig_defense_dist, width = 10, height = 3, units = "cm", dpi=600)


# plot distribution of total_contig_length
fig_total_contig_length <- ggplot(all_contigs_df, aes(x = total_contig_length)) +
  geom_density(alpha = 0.7, adjust = 5) +
  labs(x = "Total contig length (bp)", y = "Density") +
  theme_bw()

ggsave("figures/total_contig_length.png",
  fig_total_contig_length, width = 10, height = 3, units = "cm", dpi=600)

# Plot distances above median 
median(contig_length$total_contig_length)
mean(contig_length$total_contig_length)

library(ggpubr)
# plot the proteins as dots and lanM as a line
fig_proteins_lanM <- ggplot(merged_df, aes(x = midpoint_lanM, y = midpoint_query, color = type)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black") +
  labs(x = "LanM position (bp)", y = "Defense position (bp)") +
  theme_bw() +
  scale_color_manual(values = c("defense" = "red", "anti-defense" = "goldenrod")) +
  guides(color = "none") +
  stat_cor(method = "spearman", size = 3, color = "black") 

ggsave("figures/proteins_lanM_spearman.png",
  fig_proteins_lanM, width = 10, height = 8, units = "cm", dpi=600)



median_df_len <- merged_df %>% 
  filter(total_contig_length_query > median(contig_length$total_contig_length))

fig_defense_dist_median <- ggplot(median_df_len, aes(x = distance, fill = type)) +
  geom_density(alpha = 0.7) +
  labs(x = "Distance (bp)", y = "Density") +
  theme_bw() +
  scale_fill_manual(values = c("defense" = "red", "anti-defense" = "goldenrod")) +
  geom_vline(xintercept = 20000, linetype = "dotted", color = "black") +
  guides(fill = "none") 

ggsave("figures/defense_distance_all_above_median.png",
  fig_defense_dist_median, width = 10, height = 5, units = "cm", dpi=600)

# plot distances above IQ3
iq3_df_len <- merged_df %>% 
  filter(total_contig_length_query > quantile(contig_length$total_contig_length, 0.75))


fig_defense_dist_iq <- ggplot(iq3_df_len, aes(x = distance, fill = type)) +
  geom_density(alpha = 0.7, adjust=1/3) +
  labs(x = "Distance (bp)", y = "Density") +
  theme_bw() +
  scale_fill_manual(values = c("defense" = "red", "anti-defense" = "goldenrod")) +
  geom_vline(xintercept = 200000, linetype = "dotted", color = "black") +
  guides(fill = "none")

ggsave("figures/defense_distance_all_above_iq3.png",
  fig_defense_dist_iq, width = 10, height = 6, units = "cm", dpi=600)

# plot density of total_contig_length and density of distance on same plot
fig_total_contig_length_dist <- ggplot(contig_length, aes(x = total_contig_length)) +
  geom_density(alpha = 0.7, adjust = 5) +
  geom_density(data = merged_df, aes(x = distance), alpha = 0.7, color = "red", adjust = 5) +
  labs(x = "Length (bp)", y = "Density") +
  theme_bw()

ggsave("figures/total_contig_length_dist.png",
  fig_total_contig_length_dist, width = 10, height = 8, units = "cm", dpi=600)


# plot distribitopn of start_query and midpoint_lanM on same plot
iqr <- merged_df %>%
  filter(total_contig_length_lanM > median(merged_df$total_contig_length_lanM)) %>%
  left_join(genus_df) %>%
  group_by(family) %>%
  filter(n() > 30)
  

fig_start_midpoint <- ggplot(iqr ,aes(x = start_query)) +
  geom_density(alpha = 0.7, adjust=1) +
  geom_density(data = iqr, 
    aes(x = midpoint_lanM), alpha = 0.7, color = "red", adjust=2) +
  labs(x = "Position (bp)", y = "Density") +
  theme_bw() + 
  facet_wrap(~order, scales = "free_y")


# Define a window around midpoint_lanM (e.g., ±5000 bp)
range <- 200000

# Count the number of Protein_acc_query near midpoint_lanM
counts <- merged_df %>%
  mutate(is_near = abs(midpoint_query - midpoint_lanM) <= range) %>%
  group_by(Nucleotide_acc) %>%
  summarise(protein_count = sum(is_near, na.rm = TRUE)) %>%
  arrange(desc(protein_count))

# Get the top Nucleotide_acc with the most Protein_acc_query near midpoint_lanM
top_nucleotide_acc <- counts %>%
  slice_max(order_by = protein_count, n = 1)

print(top_nucleotide_acc)


S_pyogenes <- merged_df %>%
  filter(Nucleotide_acc == "NZ_CAAJDZ010000002.1")
  
S_pyogenes_plot = ggplot(S_pyogenes) +
  geom_density(aes(x = start_query), color = "black", fill = "gray", alpha = 0.5, adjust = 0.5) + 
  geom_rug(aes(x = start_query), color = "black", alpha = 0.5) +
  geom_rug(aes(x = 314572), color = "red", size = 1.5) +
  geom_vline(xintercept = 314572, linetype = "dashed", color = "red") +
  labs(x = "Position (kbp)", y = "Density", 
    title = expression(italic("S. pyogenes") ~ "NZ_CAAJDZ010000002.1")) +
  theme_bw() +
  scale_x_continuous(labels = function(x) x / 1000)

ggsave("figures/S_pyogenes.png",
  S_pyogenes_plot, width = 10, height = 6, units = "cm", dpi=600)


# Filter the data for unique genera, only taking Streptococcus and Bacillus
genera_data <- merged_df %>%
  filter(!is.na(genus)) %>%
  filter(genus %in% c("Streptococcus", "Bacillus","Enterococcus")) %>%
  distinct(Nucleotide_acc, genus, start_query, midpoint_lanM)

# Create a plot with densities and start positions for each genus
genera_plot <- ggplot(genera_data, aes(x = start_query, color = genus, fill = genus)) +
  geom_density(alpha = 0.3, adjust = 0.5) +
  geom_rug(aes(x = start_query),color="grey", size = 1.2) +
  geom_rug(aes(x = midpoint_lanM), size = 1.2) +
  geom_vline(aes(xintercept = midpoint_lanM, color = genus), linetype = "dashed") +
  labs(
    x = "Position (kbp)",
    y = "Density",
    title = "Density Plot of Start Positions by Genus with LanM Positions"
  ) +
  theme_bw() +
  scale_x_continuous(labels = function(x) x / 1000) +
  theme(legend.position = "none")

# Save the plot
ggsave("figures/genera_density_plot.png",
  genera_plot, width = 12, height = 8, units = "cm", dpi = 600)


########################################
# plot distance from LanM for 40kb from systems
########################################
defense_df_proteins <- do.call(rbind, lapply(defense_files, function(file) {
  df <- fread(file, sep = "\t", header = TRUE)
  df$Nucleotide_acc <- basename(dirname(file))
  return(df)
})) %>% 
  mutate(sys_begin_number = as.numeric(gsub(".*_(\\d+)$", "\\1", sys_beg))) 
defense_df_proteins$Protein_acc <- sub(".*_prot_(.*)_\\d+$", "\\1", defense_df_proteins$sys_beg)

defense_whole_contig = defense_df_proteins %>%
  dplyr::select(Nucleotide_acc, activity, type, subtype)%>%
  left_join(., contig_length)


fwrite(defense_whole_contig, "tables/defense_whole_contig.tsv", sep = "\t") 

defense_only_df = subset(defense_df_proteins, activity == "Defense")
n_distinct(defense_only_df$Nucleotide_acc)

antidefense_only = subset(defense_df_proteins, activity == "Antidefense")
n_distinct(antidefense_only$Nucleotide_acc)


# Create a dataframe with relevant columns and add a 'type' column
def_lanm_dist <- LanM_40k_window %>%
  dplyr::select(Nucleotide_acc, Protein_acc, start, end) %>%
  distinct() %>%
  mutate(type = case_when(
    Protein_acc %in% defense_only_df$Protein_acc ~ "defense",
    Protein_acc %in% LanM_proteins ~ "lanM",
    Protein_acc %in% antidefense_only$Protein_acc ~ "anti-defense",
    TRUE ~ NA_character_
  ))

# Compute midpoints
def_lanm_dist <- def_lanm_dist %>%
  mutate(midpoint = (start + end) / 2)

# Split into lanM, defense, and antidefense genes
lanM_df <- def_lanm_dist %>% filter(type == "lanM")
defense_df <- def_lanm_dist %>% filter(type == "defense")
antidefense_df <- def_lanm_dist %>% filter(type == "anti-defense")

# Merge to find distances (inner join on Nucleotide_acc)
merged_defense_df <- inner_join(defense_df, lanM_df, by = "Nucleotide_acc", suffix = c("_def", "_lanM"))
merged_antidefense_df <- inner_join(antidefense_df, lanM_df, by = "Nucleotide_acc", suffix = c("_antidef", "_lanM"))

# Compute absolute distances
merged_defense_df <- merged_defense_df %>%
  mutate(distance = abs(midpoint_def - midpoint_lanM))

merged_antidefense_df <- merged_antidefense_df %>%
  mutate(distance = abs(midpoint_antidef - midpoint_lanM))

# Combine both dataframes for plotting
merged_df <- bind_rows(
  merged_defense_df %>% mutate(type = "defense"),
  merged_antidefense_df %>% mutate(type = "anti-defense")
)
n_distinct( merged_defense_df$Nucleotide_acc) # 149
n_distinct(merged_antidefense_df$Nucleotide_acc) # 149
n_distinct(merged_antidefense_df$Nucleotide_acc) # 149

# Plot distances
fig5_c = ggplot(merged_df, aes(x = distance, fill = type)) +
  geom_density(alpha = 0.7) +
  labs(
       x = "Distance (bp)",
       y = "Density") +
  theme_bw() +
  xlim(0, 30000) +
  scale_fill_manual(values = c("defense" = "red", "anti-defense" = "goldenrod")) +
  theme(legend.position = c(0.95, 0.95), legend.justification = c(1, 1)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  annotate("text", x = 0, y = Inf, label = "←distance→", vjust = 5, hjust = -0.1, color = "black", size=3)

ggsave("figures/defense_distance.png",
  width = 10, height = 6, units = "cm", dpi=600)


########################################
# Figure 5 (RM-system type elements)
########################################
defense_df_subset_40k %>%
  left_join(., LanM_40k_window %>% dplyr::select(Nucleotide_acc, genus,species)) %>%
  distinct() %>%
  group_by(Nucleotide_acc, genus,species) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))  %>% View



defense_focus = c(#"NZ_CP019655.1",
  "NZ_EQ973330.1","NZ_CAAIOF010000003.1",
  "NZ_PIEU01000022.1","NZ_LN890331.1",
  "NZ_PJTG01000002.1","NZ_LAWY01000037.1",
  "NZ_GG670360.1","NZ_POMF01000052.1",
  "NZ_KI669414.1","NZ_LUXJ01000054.1","NZ_AXDY01000014.1"
  )

defense_df_subset_close %>% 
  filter(Nucleotide_acc %in% defense_focus) %>%
  filter(Nucleotide_acc %in% "NZ_LN890331.1") 

LanM_40k_window %>%
  filter(Nucleotide_acc %in% "NZ_LN890331.1") %>%
  dplyr::select(Protein_acc,pfam,name,desc,start,end) %>%
  distinct() %>%
  arrange(start) %>% View

blast <- read_links("tables/ava_nr2.o6") %>%
  filter(seq_id %in% island_focus & seq_id2 %in% island_focus) %>%
  filter(seq_id != seq_id2)

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867")
dna_def_pfam = c("PF20473","PF01420","PF04851","PF02384","PF08463","PF08463","PF00270","PF01938","PF13649","PF07669","PF07669","PF20473",
  "PF02384","PF13337", "PF08665","PF05016")
dna_def_prot = c("WP_110028332.1","WP_110028334.1")
defense_proteins = defense_df_subset_close$Protein_acc


LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")

fig_5_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% defense_focus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = Nucleotide_acc, 
                strand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% dna_def_pfam ~ "red",
    Protein_acc %in% dna_def_prot ~ "red",
    Protein_acc %in% defense_proteins ~ "red",
  )) %>%
  fill(color, .direction = "downup") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) 

fig_5_plot <- gggenomes(fig_5_df) +
  geom_seq(color="black") +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3, color = "black", fontface = "italic") +
  geom_seq_label(size = 2, vjust = 2, color = "black", fontface = "bold") +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

fig_5_plot 

  
# plot figure 5
ggsave("figures/figure_5.png", 
  plot = fig_5_plot, 
  width = 24, height = 14, units = "cm", dpi=300)

svg("figures/figure_5.svg", width = 24 / 2.54, height = 14 / 2.54)  # Convert cm to inches
print(fig_5_plot)
dev.off()



########################################
# Figure Supplementary (toxin anti-toxin)
########################################
focus_PF05016 <- c( "NZ_QRCT01000014.1","NZ_FOUV01000039.1",
  "NZ_QRCT01000014.1", "NZ_JXYX01000009.1",
  "NZ_KB217478.1", "NZ_UGJA01000003.1"
  )

fig_S5  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% focus_PF05016) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                strand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% dna_def_pfam ~ "red",
    Protein_acc %in% dna_def_prot ~ "red",
  )) %>%
  fill(color, .direction = "downup") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) 

fig_S5_plot <- gggenomes(fig_S5) +
  geom_seq(color="black") +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3, color = "black", fontface = "italic") +
  geom_seq_label(size = 2, vjust = 2, color = "black", fontface = "bold") +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("figures/figure_S5.png", 
  plot = fig_5_plot,
  width = 14, height = 8, units = "cm", dpi=300)





########################################
# Figure 6 (competence)
########################################
competence_focus = c("NC_018081.1",
"NZ_CNVF02000013.1",
"NZ_NJFO02000003.1",
"NZ_JH792105.1",
"NZ_FUXA01000008.1",
"NZ_PIJH01000018.1",
"NZ_LMBZ01000008.1"
)

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867","PF07730")
competence_pfam = c("PF18146","PF00154","PF05952","PF06133","PF12072","PF02464","PF12072","PF05389")

fig_6_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% competence_focus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                strand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% competence_pfam ~ "red",
    Protein_acc %in% competence_focus ~ "red",
  )) %>%
  fill(color, .direction = "downup") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) %>%
  mutate(bin_id = case_when(Nucleotide_acc== "NZ_LMBZ01000008.1" ~ "Solibacillus cecembensis",
    TRUE ~ bin_id))

fig_6_df %>%
  filter(grepl("NZ_NJFO", Nucleotide_acc))  %>%
  arrange(desc(start)) %>% View()

fig_6_plot <- gggenomes(fig_6_df) %>%
  focus(pfam %in% lanthipeptide_pfam, .expand = c(10000,12000)) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3, color = "black", fontface = "italic") +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

fig_6_plot

# fig_6_plot_annotated <- fig_6_plot +
#   annotate("text", 
#     x = 18000, y = 7.4, 
#     label = "ComX", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 12000, y = 7.4,
#     label = "LanM-type BGC", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 21000, y = 6.36,
#     label = "CinA   RecA", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 25000, y = 6.36,
#     label = "ComK", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 8000, y = 5.36,
#     label = "CinA, RecA", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 28500, y = 4.36,
#     label = "MecA", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 19000, y = 3.36,
#     label = "ComX", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 9000, y = 2.36,
#     label = "ComX", color = "black", size = 2, fontface = "bold") +
#   annotate("text",
#     x = 4500, y = 1.36,
#     label = "RecA, CinA", color = "black", size = 2, fontface = "bold")

# ggsave("figures/figure_6.png", 
#   plot = fig_6_plot, 
#   width = 18, height = 10, units = "cm")

ggsave("figures/figure_6.png", 
  plot = fig_6_plot, 
  width = 18, height = 10, units = "cm", dpi=600)

svg("figures/figure_6.svg", width = 18 / 2.54, height = 10 / 2.54)  # Convert cm to inches
print(fig_6_plot)
dev.off()


########################################
# Figure Supplementary (sugar) (S4)
########################################
LanM_40k_window
sugarfocus = c(
  "NZ_CP028837.1",
  "NC_019042.1",
  "NZ_CP008926.1",
  "NZ_AODG01000003.1"
)

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867","PF07730")

sugar_pfam2 = LanM_40k_window %>%
  filter(grepl("PTS", desc))  %>% 
  dplyr::select(pfam, desc) %>% distinct
sugar_pfam = c("PF00834","PF02502","PF00294","PF00294", sugar_pfam2$pfam)
sugar_prot = c("WP_086950606.1",
"WP_086950605.1")

figure_S4_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% sugarfocus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                stand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% sugar_pfam ~ "red",
    Protein_acc %in% sugar_prot ~ "red",
  )) %>%
  fill(color, .direction = "downup") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) %>%
  mutate(bin_id = case_when(Nucleotide_acc== "NZ_LMBZ01000008.1" ~ "Solibacillus cecembensis",
    TRUE ~ bin_id))

figure_S4_plot <- gggenomes(figure_S4_df) %>%
  focus(pfam %in% lanthipeptide_pfam, .expand = c(10000,12000)) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3, color = "black", fontface = "italic") +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("figures/figure_S4.png", 
  plot = figure_S4_plot, 
  width = 16, height = 7, units = "cm", dpi=600)



########################################
# Table Supplementary (amrfinder plus)
########################################
# read in all .tsv files from "data/amrfinder_out"
amrfinder_files <- list.files(path = "data/amrfinder", 
  pattern = "*.tsv", full.names = TRUE)
# lapply and fread all files into a single df
amrfinder_df <- do.call(rbind, lapply(amrfinder_files, fread, sep = "\t", header = TRUE)) %>%
  janitor::clean_names() %>%
  dplyr::rename(Nucleotide_acc = "contig_id")

LanM_windows = LanM_40k_window %>% 
  dplyr::select(Nucleotide_acc, unique_window, start_window, end_window,
    superkingdom,phylum,class,order,family,genus,species) %>%
  distinct()

amr_df = left_join(LanM_windows, amrfinder_df, 
  by = "Nucleotide_acc",
  relationship = "many-to-many") %>%
  group_by(Nucleotide_acc, unique_window) %>%
  # filter start so it is between start_window and end_window
  filter(start >= start_window & start <= end_window)

# intersect amrfinder_df with LanM_40k_window based on Nulcoeitde_acc, start and stop
# the start from amrfinder must lie within the start_window and end_window of LanM_40k_window
# and this needs to be unique for each nucleotide_acc


########################################
# Acetyltransferase
########################################

acetyl_accessions = LanM_40k_window %>% 
  filter(grepl("GNAT", desc)) %>%
  dplyr::select(Protein_acc) %>%
  distinct()

fwrite(acetyl_accessions, "data/acetyltransferases/acetyltransferases.txt",
  col.names = FALSE, row.names = FALSE, quote = FALSE)

acetyl_pfam = LanM_40k_window %>% 
  filter(grepl("GNAT", desc)) %>%
  dplyr::select(pfam) %>%
  distinct()

########################################
# Acetyltransferase
########################################
## install the latest version from github



########################################
# Subsample
########################################
rare_df <- LanM_40k_window %>% 
  ungroup() %>%
  filter(grepl("PF", pfam)) %>%
  dplyr::select(pfam, desc, Nucleotide_acc) %>%
  distinct() 


library(dplyr)
library(ggplot2)
library(ggpubr)
# Function to calculate rarefaction curve for a given Pfam, counting occurrences
rarefaction_curve_for_pfam <- function(df, pfam_of_interest, max_sample_size) {
  rarefaction_results <- data.frame(sample_size = integer(), pfam_count = integer(), pfam = character())
  
  for (i in 1:max_sample_size) {
    # Randomly sample 'i' accessions without replacement
    sampled_accessions <- sample(unique(df$Nucleotide_acc), i, replace = FALSE)
    
    # Filter for the Pfam of interest and the sampled accessions
    sampled_pfams <- df %>% filter(Nucleotide_acc %in% sampled_accessions, pfam == pfam_of_interest)
    pfam_count <- nrow(sampled_pfams)  # Count the occurrences of the Pfam in the sample
    
    # Store the result
    rarefaction_results <- rbind(rarefaction_results, data.frame(sample_size = i, pfam_count = pfam_count, pfam = pfam_of_interest))
  }
  
  return(rarefaction_results)
}

set.seed(123)  # Set seed for reproducibility
random_pfam <- sample(unique(LanM_40k_window$pfam), 3)


pfams_of_interest <- c(
  "PF00365", # negative control
  random_pfam, # negative control
  "PF02687", "PF00005", "PF02687", # transport
  "PF03412", # peptidase
  # "PF10412", T4SS
  "PF05147", # LanM
  "PF05016", # TA ParDE
  "PF06414", # Zeta toxin
  "PF04221", # antitoxin Rel
  "PF00583", # acetyltransferase
  "PF15731", # antitoxin
  "PF01420", # RM
  "PF07669", # Eco57I restriction-modification methylase
  "PF04313", # Type I restriction enzyme R protein N terminus (HSDR_N)
  "PF00589", # phage integrase
  "PF00665", # recombinase
  "PF13102", # phage integrase
  "PF00154" # recA
  )

all_rarefaction_data <- data.frame()
for (pfam in pfams_of_interest) {
  rarefaction_data <- rarefaction_curve_for_pfam(rare_df, pfam, max_sample_size = 100)
  all_rarefaction_data <- rbind(all_rarefaction_data, rarefaction_data)
}


# use spearman as its non-parametric and can assess monotonic relationships
corr_plot = ggplot(all_rarefaction_data, aes(x = sample_size, y = pfam_count)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color="red") +
  labs(title = "",
       x = "Number of Accessions Sampled",
       y = "Pfam Count (Occurrences)") +
  theme_classic() +
  theme(text = element_text(size = 10)) +
  facet_wrap(~pfam, scales = "free_y") +
  stat_cor( method = "spearman", size = 3, color = "black")
  

ggsave(plot = corr_plot,
  filename = "figures/rarefaction_curves.png", 
  width = 22, 
  dpi = 600,
  height = 20, 
  units = "cm")


########################################
# Stratify by genus
########################################
top_families <- LanM_40k_window %>% 
  group_by(family) %>% 
  summarise(n_accessions = n_distinct(Nucleotide_acc)) %>% 
  arrange(desc(n_accessions)) 


plot_rarefaction_by_genus <- function(genus_name, pfams, max_sample_size = 20) {
  df_family <- LanM_40k_window %>% 
    ungroup() %>%
    filter(grepl("PF", pfam)) %>%
    dplyr::select(pfam, desc, Nucleotide_acc, family, genus) %>%
    distinct() %>%
    group_by(family) %>%
    filter(n_distinct(Nucleotide_acc) > 10) %>%
    ungroup() %>%
    filter(family == genus_name)
  
  data_family <- do.call(rbind, lapply(pfams, function(pfam) {
    rarefaction_curve_for_pfam(df_family, pfam, max_sample_size)
  }))
  
  ggplot(data_family, aes(x = sample_size, y = pfam_count)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
    labs(
      title = paste("Rarefaction curve for", genus_name),
      x = "Number of Accessions Sampled",
      y = "Pfam Count (Occurrences)"
    ) +
    theme_classic() +
    theme(text = element_text(size = 10)) +
    facet_wrap(~pfam, scales = "free_y") +
    stat_cor(method = "spearman", size = 3, color = "black")
}

unique(LanM_40k_window$family) %>% sort()
corr_plot_Enterococcus <- plot_rarefaction_by_genus("Enterococcaceae", pfams_of_interest, max_sample_size = 27)
corr_plot_Streptococcus <- plot_rarefaction_by_genus("Streptococcaceae", pfams_of_interest, max_sample_size = 135)
corr_plot_Bacillus <- plot_rarefaction_by_genus("Bacillaceae", pfams_of_interest, max_sample_size = 141)
corr_plot_Streptomyces <- plot_rarefaction_by_genus("Streptomycetaceae", pfams_of_interest, max_sample_size = 44)
corr_plot_Paenibacillus <- plot_rarefaction_by_genus("Paenibacillaceae", pfams_of_interest, max_sample_size = 27)


# Determine the top 10 genera based on the number of distinct accessions for pfams_of_interest
top_10_family <- LanM_40k_window %>%
  filter(pfam %in% pfams_of_interest) %>%
  group_by(family) %>%
  summarise(total_accessions = n_distinct(Nucleotide_acc)) %>%
  arrange(desc(total_accessions)) %>%
  filter(!is.na(family)) %>%
  head(10) %>%
  pull(family)

# Plot the count of each pfam per genus, faceted by pfam and including only the top 10 genera with free y-axis scales
pfam_count_per_genus <- LanM_40k_window %>%
  filter(pfam %in% pfams_of_interest, genus %in% top_10_family) %>%
  group_by(genus, pfam) %>%
  summarise(n_accessions = n_distinct(Nucleotide_acc)) %>%
  ggplot(aes(x = genus, y = n_accessions, fill = family)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  facet_wrap(~ pfam, scales = "free_y") +
  labs(
    title = "Number of Accessions per Genus with Pfam",
    x = "Genus",
    y = "Number of Accessions"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/pfam_count_per_family.png", 
  plot = pfam_count_per_genus, 
  width = 24, height = 16, units = "cm", dpi=600)
