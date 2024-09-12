########################################
## Packages - don't package shame people
########################################
.libPaths(c("/data/san/data0/users/david/rstudio/packages", .libPaths()))
newlib <- "/data/san/data0/users/david/rstudio/packages"



packages <- c("cowplot","data.table", "multidplyr", "readr", "formattable", "fs", "dplyr", "ggplot2", "purrr", "ggthemes", "BiocManager", "gplots", "gridExtra", "grid", "forcats", "tidyr", "dtplyr", "topGO", "SparseM", "Biostrings", "GenomicRanges", "seqinr", "stringr", "readxl", "thacklr", "gggenomes")

load_packages <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org", lib = newlib)
      library(package, character.only = TRUE)
    }
  }
}
load_packages(packages)


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
setwd("/data/san/data0/users/david/taxonomizer") #SET UP DATABASE ONCE #
sqlFile <- "/data/san/data0/users/david/taxonomizer/accessionTaxa.sql"
  
append_taxonomy <- function(df, sqlFile) {
  df$taxid <- taxonomizr::accessionToTaxa(df$subject.x, sqlFile)
  taxa <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  for(taxon in taxa) {
    df[[taxon]] <- getTaxonomy(df$taxid, sqlFile, desiredTaxa = taxon)
  }
  return(df)
}


setwd("/data/san/data0/users/david/intelligence/code")

########################################
# Load in rodeo table for class II lanthipeptide
########################################
LanM_40k_window = fread("data/processed/LanM_paper_40k_window.tsv")
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
LanM_phyla <- LanM_40k_window %>% 
  dplyr::select(phylum)


########################################
# Load in PF00365_df
########################################
PF00365_df <- fread("data/negative_dataset/PF00365_df.tsv", sep = "\t") 

n_distinct(PF00365_df$Nucleotide_acc)
n_distinct(LanM_40k_window$Nucleotide_acc)


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

LanM_table <- LanM_40k_window %>% 
  dplyr::select(Nucleotide_acc, pfam) %>% 
  distinct()
LanM_table$group <- "lanII bacteriocin"
LanM_table$id <- "lan"


ctrl_v_lanm = rbind(PF00365_table, LanM_table) 
ctrl_df = rbind(PF00365_table)

write_tsv(ctrl_v_lanm, 
  "data/processed/ctrl_v_LanM_vs_Pks.tsv")
write_tsv(ctrl_df, 
  "data/processed/LanM_ctrl_df_pfams_nucs.tsv")




########################################
# Plot what the control df looks like
########################################
library(ggsci)
set.seed(1)
fig_1a <- PF00365_df %>%
filter(Nucleotide_acc %in% PF00365_table$Nucleotide_acc) %>%
  group_by(Nucleotide_acc) %>%
  dplyr::summarise(unique = n_distinct(locus_tag)) %>%
  ggplot(aes(x = "All Genomes", y = unique)) +
  geom_jitter(width = 0.3, size = 0.3, alpha = 0.3) +
  geom_violin(alpha = 0.7, fill="goldenrod") +
  theme_bw() +
  labs(x = "", y = "Protein Count", title = "PF00365 control") +
   theme(
    plot.title = element_text(size = 12),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10)   # Adjust size as needed for y-axis labels
  ) 
ggsave(fig_1a, filename = "figures/PF00365_control_window_protein_count_LanM.png", 
    width = 6, height = 10, units = "cm")

########################################
# Plot what the LanM df looks like
########################################
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
fig_1b <- LanM_40k_window %>%
  group_by(unique_window) %>%
  summarise(unique = n_distinct(Protein_acc)) %>%
  ggplot(aes(x = "All Genomes", y = unique)) +
  geom_jitter(width = 0.3, size = 0.3, alpha = 0.5) +
  geom_violin(alpha = 0.7, fill="red") +
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
ggsave(fig_1b, file="figures/LanM_bacteriocin_window_protein_count.png",
  width = 6, height = 10, units = "cm")


fig_1c <- LanM_40k_window %>%
  dplyr::select(genus, Nucleotide_acc,start_window) %>% 
  filter(!is.na(genus)) %>%
  distinct() %>%
  group_by(genus) %>%
  dplyr::summarise(count = n()) %>% 
  filter(count > 10) %>%
  arrange(desc(count)) %>%
  ggplot(., aes(x = reorder(genus, count), y = count)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("")

ggsave(fig_1c, filename = "figures/LanM_paper_40k_window_genus_count.png",
    height = 4, width = 4)



# Figure 1
cowplot_figure_1_top = cowplot::plot_grid(fig_1a, fig_1b,
  labels = c("(a)", "(b)"), 
  ncol = 2,
  label_size = 10)
cowplot_figure_1_bottom = cowplot::plot_grid(fig_1c, 
  labels = c("(c)"), 
  ncol = 1,
  label_size = 10)
cowplot_figure_1 <- plot_grid(cowplot_figure_1_top, 
  cowplot_figure_1_bottom, ncol = 1, 
  rel_heights = c(0.5, 0.5))  

# send final figure to github repo
ggsave2("figures/figure_1.png", 
  cowplot_figure_1, width = 12, height = 20, units = "cm", dpi=300)
ggsave2("figures/figure_1.tiff", 
  cowplot_figure_1, width = 12, height = 20, units = "cm", dpi=300)


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

  # Select the relevant columns from df and add a group column
  x <- x %>%
    dplyr::select(Nucleotide_acc, pfam) %>%
    distinct()
  x$group <- control_name

  # make contingency table for pfams
  k <- rbind(LanM_table, x) %>% 
    distinct() 
  # Filter out rows with complete cases
  k <- k[complete.cases(k), ]  
  # Replace NA values
  k[is.na(k)] <- "Unknown"

  contingency_table <- table(k$group, k$pfam)
  total_df <- sum(contingency_table[1, ])
  total_control <- sum(contingency_table[2, ])
  
  cat("\nTotal for lanII bacteriocin:", total_df)
  cat("\nTotal for", control_name, ":", total_control, "\n")
  
  return(contingency_table)
}

# Apply the function to each data frame
cont_control_PF00365 <- create_contingency_table(PF00365_table, LanM_40k_window, "PF00365")


########################################
# Counts vs a Random distribution of counts i.e. is there more 
########################################

LanMin = LanM_40k_window %>%
	dplyr::select(Nucleotide_acc, pfam) %>%
	mutate(group = "lanthipeptide",
		id = "lanthipeptide II") %>%
	filter(grepl("^PF", pfam)) # only take Pfams

tablex = rbind(LanMin, PF00365_table) %>% dplyr::select(-id)
fwrite(tablex, "data/processed/lanthipeptide_pfam_table.tsv")

# check point
tablex = fread("data/processed/lanthipeptide_pfam_table.tsv")


plot_size_of_pfam_pools_df = tablex %>%
  group_by(group) %>%
  summarize(total = n())
# these are unequal in size

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
counts <- counts %>%
  mutate(
    m_prime_lanthipeptide = total_genes * f_lanthipeptide,
    s_prime_lanthipeptide = sqrt(total_genes * f_lanthipeptide * (1 - f_lanthipeptide)),
    m_prime_control = total_genes * f_control,
    s_prime_control = sqrt(total_genes * f_control * (1 - f_control))
  )

# Calculate Z-scores for both groups
significant_threshold <- 0.05


counts <- counts %>%
  mutate(
    Z_score_lanthipeptide = (genes_in_lanthipeptide - m_prime_lanthipeptide) / s_prime_lanthipeptide,
    Z_score_control = (genes_in_control - m_prime_control) / s_prime_control,
    p_value = pnorm(-abs(Z_score_lanthipeptide)),
    bonferroni_p_value = pmin(1, p_value * n()),
    p_value_control = pnorm(-abs(Z_score_control)),
    bonferroni_p_value_control = pmin(1, p_value_control * n()),
    significant = bonferroni_p_value < significant_threshold
  ) %>% 
  filter(!is.na(pfam))

# View results
print(counts)
fwrite(counts, "tables/lanthipeptide_pfam_counts.tsv")

# plot histogram 
histogram_z_scores <- ggplot(counts, 
    aes(
    x = Z_score_lanthipeptide, 
    fill = significant)) +
  geom_histogram(
    bins = 100, 
    alpha = 0.6, 
    color = "black", 
    size=0.1) +
  labs(
    x = "Z-score", 
    y = "Frequency",
    title = "Histogram of Z-scores for Pfam Families") +
  scale_fill_manual(values = c("grey", "red", "darkgreen")) +
  theme_classic(base_size = 8)

ggsave(filename = "figures/histogram_z_scores.png", 
  plot = histogram_z_scores, width = 10, height = 6, units = "cm", dpi=300)



# plot histogram 
ggplot(counts, 
    aes(
    x = (genes_in_lanthipeptide - genes_in_control)/total_genes, 
    fill = significant)) +
  geom_histogram(
    bins = 100, 
    alpha = 0.6, 
    color = "black", 
    size=0.1) +
  labs(
    x = "x", 
    y = "y",
    title = "Histogram of Z-scores for Pfam Families") +
  scale_fill_manual(values = c("grey", "red", "darkgreen")) +
  theme_classic(base_size = 8)

ratio_df = counts %>% 
  group_by(pfam) %>%
  mutate(k = (genes_in_lanthipeptide/total_genes)) %>%
  arrange(k) 



# extreme values of Z-scores highlight that the distribution is not normal
# but that is not a problem for the test
# as we are just stating that counts are orders of magnitude above the mean

########################################
# Plot expected vs observed counts
########################################
expected_vs_observed_plot = ggplot(counts) +
   geom_point(aes(
    x = m_prime_control, 
    y = m_prime_control,
    fill = "Expected",
    ), 
    alpha = 0.3,
    shape=21) +
  geom_point(aes(
    x = m_prime_control, 
    y = genes_in_control,
    fill = "Control"), 
    shape=21,
    alpha = 0.3) +
  geom_point(aes(
    x = m_prime_lanthipeptide, 
    y = genes_in_lanthipeptide, 
    fill = "Lanthipeptide"), 
    shape=21,) +
  scale_fill_manual("", values = c("Expected" = "darkgrey",
                                "Control" = "goldenrod", 
                                "Lanthipeptide" = "red")) +
  geom_abline(slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "gray") +
    coord_cartesian(clip = "off") +
  labs(x = "Expected Number of Pfams", y = "Observed Number of Pfams",
       title = "Expected vs. Observed Pfams in 'lanthipeptide'") +
  theme_classic(base_size = 8) +
  geom_text(data = subset(counts, (name == "ABC transporter" & significant == TRUE) | 
    name == "Lanthionine synthetase C-like protein"),
             aes(x = m_prime_lanthipeptide, 
             y = genes_in_lanthipeptide, 
             label = name),
             color = "black",
             fontface = "bold",
             size=2,
             vjust = -0.8,
             position = position_nudge(y = 1.5))  + # Labels for "ABC transporter"
    geom_text(data = subset(counts, (name == "ABC transporter" & significant == TRUE) | 
    name == "Lanthionine synthetase C-like protein"), 
             aes(x = m_prime_control, y = genes_in_control, label = name),
             color = "black",
             fontface = "bold",
             size=2,
             vjust = -0.8,
             position = position_nudge(y = 1.5)) 

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/expected_vs_observed_plot.png", 
  plot = expected_vs_observed_plot, width = 10, height = 7, units = "cm")


subset_count = counts %>%
  filter(m_prime_lanthipeptide < 500 &  m_prime_control < 500) %>%
  filter(genes_in_lanthipeptide > m_prime_lanthipeptide) %>%
  filter(bonferroni_p_value < 0.05)

expected_vs_observed_plot_subset = ggplot(subset_count) +
  geom_point(aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide),
    shape=21,
    color = "black",
    # alpha = (1- p_value),
    fill="red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Expected Number of Pfams", y = "Observed Number of Pfams", 
       title = "Expected vs. Observed Pfams in 'lanthipeptide'")  +
  theme_classic(base_size = 8) +
  geom_label(data = subset(subset_count, (name == "Competence-damaged protein") | 
    desc_pfam == "Methylase_S" | desc_pfam == "Methyltransf_25" |
    desc_pfam == "Mersacidin"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
            #  fontface = "italic",
             size=2,
            #  hjust=-.011,
            #  vjust = -0.8,
             position = position_nudge(y = 1.5)) +
  geom_label(data = subset(subset_count, desc_pfam == "AbiJ_NTD3"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
             size=2,
             position = position_nudge(y = 120)) +
  geom_label(data = subset(subset_count, desc_pfam == "MqsA_antitoxin"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
             size=2,
             position = position_nudge(y = 180)) +
  geom_label(data = subset(subset_count, desc_pfam == "Eco57I"), 
             aes(x = m_prime_lanthipeptide, y = genes_in_lanthipeptide, label = desc_pfam),
             color = "black", 
             size=2,
             position = position_nudge(y = 260))

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/expected_vs_observed_plot_subset.png", 
  plot = expected_vs_observed_plot_subset, width = 10, height = 7, units = "cm")


# Figure 2
cowplot_figure_2_top = cowplot::plot_grid(fig_1a, fig_1b,
  labels = c("(a)", "(b)"), 
  ncol = 2,
  label_size = 10)
cowplot_figure_2_bottom = cowplot::plot_grid(fig_1c, 
  labels = c("(c)"), 
  ncol = 1,
  label_size = 10)
cowplot_figure_2 <- plot_grid(cowplot_figure_2_top, 
  cowplot_figure_2_bottom, ncol = 1, 
  rel_heights = c(0.5, 0.5))


ggsave2("/data/san/data0/users/david/intelligence/code/figures/figure_2.png", 
  cowplot_figure_1, width = 12, height = 20, units = "cm", dpi=300)




########################################
# Q-Q plot
########################################
# Filter data for lanthipeptide group and count Pfam families
lanthipeptide_counts <- table(tablex %>% filter(group == "lanthipeptide") %>% .$pfam)
control_counts <- table(tablex %>% filter(group == "control") %>% .$pfam)
lanthipeptide_counts <- as.numeric(lanthipeptide_counts)
control_counts <- as.numeric(control_counts)


qqplot(lanthipeptide_counts, control_counts, main = "Q-Q Plot", xlab = "Lanthipeptide", ylab = "Control")
abline(0, 1, col = "red")
length(lanthipeptide_counts)
length(control_counts)

# Shapiro-Wilk test for normality
shapiro_test_lanthipeptide <- shapiro.test(lanthipeptide_counts)
shapiro_test_control <- shapiro.test(sample(control_counts, size = 4000) )

# Print the p-values
print(paste("Lanthipeptide p-value:", shapiro_test_lanthipeptide$p.value))
print(paste("Control p-value:", shapiro_test_control$p.value))

# If the lengths are not the same, pad the shorter one with zeros
max_length <- max(length(lanthipeptide_counts), length(control_counts))

lanthipeptide_counts <- c(lanthipeptide_counts, rep(0, max_length - length(lanthipeptide_counts)))
control_counts <- c(control_counts, rep(0, max_length - length(control_counts)))

# Perform the Q-Q plot
qqplot(lanthipeptide_counts, control_counts, main = "Q-Q Plot", xlab = "Lanthipeptide", ylab = "Control")
abline(0, 1, col = "red")


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


########################################
# stats on data
########################################
def = fread("tables/phage_defense.txt", header=FALSE, 
  col.names = c("pfam"))
pfams_proteins = LanM_40k_window %>% 
  dplyr::select(-Query) %>%
  distinct()

# what are the highest pfams
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
melted_data <- melt(tbl)
colnames(melted_data) <- c("group", "pfam", "count")
melted_expected <- melt(expected_mc)
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
ratio_plot_observed_vs_expected <- ggplot(comparison, aes(x=reorder(pfam, -ratio), y=ratio, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("goldenrod","red")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x.top = element_blank(),
  	axis.text.x.bottom = element_blank(),
    axis.title.x = element_text()
    ) +  
  theme_classic() + 
  xlab("Pfam") +
  labs(title="Ratio of Observed to Expected Counts")

ggsave(filename = "/data/san/data0/users/david/intelligence/figures/ratio_plot_observed_vs_expected.png", 
	plot = ratio_plot_observed_vs_expected, 
  width = 15, height = 10, units = "cm")

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
  rownames_to_column(var = "pfam")
 
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

# fwrite(PF00365_fishers_LanM_rodeo, "/data/san/data0/users/david/intelligence/tables/PF00365_fishers_LanM_rodeo.tsv", sep = "\t")
fwrite(PF00365_fishers_LanM_rodeo, "/data/san/data0/users/david/intelligence/code/PF00365_fishers_LanM_rodeo.tsv", sep = "\t")



########################################
# Volcano Plot
########################################
library(ggrepel)

vplot1 = PF00365_fishers_LanM_rodeo %>%
  filter(adjusted_p_value < 0.05) %>%
  filter(proportion_diff > 0) %>%
  left_join(., (dplyr::select(LanM_40k_window, desc, pfam) %>% 
  distinct())) %>%
  arrange(desc(adjusted_p_value)) %>% 
  filter(control == "PF00365" | control == "PF01225") %>%
  ggplot(aes(x = log10(proportion_diff), y = -log10(adjusted_p_value))) +
  geom_jitter(
    alpha = 0.9,
    size = 1, 
    aes(color=control)) +
  geom_label_repel(data = . %>% filter(desc %in% c(
		"Lanthionine synthetase C-like protein",
		"Type-A lantibiotic",
     "Competence-damaged protein",
     "Bacillus competence pheromone ComX",
		"ABC transporter")), 
            aes(label = desc), vjust = 1, hjust = 1, size = 2) +
  xlab("log10(Proportion difference)") +
  ylab("-log10(Adjusted P-Value)") +
  ggtitle("Volcano Plot: LanM-associated Pfams") +
  scale_color_manual(values = c("black", "#a22d2d")) +
  theme_classic(
    base_size = 8, 
    base_family = "arial") +
  theme(legend.position = "none",
    plot.title = element_text(size = 10), 
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  coord_cartesian(clip="off")

ggsave("/data/san/data0/users/david/intelligence/figures/LanM_controlPk_volcano_fishers_prop_diff_2controls_pfam.png", 
	vplot1, width = 10, height = 7, units = "cm", dpi=300)

vplot2 = PF00365_fishers_LanM_rodeo %>% 
  filter(proportion_diff < 0.0015) %>%
  filter(adjusted_p_value < 0.05) %>% 
  filter(proportion_diff > 0) %>%
  ggplot(aes(x = log10(proportion_diff), y = -log10(adjusted_p_value))) +
  geom_jitter(alpha = 0.9, 
    shape=21,
    size = 1, 
    aes(color=control,
      fill= "#d4d4d4")) +
    geom_label_repel(data = . %>% filter(desc %in% c(
      "Zeta toxin",
      "Competence-damaged protein",
      "Type I restriction modification DNA specificity domain",
      "Type I restriction and modification enzyme - subunit R C terminal",
      "Competence-damaged protein")), 
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

ggsave("/data/san/data0/users/david/intelligence/figures/LanM_controlPk_volcano_fishers_prop_diff_pfam_zoomed.png", 
	vplot2, width = 10, height = 7, units = "cm", dpi=300)


# cow plot both volcano plots
cowplot_x1 = cowplot::plot_grid(vplot1, vplot2, 
  labels = c("(a)", "(b)"), 
  ncol = 2,
  label_size = 10)

ggsave2("/data/san/data0/users/david/intelligence/figures/LanM_controlPk_volcano_fishers_prop_diff_pfam_cowplot.png", 
  cowplot_x1, width = 16, height = 7, units = "cm", dpi=300)


##################################
# Random Forrest for PFAMs vs PF00365
#################################
library(randomForest)
library(ranger)
library(caret)
library(ROCR)
library(pROC)


pfam_descs = LanM_40k_window %>% 
	dplyr::select(pfam, desc) %>% 
	group_by(pfam, desc) %>%
	summarise(n = n()) %>% arrange(desc(n)) %>%
	as.data.frame()

rm_these = pfam_descs %>%
	dplyr::slice(1:50) %>%
	dplyr::select(pfam) %>%
	as.vector()


top_50_pfams_per_group <- tablex %>%
  group_by(group, pfam) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(group, desc(count)) %>%
  group_by(group) %>%
  slice_head(n = 50) %>%
  ungroup()

# Save the top 100 Pfams to a dataframe called "remove_these"
remove_these <- top_50_pfams_per_group %>%
  arrange(desc(count)) %>%
  slice_head(n = 100)

# Print the "remove_these" dataframe
print(remove_these)

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Run RF Model for 5 seeds 
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# make circ_v_control_wo_query a hot matrix 
control <- circ_v_control_wo_query %>%



perform_iterations <- function(circ_v_control_wo_query, num_iterations, seed_values) {
  all_roc_data <- list()
  feature_importance_list <- list()
  for (seed in seed_values) {
    set.seed(seed)
    
    training_rows <- caret::createDataPartition(circ_v_control_wo_query$group, p = 0.7, list = FALSE)
    training_data <- circ_v_control_wo_query[training_rows, ]
    testing_data <- circ_v_control_wo_query[-training_rows, ]
    
    rf_model <- ranger::ranger(group ~ pfam, data = training_data,
                               importance = "impurity", num.threads = 32,
                               probability = TRUE)
    
    predictions <- predict(rf_model, data = testing_data, num.threads = 32, type = "response")
    prob_control <- predictions$predictions[, "control"]
    
    roc_obj <- pROC::roc(testing_data$group, prob_control)
    roc_df <- data.frame(
      FPR = roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Seed = seed
    )
    auc_val <- pROC::auc(roc_obj)
    
    all_roc_data[[as.character(seed)]] <- list(roc_df = roc_df, auc_val = auc_val)
    feature_importance_list[[as.character(seed)]] <- rf_model$variable.importance
  }
  
  return(all_roc_data)
}

circ_v_control_df <- tablex
circ_v_control_df <- na.omit(tablex)
circ_v_control_df <- circ_v_control_df %>% filter(!pfam %in% remove_these$pfam)
circ_v_control_df$pfam <- as.factor(circ_v_control_df$pfam)
circ_v_control_df$group <- as.factor(circ_v_control_df$group)

values_to_filter <- c("hypothetical", "PF00365", "PF00005",
                      "PF01225", "PF08245", "PF02875", "PF01225",
                      "PF13241", "PF04101", "PF03033", "PF00905",
                      "PF03717", "PF06271", "PF08245", "PF02875",
                      "PF03033", "PF00905", "PF03717", "PF01795",
                      "PF01098", "PF00091", "PF12327", "PF00953",
                      "PF08478", "PF07478", "PF13535", "PF14450",
                      "PF06723", "PF01820", "PF09221", rm_these)

circ_v_control_wo_query <- circ_v_control_df %>%
  filter(!(pfam %in% values_to_filter))

seed_values <- c(123, 456, 789, 987, 654)  # Add your desired seed values
roc_data <- perform_iterations(circ_v_control_wo_query, num_iterations = 5, seed_values = seed_values)

# Create a combined ROC curve plot
all_roc_df <- do.call(rbind, lapply(roc_data, function(x) x$roc_df))
average_auc <- data.frame(Seed = names(roc_data), AUC = sapply(roc_data, function(x) x$auc_val))
average_auc_subset <- average_auc[average_auc$Seed == 123, ]
plot_ROC = ggplot(all_roc_df, aes(x = FPR, y = TPR, color = as.factor(Seed))) +
  geom_line() +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Seed") +
  ggtitle("ROC Curves - Pks") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),  # Adjust size as needed
    axis.title.x = element_text(size = 10),  # Adjust size as needed
    axis.title.y = element_text(size = 10),  # Adjust size as needed
    axis.text.x = element_text(size = 10),  # Adjust size as needed for x-axis labels
    axis.text.y = element_text(size = 10),   # Adjust size as needed for y-axis labels
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  scale_x_reverse() +
  geom_text_repel(data = average_auc_subset, aes(x = 0.8, y = 0.1 * as.numeric(as.factor(Seed)), 
                                              label = paste("Avg. AUC = ", round(AUC, 3))), 
                size = 3, hjust = -0.4,
                vjust=-20,
                segment.color = NA, 
                color = "black") +
  scale_color_brewer(palette = "Reds")

ggsave("/data/san/data0/users/david/intelligence/figures/LanM_bacteriocin_RF_predict_PKS_newDF.png", 
       plot_ROC, width = 10, height = 10, units = "cm")


# Extract and combine feature importance from all iterations
feature_importance <- do.call(rbind, lapply(results$feature_importance, function(x) as.data.frame(t(x))))
feature_importance <- feature_importance %>%
  rownames_to_column(var = "pfam") %>%
  gather(key = "Seed", value = "Importance", -pfam)

# Calculate average importance for each pfam
average_importance <- feature_importance %>%
  group_by(pfam) %>%
  summarise(Avg_Importance = mean(Importance, na.rm = TRUE)) %>%
  arrange(desc(Avg_Importance))

# Print or visualize the top contributing Pfams
top_contributing_pfams <- average_importance %>%
  slice_head(n = 10)  # Adjust the number as needed to see more or fewer Pfams

print(top_contributing_pfams)




# ~~~~~~~~~~~~~~~~~~~~~~~~
# VOLCANO PLOT
# ~~~~~~~~~~~~~~~~~~~~~~~~
cont_control_PF00365


# assuming `contingency_table` is your table
diff_pfam_PF00365 <- cont_control_PF00365[1,] - cont_control_PF00365[2,]
abs_diff_PF00365 <- abs(diff_pfam_PF00365)

# sorting to find the most significant differences
top_diff_PF00365 <- sort(abs_diff_PF00365, decreasing = TRUE)[1:100]


# ~~~~~~~~~~~~~~~~~~~~~~~~
# PFAM TO GO
# ~~~~~~~~~~~~~~~~~~~~~~~~
ctrl_v_lanm
library(ragp)
# Specify the Bioconductor version when installing
install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/topGO_2.52.0.tar.gz", repos = NULL, type = "source")
library(topGO)

ctrl_v_lanm 
circular_DF_GO <- pfam2go(data_pfam = ctrl_v_lanm, pfam = "pfam")
contingency_table_GO <- table(circular_DF_GO$id, circular_DF_GO$GO_acc)


# Number of top GO terms to select
top_n <- 30

# Find top N GO terms for each group
top_GO_terms <- circular_DF_GO %>%
  filter(!is.na(GO_name)) %>%
  filter(GO_name != "hypothetical") %>%
  filter(GO_name != "") %>%
  dplyr::select(GO_name, group, Nucleotide_acc, id) %>%
  group_by(id, GO_name) %>%
  summarise(total = n(), .groups = "drop") %>%
  group_by(id) %>%
  arrange(desc(total)) %>%
  top_n(top_n, total) 


# ~~~~~~~~~~~~~~~~~~~~~~~~
# GO (1) Biological Process
# ~~~~~~~~~~~~~~~~~~~~~~~~

library(GO.db)

go_ids <- as.character(circular_DF_GO$GO_acc)
# Get the GO categories (ontologies) corresponding to these IDs
go_categories <- select(GO.db, keys = go_ids, columns = "ONTOLOGY", keytype = "GOID")
colnames(go_categories) <- c("GO_acc", "category")
go_categories = distinct(go_categories)

# filter our proteins with pfams
values_to_filter

circular_DF_GO_filtered <- circular_DF_GO %>% filter(!pfam %in% values_to_filter)
circular_DF_GO_filtered = circular_DF_GO_filtered %>% left_join(., go_categories, by = "GO_acc") 

# Find top N GO terms for each group
top_GO_terms_BP <- circular_DF_GO_filtered %>%
  filter(category == "BP") %>%
  filter(!is.na(GO_name)) %>%
  filter(GO_name != "hypothetical") %>%
  filter(GO_name != "") %>%
  dplyr::select(GO_name, id, Nucleotide_acc) %>%
  group_by(id, GO_name) %>%
  summarise(total = n(), .groups = "drop") %>%
  group_by(id) %>%
  arrange(desc(total)) %>%
  top_n(top_n, total) 

# Plot
GO_BP_duf95vsPks = ggplot(top_GO_terms_BP, aes(x = reorder(GO_name, total), y = total, fill = id)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "", y = "Count", fill = "Group") +
  theme_bw() +
  scale_fill_manual(values=c("lan" = "#d02727", "PF00365" = "goldenrod", "PF0090" = "#686a68", "PF01225" = "#4c9f7b")) +
  ggtitle(paste("Top", top_n, "GO: Biological Process"))
ggsave("/data/san/data0/users/david/intelligence/figures/LanM_GO_BP_v_pks.png", GO_BP_duf95vsPks, width = 20, height = 26, units = "cm")




# ~~~~~~~~~~~~~~~~~~~~~~~~
# correlational networking
# ~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("cooccur")
library(igraph)
library(ggraph)
library(cooccur)
library(reshape2)


matrix <- dcast(tablex, Nucleotide_acc~pfam)
matrix <- column_to_rownames(matrix, var = "Nucleotide_acc")

co <- cooccur(matrix, spp_names = TRUE)
matrix[1:5, 1:5]



# ~~~~~~~~~~~~~~~~~~~~~~~~
# correlational networking
# ~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("GGally")
library(GGally)

tablex_transformed <- tablex %>%
  dplyr::count(Nucleotide_acc, pfam) %>%
  pivot_wider(names_from = pfam, values_from = n, values_fill = list(n = 0))





# ~~~~~~~~~~~~~~~~~~~~~~~~
# Look for patterns of systems    1) PHAGE DEFENSE SYSTEMS
# ~~~~~~~~~~~~~~~~~~~~~~~~
ctrl_v_lanm = fread("/data/san/data1/users/david/intelligence_database/ctrl_v_lanm.tsv")


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
# Figure 4 
########################################
LanM_40k_window
defense_focus = c(#"NZ_CP019655.1",
"NZ_EQ973330.1","NZ_CAAIOF010000003.1","NZ_PIEU01000022.1","NZ_LN890331.1","NZ_PJTG01000002.1",
"NZ_LAWY01000037.1"
  )

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867")
dna_def_pfam = c("PF20473","PF01420","PF04851","PF02384","PF08463","PF08463","PF00270","PF01938","PF13649","PF07669","PF07669","PF20473",
  "PF02384","PF13337", "PF08665")
dna_def_prot = c("WP_110028332.1","WP_110028334.1")
library(dplyr)
library(tidyr)

fig_4_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% defense_focus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                stand = dir,
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

fig_4_plot <- gggenomes(fig_4_df) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3) +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("/data/san/data0/users/david/intelligence/figures/fig_4_plot.png", 
  plot = fig_4_plot, 
  width = 18, height = 10, units = "cm")

########################################
# Figure 5
########################################
LanM_40k_window
competence_focus = c("NC_018081.1",
"NZ_CNVF02000013.1",
"NZ_NJFO02000003.1",
"NZ_JH792105.1",
"NZ_FUXA01000008.1",
# "NC_009328.1",
"NZ_PIJH01000018.1",
"NZ_LMBZ01000008.1"
)

lanthipeptide_pfam = c("PF05147","PF13575","PF14867","PF03412","PF00082","PF16934","PF14867","PF07730")
competence_pfam = c("PF18146","PF00154","PF05952","PF06133","PF12072","PF02464","PF12072","PF05389")


fig_5_df  <- LanM_40k_window %>%
  filter(Nucleotide_acc %in% competence_focus) %>%
  mutate(color = "grey88") %>%
  dplyr::rename(seq_id = unique_window, 
                stand = dir,
                start = start,
                end = end,
                bin_id = species) %>%
  group_by(Protein_acc) %>%
  mutate(color = case_when(
    pfam %in% lanthipeptide_pfam ~ "goldenrod",
    pfam %in% competence_pfam ~ "red",
    Protein_acc %in% competence_focus ~ "red",
  )) %>%
  fill(color, .direction = "down1up") %>%  # Fill the color for all rows within each group
  ungroup() %>%
  mutate(color = ifelse(is.na(color), "grey88", color)) %>%
  mutate(bin_id = case_when(Nucleotide_acc== "NZ_LMBZ01000008.1" ~ "Solibacillus cecembensis",
    TRUE ~ bin_id))

fig_5_plot <- gggenomes(fig_5_df) %>%
  focus(pfam %in% lanthipeptide_pfam, .expand = c(10000,12000)) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3) +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("/data/san/data0/users/david/intelligence/figures/fig_5_plot.png", 
  plot = fig_5_plot, 
  width = 18, height = 10, units = "cm")


########################################
# Figure 6
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

fig_6_df  <- LanM_40k_window %>%
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

fig_6_plot <- gggenomes(fig_6_df) %>%
  focus(pfam %in% lanthipeptide_pfam, .expand = c(10000,12000)) +
  geom_seq() +
  geom_gene(aes(fill = color),
    size = 4) +
  scale_fill_identity() +
  coord_cartesian(clip = 'off') +
  scale_x_continuous(expand = c(0.2,0)) +
  geom_bin_label(size = 3) +
  geom_seq_label(size = 2, vjust = 1.5) +
  scale_x_continuous(expand=c(0.4,0.7,0.01,0.7)) +
  theme(legend.position="right") 

ggsave("/data/san/data0/users/david/intelligence/figures/fig_6_plot.png", 
  plot = fig_6_plot, 
  width = 18, height = 10, units = "cm")
