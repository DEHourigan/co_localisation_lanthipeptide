########################################
## Packages - don't package shame people
########################################
.libPaths(c("/data/san/data0/users/david/rstudio/packages", .libPaths()))
newlib <- "/data/san/data0/users/david/rstudio/packages"
packages <- c("cowplot","data.table", "multidplyr", "readr", "formattable", "fs", "dplyr", "ggplot2", "purrr", "ggthemes", "BiocManager", "gplots", "gridExtra", "grid", "forcats", "tidyr", "dtplyr", "topGO", "SparseM", "Biostrings", "GenomicRanges", "seqinr", "stringr", "readxl", "thacklr", "gggenomes","ggrepel")

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
# function to append taxonomy given a nucleotide accession
########################################
library(taxonomizr)
# wget was used in terminal ## getAccession2taxid(baseUrl='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/')
# setwd("/data/san/data0/users/david/taxonomizer") #SET UP DATABASE ONCE #
sqlFile <- "/data/san/data0/users/david/taxonomizer/accessionTaxa.sql"
  
append_taxonomy_control <- function(df, sqlFile) {
  df$taxid <- taxonomizr::accessionToTaxa(df$Nucleotide_acc, sqlFile)
  taxa <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  for(taxon in taxa) {
    df[[taxon]] <- getTaxonomy(df$taxid, sqlFile, desiredTaxa = taxon)
  }
  return(df)
}

append_taxonomy <- function(df, sqlFile) {
  df$taxid <- taxonomizr::accessionToTaxa(df$Nucleotide_acc, sqlFile)
  taxa <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  for(taxon in taxa) {
    df[[taxon]] <- getTaxonomy(df$taxid, sqlFile, desiredTaxa = taxon)
  }
  
  return(df)
}

########################################
## Control dataset taxonomy
########################################
control_x <- fread("/data/san/data2/users/david/co_localisation/data/processed/LanM_ctrl_df_pfams_nucs.tsv")
control_x <- append_taxonomy_control(control_x, sqlFile)
control_x <- control_x %>% dplyr::select(-pfam,id) %>%
	filter(group == "control") %>%
	distinct()

# plot the number of unique Nucleotide accs per phylum
control_x %>%
	group_by(genus) %>%
	summarise(n = n()) %>%
	arrange(desc(n)) %>%
	top_n(10, n) %>%
	ggplot(aes(x = reorder(genus, n), y = n)) +
	geom_bar(stat = "identity") +
	coord_flip() +
	labs(title = "Top 10 Genera by Number of Unique Nucleotide Accs",
		x = "Genus",
		y = "Number of Unique Nucleotide Accs") +
	theme_classic(base_size = 8) +
	theme(legend.position = "none")

########################################
pfam_descs <- fread("tables/pfam_desc.tsv", col.names = c("pfam", "clan","clan_name","desc_pfam","name"))
def = fread("tables/phage_defense.txt", header=FALSE, 
  col.names = c("pfam")) %>%
  left_join(., pfam_descs, by = c("pfam" = "pfam")) %>%
  filter(!grepl("CL0023", clan)) # remove general ATPase




LanM_40k_window = fread("data/processed/LanM_paper_40k_window.tsv")
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
LanM_phyla <- LanM_40k_window %>% 
  dplyr::select(phylum)
LanM_contigs <- LanM_40k_window %>%
  dplyr::select(Nucleotide_acc) %>%
  distinct()


lanM_long <- fread("data/processed/LanM_long_format.tsv") %>%
  as.data.frame() %>%
  mutate(lantype = "lanM") 
colnames(lanM_long) 


LanM_pfam = c("PF05147","PF13575")
transport_pfams = c("PF00005","PF07690","PF00664","PF02687")
peptidase_pfams = c("PF03412", "PF00082")
core_pfams = c("PF16934", "PF04604", "PF14867")


new_nucleotide_acc = lanM_long %>%
	filter(!Nucleotide_acc %in% LanM_contigs$Nucleotide_acc) %>%
	group_by(Nucleotide_acc, Protein_acc) %>%
	# make sure each group contains the positive_pfams
	filter(all(LanM_pfam %in% pfam)) %>%
	dplyr::select(Nucleotide_acc) %>%
	distinct()

new_nucleotide_df = lanM_long %>%
	filter(Nucleotide_acc %in% new_nucleotide_acc$Nucleotide_acc) %>%
	group_by(Nucleotide_acc) %>%
	filter(
		any(pfam %in% transport_pfams) &
		any(pfam %in% peptidase_pfams) &
		any(pfam %in% core_pfams)
	) %>%
	ungroup()

n_distinct(LanM_40k_window$Nucleotide_acc)
n_distinct(new_nucleotide_df$Nucleotide_acc)

percentages_table = new_nucleotide_df %>%
	dplyr::select(Nucleotide_acc, pfam, desc) %>%
	distinct() %>%
	group_by(pfam, desc) %>%
	summarise(n = n()) %>%
	arrange(desc(n)) %>%
	# mutate the percentage of total Nucleotide accs that have each pfam
	mutate(perc = n / n_distinct(new_nucleotide_df$Nucleotide_acc) * 100) 

new_nucleotide_pfam = new_nucleotide_df %>%
	dplyr::select(Nucleotide_acc, pfam, desc) %>%
	distinct()

# Get all the pfams that are acetyltransferases
acetyltransferase_pfams <- new_nucleotide_df %>%
	filter(grepl("Acetyltransferase", desc)) %>%
	dplyr::select(pfam) %>%
	distinct() %>%
	pull(pfam) %>%
	unique()

# plot bar chart of percentages percentages_table for top 20 pfams
percentages_table %>%
	filter(pfam %in% acetyltransferase_pfams) %>%
	ggplot(aes(x = reorder(pfam, perc), y = perc)) +
	geom_bar(stat = "identity") +
	coord_flip() +
	labs(
		title = "Percentage of Nucleotide Accs with Acetyltransferase Pfams",
		x = "Pfam",
		y = "Percentage of Nucleotide Accs"
	)
ggsave("figures/figure_SX_Aceyltransferase_validation.png", 
	width = 12, height = 10, units = "cm", dpi = 300)


cleaned_data <- new_nucleotide_pfam %>%
  filter(!is.na(pfam) & pfam != "")

pfam_pairs <- cleaned_data %>%
  group_by(Nucleotide_acc) %>%
  summarise(pfam_combinations = list(combn(pfam, 2, simplify = FALSE))) %>%
  unnest(pfam_combinations)


# Count co-occurrence of each Pfam pair
co_occurrence_counts <- pfam_pairs %>%
  unnest_wider(pfam_combinations, names_sep = "_") %>%
  group_by(Pfam1 = pfam_combinations_1, Pfam2 = pfam_combinations_2) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

# Randomize Pfams within Nucleotide_acc
set.seed(123)
randomized_data <- cleaned_data %>%
  group_by(Nucleotide_acc) %>%
  mutate(shuffled_pfam = sample(pfam)) %>%
  ungroup()

randomized_pairs <- randomized_data %>%
  group_by(Nucleotide_acc) %>%
  summarise(random_combinations = list(combn(shuffled_pfam, 2, simplify = FALSE))) %>%
  unnest(random_combinations)

randomized_counts <- randomized_pairs %>%
  unnest_wider(random_combinations, names_sep = "_") %>%
  group_by(Pfam1 = random_combinations_1, Pfam2 = random_combinations_2) %>%
  summarise(random_count = n(), .groups = "drop")

comparison <- co_occurrence_counts %>%
  left_join(randomized_counts, by = c("Pfam1", "Pfam2")) %>%
  mutate(random_count = replace_na(random_count, 0)) %>%
  mutate(p_value = (random_count + 1) / (sum(random_count, na.rm = TRUE) + 1)) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

significant_pairs <- comparison %>%
  filter(adj_p_value < 0.05) %>%
  filter(Pfam1 %in% core_pfams) %>%
  left_join(., pfam_descs, by = c("Pfam2" = "pfam")) 

significant_pairs %>%
	filter(Pfam2 %in% acetyltransferase_pfams) %>% View

plot_counts_vs_random <- ggplot(significant_pairs,
	aes(
		x = random_count, 
		y = count, 
		fill = factor(adj_p_value < 0.05, levels = c(FALSE, TRUE), labels = c("Not Significant", "Significant")))) +
	geom_point(
		alpha = 0.3,
		shape = 21,
		color = "black") +
	labs(title = "Scatter Plot: Observed vs. Randomized Co-Occurrence Counts", x = "Observed Count", y = "Randomized Count") +
	theme_classic(base_size = 8) +
	theme(legend.title = element_blank(),
		legend.position = c(0.25, 0.80)) +
	geom_text_repel(aes(label = ifelse(Pfam2 %in% acetyltransferase_pfams, desc_pfam, "")),
		vjust = -0.5,
		hjust = 1,
		size = 2.5,
		max.overlaps = getOption("ggrepel.max.overlaps", default = 25),
		color = "black") +
	geom_text_repel(aes(label = ifelse(Pfam2 %in% def$pfam, desc_pfam, "")),
		vjust = -0.5,
		hjust = 1,
		size = 2.5,
		max.overlaps = getOption("ggrepel.max.overlaps", default = 25),
		color = "red") +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")

ggsave(plot_counts_vs_random, 
  file="figures/figure_SX_validation.png",
  width = 12, height = 10, units = "cm", dpi = 300)





bacteriocin_pfams = c("PF14867", "PF13537", "PF13575", "PF12730", "PF12730", "PF02441", "PF01944",
                 "PF18218", "PF01320", "PF00072", "PF00486", "PF00005", "PF02518", "PF00512",
                 "PF03857", "PF08951", "PF05147", "PF00881", "PF06182", "PF00296", "PF05402",
                 "PF15565", "PF15586", "PF01719", "PF02794", "PF04055", "PF13165", "PF15007",
                 "PF11083", "PF13471", "PF00733", "PF02540", "PF05147", "PF02052", "PF04604",
                 "PF08130")

  
genomad_df = fread("/data/san/data0/databases/genomad_v1.7/genomad_db/genomad_marker_metadata.tsv")
genomad_df$ANNOTATION_ACCESSIONS %>% unique()
genomad_pfams <- unique(unlist(str_extract_all(genomad_df$ANNOTATION_ACCESSIONS, "PF\\d{5}")))

dram_amg_df = fread("/data/san/data0/databases/DRAM/20240927/amg_database.20240927.tsv")
dram_amg_df$PFAM



########################################
## Liklihood of finding a defense system in a randomised dataset
########################################
set.seed(420)
input_dir <- "/data/san/data2/users/david/co_localisation/data/contigs_control_files"
output_dir <- "/data/san/data2/users/david/co_localisation/data/contigs_control_sampled"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
sample_size <- 30
iterations <- 200

sample_consecutive_rows <- function(df, n) {
  if (nrow(df) < n) {
    stop("Data frame has fewer rows than the sample size.")
  }
  start_index <- sample(1:(nrow(df) - n + 1), 1)
  sampled_df <- df[start_index:(start_index + n - 1), , drop = FALSE]
  return(sampled_df)
}

for (i in 1:iterations) {
  file <- sample(tsv_files, 1)
  df <- fread(file, header = FALSE, sep = "\t")
  sampled_df <- sample_consecutive_rows(df, sample_size)
  base_name <- gsub(".*/|\\.tsv$", "", file)
  output_file <- file.path(output_dir, paste0(base_name, "_sampled_", i, ".tsv"))
  fwrite(sampled_df, output_file, sep = "\t", col.names = FALSE, row.names = FALSE)
}

########################################
## Liklihood of finding a defense system 40kb window of lanthipeptide
########################################
# the below command was used to get the defense systems in the random dataset
# cut -f1 -d' ' *  | grep -f - ../defense_control/*/*systems*.tsv  -h  | cut -f1,2 -d_ | sort | uniq | wc | cut -f1
# Proportions of occurrence in each dataset

iteration_occurrences <- 6 / 200
test_dataset_occurrences <- 186 / 1412


observed <- matrix(c(194, 6, 1226, 186), nrow = 2, byrow = TRUE)

test_result <- fisher.test(observed)


proportions <- c(iteration_occurrences, test_dataset_occurrences)
names(proportions) <- c("Random 40kb", "Lanthipeptide 40kb")

proportions_df <- data.frame(
	Dataset = c("Random 40kb", "Lanthipeptide 40kb"),
	Proportion = proportions)

random_40kb_region <- ggplot(proportions_df, aes(x = Dataset, y = Proportion, fill = Dataset)) +
	geom_bar(stat = "identity", color = "black") +
	scale_fill_manual(values = c("Random 40kb" = "goldenrod", "Lanthipeptide 40kb" = "red")) +
	labs(y = "Proportion", x = "Dataset") +
	theme_classic() +
	theme(legend.position = "none") +
	ylim(0, 0.20) +
	annotate("text", x = 1.5, y = 0.2, label = paste0("Fisher's p-value: ", signif(test_result$p.value, 3)), size = 4) +
	annotate("text", x = 1, y = test_dataset_occurrences + 0.01, label = "~5x", size = 4)

ggsave(random_40kb_region, 
	file = "figures/figure_SX_random_40kb.png",
	width = 10, height = 8, units = "cm", dpi = 600)


########################################
# Liklihood of seeing a pfam enriched by chance
########################################
set.seed(420)

# read in enriched pfams from Pfk
control_lanthipeptide_pfams <- fread("tables/Supplementary_table_S1.csv")
lanthipeptide_pfams = fread("data/processed/ctrl_v_LanM_vs_Pks.tsv") %>%
	filter(grepl("^PF", pfam)) %>%
	dplyr::select(pfam) %>%
	distinct()

pfam_counts <- control_lanthipeptide_pfams %>% 
  dplyr::select(pfam, genes_in_lanthipeptide) %>%
  filter(pfam %in% lanthipeptide_pfams$pfam)

# Generate random counts by sampling PFAMs with replacement
random_counts <- table(sample(pfam_counts$pfam, sum(pfam_counts$genes_in_lanthipeptide), replace = TRUE))


pfam_counts <- merge(
  pfam_counts, 
  data.table(pfam = names(random_counts), expected = as.integer(random_counts)), 
  by = "pfam", 
  all.x = TRUE
)

pfam_counts[, expected := fifelse(is.na(expected), 0, expected)]
pfam_counts[, sd_random := sqrt(expected)]
pfam_counts[, Z_score := (genes_in_lanthipeptide - expected) / sd_random]
pfam_counts[, p_value := 2 * pnorm(-abs(Z_score))]
pfam_counts[, bonferroni_p_value := p.adjust(p_value, method = "bonferroni")]
pfam_counts[, significant := bonferroni_p_value < 0.05]
pfam_counts[, odds_ratio := (genes_in_lanthipeptide / expected)]

pfam_counts_vs_chance = left_join(pfam_counts, pfam_descs, by = c("pfam" = "pfam"))
fwrite(pfam_counts_vs_chance,
	 "tables/pfam_counts_vs_chance.tsv", sep = "\t", row.names = FALSE)

# plot 
random_pfams_histogram = ggplot(pfam_counts_vs_chance, aes(x = expected)) +
	geom_histogram(binwidth = 1, fill = "red", color = "black", alpha = 0.7) +
		labs(
		x = "Expected Pfam Counts",
		y = "Frequency") +
	theme_classic(base_size = 8) +
	theme(legend.position = "none")


ggsave(random_pfams_histogram, 
	file = "figures/figure_SX_random_pfams_histogram.png",
	width = 8, height = 6, units = "cm", dpi = 600)

########################################
# Enriched functions DRAM and mobileOGdb
########################################
mobileogdb = fread("/data/san/data2/users/david/co_localisation/mobileOG-db_workdir/pfam_assignments.tsv",
	col.names = c("name", "desc_pfam", "pfam", "evalue", "score")) %>%
	filter(evalue < 1e-30)
dram = fread("/data/san/data0/databases/DRAM/20240927/amg_database.20240927.tsv")
doron = fread("tables/doron_phage_defense_pfams.tsv")
pfk_enriched = fread("tables/Supplementary_table_S2.tsv") %>%
	filter(adjusted_p_value < 0.05)

pfam_counts_vs_chance_or2 <- pfam_counts_vs_chance %>%
	filter(significant == TRUE & odds_ratio > 2) %>%
	mutate(
		in_mobileogdb = sapply(pfam, function(x) any(grepl(x, mobileogdb$pfam))),
		in_dram = sapply(pfam, function(x) any(grepl(x, dram$PFAM))),
		in_doron = sapply(pfam, function(x) any(grepl(x, doron$Family))),
		in_pfk_enriched = sapply(pfam, function(x) any(grepl(x, pfk_enriched$pfam))),
	)

# count significant pfams in each category
pfam_tally <- pfam_counts_vs_chance_or2 %>%
	summarise(
		count_in_mobileogdb = sum(in_mobileogdb == TRUE, na.rm = TRUE),
		count_in_dram = sum(in_dram == TRUE, na.rm = TRUE),
		count_in_doron = sum(in_doron == TRUE, na.rm = TRUE)
	) %>%
	pivot_longer(
		cols = everything(), 
		names_to = "category", 
		values_to = "count"
	) %>%
	mutate(
		main_category = case_when(
			str_detect(category, "mobileogdb") ~ "MobileOG-db",
			str_detect(category, "dram") ~ "DRAM",
			str_detect(category, "doron") ~ "Doron et al."
		)
	)

# plot the counts of pfams in each category
bar_plot_enriched_functions <- pfam_tally %>%
	ggplot(aes(x = main_category, y = count, fill = main_category)) +
	geom_bar(stat = "identity", color = "black", alpha = 0.7) +
	labs(
		x = "Category",
		y = "Count",
		fill = "Category"
	) +
	scale_fill_manual(values = c(
		"MobileOG-db" = "goldenrod",
		"DRAM" = "grey88",
		"Doron et al." = "red"
	)) +
	theme_classic(base_size = 8) +
	theme(legend.position = "none")

ggsave(bar_plot_enriched_functions,
	file = "figures/figure_SX_enriched_functions.png",
	width = 8, height = 6, units = "cm", dpi = 600)

# count in_pfk_enriched = TRUE 
pfam_counts_vs_chance_or2 %>%
	filter(in_pfk_enriched == TRUE) 

########################################
# Volcano plot
########################################
pfam_labels_2 <- c(
  "Competence-damaged protein", 
  "Methylase_S",
  "Methyltransf_25", 
  "Lantibiotic_a",
  "Peptidase_S8",
  "AbiJ_NTD3", 
  "MqsA_antitoxin", 
  "Phage_integrase",
  "EcoR124_C",
  "rve",
  "DEAD",
  "CinA",
  "ResIII",
  "Recombinase",
  "MobC_2",
  "Zeta_toxin",
  "ParE_toxin",
  "Methyltransf_12",
  "Bacteriocin_IIc",
  "ABC_tran",
  "Acetyltransf_6",
  "Acetyltransf_7",
  "rve_3",
  "Eco57I")

pfam_counts_10_pfam_counts <- pfam_counts_vs_chance %>% 
	mutate(significant = bonferroni_p_value < 0.05) %>% # Ensure 'significant' column exists
	filter(genes_in_lanthipeptide > 10)

# Create a volcano plot with odds ratio on x-axis and -log10(p-value) on y-axis
pfam_volcano_chance <- pfam_counts_vs_chance %>%
	ggplot(aes(x = log2(odds_ratio), y = -log10(p_value), fill = ifelse(significant & odds_ratio > 1, "red", "grey88"))) +
	geom_point(alpha = 0.5, shape = 21, color = "black") +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
	geom_vline(xintercept = log2(2), linetype = "dashed", color = "black") +
	ggrepel::geom_text_repel(aes(label = ifelse(significant & desc_pfam %in% pfam_labels_2, desc_pfam, "")), 
		size = 2.5, max.overlaps = Inf, color = "black", nudge_x = 0.1, nudge_y = 0.1, force = 10) +
	labs(
		x = "Log2(Odds Ratio)",
		y = "-log10(p-value)") +
	scale_fill_identity() +
	theme_classic(base_size = 8) +
	theme(legend.position = "none") +
	ggside::ggside(x.pos = "top", y.pos = "right") +
	ggside::geom_xsidedensity(aes(y = after_stat(density)), alpha = 0.5, fill = "red") +
	ggside::geom_ysidedensity(aes(x = after_stat(density)), alpha = 0.5, fill = "red") +
	coord_cartesian(clip = "off") +
	xlim(-1, 10) 

ggsave(
	pfam_volcano_chance,
	file = "figures/figure_SX_pfam_counts_vs_chance.png",
	width = 14, height = 10, units = "cm", dpi = 300)