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


########################################
# Load in rodeo table for class II lanthipeptide
########################################
LanM_40k_window = fread("/data/san/data0/users/david/intelligence/lanthipeptide_rodeo_2020paper/LanM_paper_40k_window.tsv")
LanM_40k_window$unique_window <- paste(LanM_40k_window$Nucleotide_acc, LanM_40k_window$start_window, sep = "_")
LanM_phyla <- LanM_40k_window %>% 
  dplyr::select(phylum)


########################################
# Load in pfam descriptions
########################################
pfam_descs = fread("/data/san/data0/users/david/intelligence/tables/pfam_desc.tsv",
  col.names = c("pfam","CL","desc_CL","desc_pfam","name"))


########################################
# Load in control table and append taxonomy to it 
########################################
refseq_tsv_wPF00365_40kb = fread("/data/san/data0/users/david/intelligence/negative_dataset/databases/refseq_tsv_wPF00365_40kb_window.tsv")

## databases from mmseqs (clustering of all proteins within the refseq whole genome). 
refseq_pfam_annot <- fread("/data/san/data1/users/david/intelligence_database/rep_v_pfam_searchout.m8")
colnames(refseq_pfam_annot) <- c("locus_tag", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# left join this with the refseq_tsv_wPF00365_40kb_window.tsv
# and we can get pfam per protein in the window to mirror the rodeo output
PF00365_df <- left_join(refseq_tsv_wPF00365_40kb,
  refseq_pfam_annot, 
  by = "locus_tag",
  relationship = "many-to-many") 

# append taxonomy to the control dataset
PF00365_df <- append_taxonomy(PF00365_df, sqlFile) %>% 
  ungroup()

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
# Just get phylum that have Lans
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
  "/data/san/data1/users/david/intelligence_database/ctrl_v_LanM_vs_Pks.tsv")
write_tsv(ctrl_df, 
  "/data/san/data1/users/david/intelligence_database/LanM_ctrl_df_pfams_nucs.tsv")




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
ggsave(fig_1a, filename = "/data/san/data0/users/david/intelligence/figures/PF00365_control_window_protein_count_LanM_september23.png", 
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
ggsave(fig_1b, file="/data/san/data0/users/david/intelligence/figures/LanM_bacteriocin_window_protein_count.png",
  width = 6, height = 10, units = "cm")


# Figure 1
cowplot_figure_1 = cowplot::plot_grid(fig_1a, fig_1b, 
  labels = c("(a)", "(b)"), 
  ncol = 2,
  label_size = 10)

ggsave2("/data/san/data0/users/david/intelligence/figures/figure_1.png", 
  cowplot_figure_1, width = 12, height = 10, units = "cm", dpi=300)


