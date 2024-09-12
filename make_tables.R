########################################
## This script creates the dataframes for the analysis
## However, this requires the user to have the Genbank in tsv format
## 
## 
## 
## Author: DEHourigan
########################################

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
## taxonomizeR - note this is a large download and install and requires to wget the database
########################################
# wget was used in terminal ## getAccession2taxid(baseUrl='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/')
library(taxonomizr)
setwd("/data/san/data0/users/david/taxonomizer") #SET UP DATABASE ONCE #
sqlFile <- "/data/san/data0/users/david/taxonomizer/accessionTaxa.sql"
# define the function to append taxonomy as columns
append_taxonomy <- function(df, sqlFile) {
  df$taxid <- taxonomizr::accessionToTaxa(df$Nucleotide_acc, sqlFile)
  taxa <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  for(taxon in taxa) { # append each taxonomic level
    df[[taxon]] <- getTaxonomy(df$taxid, sqlFile, desiredTaxa = taxon)
  }
  return(df)
}
########################################
## files required include
########################################


########################################
## files to be made include
########################################
# "LanM_paper_40k_window.tsv"
# "pfam_desc.tsv"
# "refseq_tsv_wPF00365_40kb_window.tsv"


########################################
## Make the rodeo output long format
# This needs to be done as there are cases where multiple pfams are on a single protein.
########################################
setwd("/data/san/data0/users/david/intelligence/code")
csv_files <- dir_ls(path = "data/rodeo_out/LanM_rodeo_out",
                    recurse = T ,
                    glob = "*occur.csv",
                    type = "file")
# read in tables
LanM <- map_dfr(csv_files, 
              ~ fread(.x, fill = TRUE))
# rename descriptors with dplyr
LanM <-  dplyr::rename(LanM, PfamID4 = V20)
LanM <-  dplyr::rename(LanM, Name4 = V21)
LanM <-  dplyr::rename(LanM, Description4 = V22)
LanM <-  dplyr::rename(LanM, 'E-value4' = V23)
LanM <-  dplyr::rename(LanM, PfamID5 = V24)
LanM <-  dplyr::rename(LanM, Name5 = V25)
LanM <-  dplyr::rename(LanM, Description5 = V26)
LanM <-  dplyr::rename(LanM, 'E-value5' = V27)

## call hypotheticals
LanM <- LanM %>%
  mutate(PfamID1=replace(PfamID1, PfamID1=="", "hypothetical")) %>%
  as.data.frame()
LanM <- LanM %>%
  mutate(PfamID2=replace(PfamID2, PfamID2=="", NA)) %>%
  as.data.frame()
LanM <- LanM %>%
  mutate(PfamID3=replace(PfamID3, PfamID3=="", NA)) %>%
  as.data.frame()
LanM <- LanM %>%
  mutate(PfamID4=replace(PfamID4, PfamID4=="", NA)) %>%
  as.data.frame()
LanM <- LanM %>%
  mutate(PfamID5=replace(PfamID5, PfamID5=="", NA)) %>%
  as.data.frame()

## Long format to gather all Pfams into one column to account for multiple pfams on a single protein
number_before_long <- length(unique(LanM$Protein_acc))
long <- data.table::melt(setDT(LanM), measure = patterns("PfamID.$", "Name.$","Description.$", "E-value.$"),
          value.name = c("pfam", "name", "desc", "eval"))[,
          variable:= NULL][order(Query)] 
write_tsv(long, file = "data/processed/LanM_long_format.tsv") ## write it to a table
number_after_long <- length(unique(long$Protein_acc))
message("before and after are equal:")
length(number_after_long) == length(number_before_long)  ## make sure they are equal and no data lost

########################################
## Append taxonomy using taxonomizer
########################################
lanM_long <- fread("data/processed/LanM_long_format.tsv") %>%
  as.data.frame() %>%
  mutate(lantype = "lanM") 
lan_long = as.data.table(lanM_long) 
lanM_long_tax = lanM_long %>% 
    append_taxonomy(., sqlFile) %>% 
    as_tibble() 

########################################
## read in walker et. al 2020 paper
########################################
library(readxl)
library(janitor)
file_path <- "data/rodeo_paper.xlsx"

# Read the sheets into a list and add a 'sheet_number' column
walker <- lapply(1:4, function(i) {
  df <- read_xlsx(file_path, sheet = i)
  df <- mutate(df, sheet_number = i)
  return(df)
})

# PF13575 & PF05147
lanM_paper_40k_window <- lanM_long_tax %>% 
  lazy_dt() %>%
  filter(Nucleotide_acc %in% walker[[2]]$nucleotide_acc) %>% 
  mutate(start_window = case_when(Protein_acc %in% dougL[[2]]$"LanM accession number" ~ start - 20000),
  end_window = case_when(Protein_acc %in% dougL[[2]]$"LanM accession number" ~ start + 20000),
  ) %>%   
  group_by(Nucleotide_acc) %>%
  arrange(desc(start_window)) %>% 
  fill(start_window, .direction = "downup") %>% 
  fill(end_window, .direction = "downup") %>% as_tibble() %>% 
  filter(case_when(start > start_window & start < end_window ~ T,
    TRUE ~ FALSE)) %>% 
  mutate(window_size = abs(start_window - end_window)) %>% 
  ungroup() %>%
  dplyr::select(-lantype, Query) %>%
  distinct() %>%
  as_tibble() 

########################################
## write LanM_paper_40k_window.tsv file
########################################
write_tsv(lanM_paper_40k_window, file = "data/processed/LanM_paper_40k_window.tsv")
lanM_paper_40k_window = fread("data/processed/LanM_paper_40k_window.tsv")



########################################
## Negative dataframe, PKS
########################################
# Blast results for control, these controls are strains in refseq with hits for
# Phosphofructokinase (PF00365) , Penicillin binding protein transpeptidase domain (PF00905)
# Mur ligase family, catalytic domain (PF01225)

paths <- c(
  "data/negative_dataset/blast_out/PF00365_refseq.blast",
  "data/negative_dataset/blast_out/PF00905_refseq.blast",
  "data/negative_dataset/blast_out/PF01225_refseq.blast"
) # only PF00365 is used in this analysis

df <- map_dfr(paths, ~ {
  data <- fread(.x, fill = TRUE)
  data$ID <- sub("_refseq.blast", "", basename(.x))
  data
})

colnames(df) <- c("query", "subject", "identity", "length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "qlen", "slen", "species", "ID")
control_df <- df %>%
  filter(length > qlen * 0.9 & evalue < 1e10 & identity > 40) %>%
  as.data.frame() %>%
  distinct()

write_tsv(control_df, file = "data/negative_dataset/regions/control_vs_refseq_filtered_blast.tsv")

## blast results for bacteriocin positive genomes in refseq
bac_p_noab <- fread("/data/san/data0/users/david/intelligence/negative_dataset/positive_bacteriocin_query_set_vs_refseq_blastn_filtered_rm_reg-abc-pep.txt")
colnames(bac_p_noab) <- c("query", "subject", "identity", "length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "qlen", "slen")
bac_p <- fread("/data/san/data0/users/david/intelligence/negative_dataset/positive_bacteriocin_query_set_vs_refseq_blastn.txt")
pep_locs <- c("/data/san/data0/users/david/intelligence/rodeo/rodeo_doug_all/main_results.csv")

########################################
## PF00365 control dataset (Pks)
########################################
PF00365_control_df <- control_df %>% 
  filter(ID == "PF00365", !is.na(evalue), !is.na(identity), !is.na(length)) %>%
  group_by(subject) %>%
  arrange(desc(evalue), desc(identity), desc(length)) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  dplyr::select(subject, s_start, s_end) %>%
  mutate(Nucleotide_acc = sub("\\..*$", "", subject)) %>%
  distinct()


########################################
# READ IN ALL REFSEQ TSV
########################################
library(parallel)
refseq_tsv_locs <- dir_ls(
  path = as.character("/data/san/data0/databases/genbank/refseq/bacteria"), # read in tables as list
  recurse = T,
  glob = "*.tsv",
  type = "file"
)

num_cores <- 16 
cl <- makeCluster(num_cores) # Create a cluster
clusterEvalQ(cl, library(data.table)) # Load the data.table package on each node
clusterExport(cl, c("refseq_tsv_locs")) # Export variables to the cluster

process_file <- function(file) {
  df <- fread(file, stringsAsFactors = FALSE, header = FALSE)
  df$id <- basename(file)
  return(df)
}
result_list_all <- parLapply(cl, refseq_tsv_locs, process_file) # Apply the function to each file in parallel
stopCluster(cl) # Stop the cluster

refseq_tsv <- do.call(rbind, result_list_all) # Combine the results
colnames(refseq_tsv) <- c("locus_tag", "type", "gene_name", "ec_num", "product", "contig_id", "strain_name", "Nucleotide_acc", "prot_seq", "start", "end", "ID")

refseq_tsv %>%
  mutate(
    filename = str_extract(ID, "/[^/]+$"), # extract filename
    filename = str_remove(filename, "/"), # remove leading "/"
    filename = str_extract(filename, "^[^_]*_[^_]*") # extract everything before 2nd "_"
  )


refseq_tsv_wPF00365 <- left_join(PF00365_control_df, refseq_tsv , by = "Nucleotide_acc", relationship = "many-to-many") # append


# Write the regions around PF00365 to a file
refseq_tsv_wPF00365_40kb <-  refseq_tsv_wPF00365 %>%
  distinct() %>%
  group_by(Nucleotide_acc) %>%
  filter(start >= (s_start-20000) & start <= (s_end+20000))
write_tsv(refseq_tsv_wPF00365_40kb, "data/negative_dataset/refseq_tsv_wPF00365_40kb_window.tsv")
 

refseq_tsv_wPF00365_40kb = fread("data/negative_dataset/refseq_tsv_wPF00365_40kb_window.tsv")


########################################
# Load in pfam descriptions
########################################
pfam_descs = fread("tables/pfam_desc.tsv",
  col.names = c("pfam","CL","desc_CL","desc_pfam","name"))


########################################
# Load in control table and append taxonomy to it 
########################################
refseq_tsv_wPF00365_40kb = fread("data/processed/refseq_tsv_wPF00365_40kb_window.tsv")

## databases from mmseqs (clustering of all proteins within the refseq whole genome). 
refseq_pfam_annot <- fread("data/negative_dataset//rep_v_pfam_searchout.m8")
colnames(refseq_pfam_annot) <- c("locus_tag", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# left join this with the refseq_tsv_wPF00365_40kb_window.tsv
# and we can get pfam per protein in the window to mirror the rodeo output
PF00365_df <- left_join(refseq_tsv_wPF00365_40kb,
  refseq_pfam_annot, 
  by = "locus_tag",
  relationship = "many-to-many") 

# append taxonomy to the control dataset
PF00365_df <- append_taxonomy(PF00365_df, sqlFile) %>%
  ungroup() %>%
  mutate(Nucleotide_acc = subject.x)

fwrite(PF00365_df, "data/negative_dataset/PF00365_df.tsv", sep = "\t")
PF00365_df <- fread("data/negative_dataset/PF00365_df.tsv", sep = "\t") 
