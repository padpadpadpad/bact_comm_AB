# having a go at analysing sanger sequence data
# using the .seq data rather than the .ab1 files

# clear workspace ####
# mise::mise(vars = TRUE, console = TRUE, figs = TRUE, pkgs = TRUE)

# load in packages - need to make sure these are installed
# have to install Blast+ separately
library(rBLAST)
library(taxize)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# functions ####

# write function to do cool shit - allows for calling of blast from R on multiple sequences at once
rBlast_all <- function(seq_files, database, keep = 0.95){
  
  test <- read.table(seq_files, blank.lines.skip = TRUE, as.is = TRUE)
  seq1 <- paste(test$V1, collapse = '')
  seq2 <- BStringSet(seq1)
  names(seq2) <- basename(seq_files)
  
  seq_pred <- predict(database, seq2) %>%
    dplyr::filter(., Perc.Ident > keep*100) %>%
    mutate_at(., c('QueryID', 'SubjectID'), as.character) %>%
    mutate(seq = seq1)
    
  return(seq_pred)
}

# function to clean up output from assignTaxonomy
compile_taxonomy <- function(assignTaxonomy_object, tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")){
  temp <- data.frame(unname(assignTaxonomy_object), stringsAsFactors = FALSE)
  names(temp) <- tax_levels
  return(temp)
}

# a function to return the longest sequence that only contains CATG
clean_base_pairs <- function(sequence){
  temp <- stringr::str_split(sequence, "[^CATG]")
  temp <- data.frame(seq = unlist(temp), stringsAsFactors = FALSE)
  temp$length = nchar(temp$seq)
  return(temp[temp$length == max(temp$length),]$seq)
}

# return dataframe from taxize::classification
classification_df <- function(x, ...){
  temp <- taxize::classification(x, ...)
  temp_df <- data.frame(name = temp[[1]]$name, rank = temp[[1]]$rank, stringsAsFactors = FALSE)
  return(temp_df)
}

#download databases - ONLY NEED TO DO THIS ONCE
#download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz","16SMicrobial.tar.gz", mode='wb')
#download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz", "taxdb.tar.gz", mode='wb')
#untar("16SMicrobial.tar.gz", exdir="blast_16SMicrobial_db")
#untar("taxdb.tar.gz", exdir="blast_16SMicrobial_db")

# list files in the data folder
files <- list.files('sanger/trimmed_sanger', pattern = '.txt', full.names = TRUE)

# load database
bl <- blast(db = "sanger/blast_16SMicrobial_db/16SMicrobial")

# get all alignments for all sequences
# keeps alignments that are over 97% in common
seq_IDs <- purrr::map_df(files, rBlast_all, database = bl, keep = 0.97) %>%
  mutate(., num = 1:n())

# We are in business!!!

# Find files which have no hit in Blast
no_blast_hit <- files[! basename(files) %in% c(seq_IDs$QueryID)]
length(no_blast_hit)
# 23 out of 96 have no hit on blast

# copy the .ab1 files of all the files where no blast hit was found into a new folder - My project does not like copying files anymore...
no_blast_ab1 <- list.files('sanger/raw', pattern = '.ab1', full.names = TRUE)[tools::file_path_sans_ext(basename(list.files('sanger/raw', pattern = '.ab1', full.names = TRUE))) %in% tools::file_path_sans_ext(basename(no_blast_hit))]
# file.copy(no_blast_ab1, 'sanger/No_blast_hit', overwrite = TRUE)

# convert genbank ID into UID that can be accessed via genbank
# this takes a while as it calls the API

# filter the single top hit for each one. Still gives 290 rows as some are at the same identical percentage
seq_IDs <- group_by(seq_IDs, QueryID) %>%
  arrange(., Perc.Ident) %>%
  top_n(., 1, row_number()) %>%
  ungroup() %>%
  mutate(num = 1:n())

# work out the uid from the genbank id - given as SubjectID
# This has proved problematic as we are trying to deal with errors
# approach - nest each row in a list and call the API separately - return NA if it does not work
seq_IDs_nest <- seq_IDs %>%
  nest(., SubjectID) %>%
  mutate(., uid = map(data, possibly(genbank2uid, otherwise = NA, quiet = TRUE)))

# unnest these columns to get the dataframes out. There shall be warnings but they are ok as the correct uids are preserved
seq_IDs <- seq_IDs_nest %>%
  unnest(data, uid)

# now nest uids of each row and run the get_taxon_summary API call on each list element
seq_IDs_nest <- seq_IDs %>%
  nest(., uid) %>%
  mutate(., tax_info = map(data, possibly(ncbi_get_taxon_summary, otherwise = NA, quiet = TRUE)))

# unnest this dataframe for final dataframe of what everything is
seq_IDs <- seq_IDs_nest %>%
  unnest(tax_info) %>%
  select(., -data)

# BINGO we are in business

# keep columns that we want...
select(seq_IDs, QueryID, Perc.Ident, name) %>%
  write.csv(., 'sanger/processed.csv', row.names = FALSE)

check <- select(seq_IDs, QueryID, Perc.Ident, name) %>%
  arrange(., name)

# alternative is to use the 7 consensus sequences from alignment in geneious ####

# load in consensus sequences
consensus <- read.table('sanger/consensus/Consensus_sequences_alignment_all.fasta', stringsAsFactors = FALSE)

# get names
labels <- consensus$V1[grepl('>', consensus$V1)]
labels <- gsub('>', '', labels)

# get sequences
seqs <- consensus$V1[!grepl('>', consensus$V1)]

# change '-' to N for noise
seqs <- gsub('-', 'N', seqs)

# save new files 
for(i in 1:length(labels)){
  write(seqs[i], paste('sanger/consensus', '/', labels[i], '.seq', sep = ''))
}

# list these new files
consensus_list <- list.files('sanger/consensus', pattern = '.seq$', full.names = TRUE)

# run BLAST
consensus_IDs <- purrr::map_df(consensus_list, rBlast_all, database = bl, keep = 0.97) %>%
  mutate(., num = 1:n())

# filter the single top hit for each one. Still gives 290 rows as some are at the same identical percentage
consensus_IDs <- group_by(consensus_IDs, QueryID) %>%
  top_n(., 1, Perc.Ident) %>%
  ungroup() %>%
  mutate(num = 1:n())

# work out the uid from the genbank id - given as SubjectID
consensus_IDs_nest <- consensus_IDs %>%
  nest(., SubjectID) %>%
  mutate(., uid = map(data, possibly(genbank2uid, otherwise = NA, quiet = TRUE)))

# unnest these columns to get the dataframes out. There shall be warnings but they are ok as the correct uids are preserved
consensus_IDs <- consensus_IDs_nest %>%
  unnest(data, uid)

# now nest uids of each row and run the get_taxon_summary API call on each list element
consensus_IDs_nest <- consensus_IDs %>%
  nest(., uid) %>%
  mutate(., tax_info = map(data, possibly(ncbi_get_taxon_summary, otherwise = NA, quiet = TRUE)))

# unnest this dataframe for final dataframe of what everything is
consensus_IDs <- consensus_IDs_nest %>%
  unnest(tax_info) %>%
  select(., -data)

# BINGO we are in business

# filter consensus_IDs
consensus_IDs <- group_by(consensus_IDs, QueryID) %>%
  distinct(., name, .keep_all = TRUE) %>%
  ungroup() %>%
  select(., QueryID, Perc.Ident, name) %>%
  mutate(., num = row_number())

write.csv(consensus_IDs, 'sanger/consensus_processed.csv', row.names = FALSE)

# get higher taxonomy using taxize ####
consensus_IDs <- consensus_IDs %>%
  nest(name) %>%
  mutate(., higher_tax = map(data, possibly(classification_df, otherwise = NA, quiet = TRUE), db = 'ncbi')) %>%
  unnest(higher_tax)

# look at taxonomy levels
unique(consensus_IDs$rank)

# ranks to keep 
tax_rank <- c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

consensus_IDs <- filter(consensus_IDs, rank %in% tax_rank) %>%
  spread(., rank, name) %>%
  select(., num, QueryID, Perc.Ident, species, genus, family, order, class, phylum, superkingdom) %>%
  data.frame()

write.csv(consensus_IDs, 'sanger/consensus_processed_higher_taxonomy.csv', row.names = FALSE)

# check for differences between consensus and raw
spp_raw <- unique(seq_IDs$name)
spp_consensus <- unique(consensus_IDs$species)

spp_raw
spp_consensus

# assign Taxonomy in R using dada2 ####
# use dada2 to assign taxonomy - also see attempt in scripts/dada2_play.R

# remove 'N' from seqs - assignTaxonomy DOES NOT LIKE 'N's
seqs2 <- map(seqs, clean_base_pairs) %>%
  unlist()

# rdp_taxonomy  
rdp_taxa <- dada2::assignTaxonomy(seqs2, "sanger/ref_trainsets/rdp_train_set_16.fa", multithread = TRUE)
rdp_taxa <- dada2::addSpecies(rdp_taxa, "sanger/ref_trainsets/rdp_species_assignment_16.fa", allowMultiple = TRUE)
rdp_taxa <- compile_taxonomy(rdp_taxa) %>%
  mutate(., db = 'rdp', seq = 1:7)

# green genes_taxonomy
gg_taxa <- dada2::assignTaxonomy(seqs2, "sanger/ref_trainsets/gg_13_8_train_set_97.fa", multithread = TRUE)
gg_taxa <- compile_taxonomy(gg_taxa) %>%
  mutate(., db = 'gg',
         seq = 1:7)

# silva taxonony
silva_taxa <- dada2::assignTaxonomy(seqs2, "sanger/ref_trainsets/silva_nr_v128_train_set.fa", multithread = TRUE)
silva_taxa <- dada2::addSpecies(silva_taxa, "sanger/ref_trainsets/silva_species_assignment_v128.fa", allowMultiple = TRUE)
silva_taxa <- compile_taxonomy(silva_taxa) %>%
  mutate(., db = 'silva',
         seq = 1:7)

# bind together
consensus_bind <- bind_rows(silva_taxa, rdp_taxa, gg_taxa)
write.csv(consensus_bind, 'sanger/consensus_multi_db.csv', row.names = FALSE)


# read in each sequence
trimmed_seqs <- map_df(files, read.table, blank.lines.skip = TRUE, as.is = TRUE) %>%
  rename(., seq = V1) %>%
  mutate(file = tools::file_path_sans_ext(basename(files)),
         length = nchar(seq))

# find longest portion of each sequence without an N = Errors
trimmed_seqs <- trimmed_seqs %>%
  nest(., seq) %>%
  mutate(., clean_seq = map(data, clean_base_pairs)) %>%
  unnest(data) %>%
  unnest(clean_seq) %>%
  mutate(length_clean_seq = nchar(clean_seq)) %>%
  # keep only sequences longer than 100 base pairs
  filter(length_clean_seq > 100)

# assign taxonomy on all the samples
rdp_taxa_all <- dada2::assignTaxonomy(trimmed_seqs$clean_seq, "sanger/ref_trainsets/rdp_train_set_16.fa", multithread = TRUE)
rdp_taxa_all <- dada2::addSpecies(rdp_taxa_all, "sanger/ref_trainsets/rdp_species_assignment_16.fa", allowMultiple = 2)
rdp_taxa_all <- compile_taxonomy(rdp_taxa_all, tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
rdp_taxa_all <- cbind(rdp_taxa_all, select(trimmed_seqs, file, length_clean_seq))

write.csv(rdp_taxa_all, 'sanger/rdp_all_taxonomy.csv', row.names = FALSE)
