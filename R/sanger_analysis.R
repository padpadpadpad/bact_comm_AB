# having a go at analysing sanger sequence data
# using the .seq data rather than the .ab1 files

# load in packages - need to make sure these are installed
# have to install Blast+ separately
library(rBLAST)
library(taxize)
library(dplyr)
library(tidyr)
library(purrr)

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

# copy the .ab1 files of all the files where no blast hit was found into a new folder
no_blast_ab1 <- list.files('sanger/raw', pattern = '.ab1', full.names = TRUE)[tools::file_path_sans_ext(basename(list.files('sanger/raw', pattern = '.ab1', full.names = TRUE))) %in% tools::file_path_sans_ext(basename(no_blast_hit))]
file.copy(no_blast_ab1, 'sanger/No_blast_hit')

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

# alternative is to use the 6 consensus sequences from alignment in geneious ####

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

classification_df <- function(x, ...){
  temp <- taxize::classification(x, ...)
  temp_df <- data.frame(name = temp[[1]]$name, rank = temp[[1]]$rank, stringsAsFactors = FALSE)
  return(temp_df)
}
  
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
# Boom town
