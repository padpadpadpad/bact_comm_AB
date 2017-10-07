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
files <- list.files('sanger', pattern = '.seq', full.names = TRUE)

# load database
bl <- blast(db = "blast_16SMicrobial_db/16SMicrobial")

# get all alignments for all sequences
# keeps alignments that are over 97% in common
seq_IDs <- purrr::map_df(files, rBlast_all, database = bl, keep = 0.97) %>%
  mutate(., num = 1:n())

# We are in business!!!

# Find files which have no hit in Blast
no_blast_hit <- files[! basename(files) %in% c(seq_IDs$QueryID)]
length(no_blast_hit)
# 23 out of 96 have no hit on blast

# convert genbank ID into UID that can be accessed via genbank
# this takes a while as it calls the API

# filter the single top hit for each one. Still gives 290 rows as some are at the same identical percentage
seq_IDs <- group_by(seq_IDs, QueryID) %>%
  top_n(., 1, Perc.Ident) %>%
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

