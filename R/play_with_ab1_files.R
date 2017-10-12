# looking at the .ab1 files - trim and clean some of these files, then save the sequence as a .seq file

# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite(c('DECIPHER', 'sangerseqR', 'ShortRead', 'msa'))

# load packages
library(sangeranalyseR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(Biostrings)
library(stringr)
library(ShortRead)
library(msa)

# Primers used
# 515F - GTGYCAGCMGCCGCGGTAA
# 806R - GGACTACNVGGGTWTCTAAT
fwd_515F <- 'GTGYCAGCMGCCGCGGTAA'
rev_806R <- 'GGACTACNVGGGTWTCTAAT'

# only need a tiny match so lets reduce that sequence
fwd_515F_partial <- 'CGGTAA'
rev_806R_partial <- 'GGACTA'

# FUNCTIONS I USED ####

# function for all chromatograms
get_chromatogram <- function(seq_file, output_folder){
  temp <- sangerseqR::read.abif(seq_file)
  seq_sang <- sangerseqR::sangerseq(temp)
  sangeranalyseR::secondary.peaks(seq_sang, output.folder = output_folder, file.prefix = tools::file_path_sans_ext(basename(seq_file)))
}

# check whether primer is present in each sequence
primer_presence <- function(seq_file, fwd_primer_seq, rev_primer_seq){
  # create reverse complement sequence
  rev_primer_seq <- chartr('ATGC', 'TACG', rev_primer_seq)
  
  # load in .ab1 file and get sequence
  temp <- sangerseqR::read.abif(seq_file)
  seq_sang <- sangerseqR::sangerseq(temp)
  seq <- temp@data$PBAS.2
  
  # create dataframe
  d <- data.frame(name = tools::file_path_sans_ext(basename(seq_file)), 
                  fwd_primer_present = stringr::str_detect(fwd_primer_seq, seq),
                  rev_primer_present = stringr::str_detect(rev_primer_seq, seq), 
                  stringsAsFactors = FALSE)
  return(d)
  
}

# trim and resave trimmed sequence
trim_and_save <- function(seq_file, output_path, trim_cutoff = 1e-04){
  temp <- sangerseqR::read.abif(seq_file)
  trims <- sangeranalyseR::trim.mott(temp, cutoff = trim_cutoff)
  seq <- substring(temp@data$PBAS.2, trims$start, trims$finish)
  write(seq, paste(output_path, '/', tools::file_path_sans_ext(basename(seq_file)), '.txt', sep = ''))
}

# load in trimmed files and make into a StringSet from which we can find consensus sequences
read_and_bind <- function(trimmed_file){
  temp <- data.frame(file = basename(trimmed_file), seq = read.table(trimmed_file, stringsAsFactors = FALSE)$V1, stringsAsFactors = FALSE)
  return(temp)
}

# output folder figs ####
fig_path <- 'sanger/figs/chromatogram'

# output data folder
trimmed_path <- 'sanger/trimmed_sanger'

# list files in the data folder - .ab1 files
files <- list.files('sanger/raw', pattern = '.ab1', full.names = TRUE)

# read in a single sequence
seq_abif = read.abif(files[3])
seq_sang <- sangerseq(seq_abif)

# search for primer sequences in my returned sequences
d_prim_pres <- map_df(files, primer_presence, fwd_primer_seq = fwd_515F_partial, rev_primer_seq = rev_806R_partial)
# they are not present in any of the files
# if they were there we would start and end our sequences at the position of the primer

# get all chromatograms
map(files, get_chromatogram, output_folder = fig_path)

# trim each file, save sequence out as the raw sequence ####

# get summary data for each file
file_sum <- summarise.abi.folder('sanger/raw')
View(file_sum$summaries)
write.csv(file_sum$summaries, 'sanger/20171009_sanger_seq_qualcheck.csv')
# the trimming parameters will improve the quality of the files

# do trimming and re-save files
walk(files, trim_and_save, output_path = trimmed_path, trim_cutoff = 1e-04)

trimmed_files <- list.files(output_path, full.names = TRUE)

# bind all files together
d_trim <- map_df(trimmed_files, read_and_bind) %>%
  mutate(., seq_len = nchar(seq))

# filter out files that have < 100 base pairs
d_trim <- filter(d_trim, seq_len >= 100)
# 4 have dropped out

# make a string set
d_SS <- DNAStringSet(d_trim$seq)
names(d_SS) <- d_trim$file

# try and find consensus sequences
alignment <- msa(d_SS)

msaPrettyPrint(alignment, output="pdf",
               showNames="none", showLogo="none", askForOverwrite=FALSE)

# get consensus sequence
printSplitString <- function(x, width=getOption("width") - 1){
  starts <- seq(from=1, to=nchar(x), by=width)
  
  for (i in 1:length(starts))
    cat(substr(x, starts[i], starts[i] + width - 1), "nn")
}

printSplitString(msaConsensusSequence(alignment))

# This may be the limit of the data if I cannot get multiple sequence alignments...



