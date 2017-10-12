# analysis of sanger sequencing data with dada2

# read in 1 abif file

# load packages
library(dplyr)
library(dada2)
library(ChIPsim)

# list files in the data folder - .fastq files
files <- list.files('sanger/fastq/raw', pattern = '.fastq', full.names = TRUE)
files2 <- files[files == 'sanger/fastq/raw/96sampl.fastq']
files <- files[files != 'sanger/fastq/raw/96sampl.fastq']

temp <- sangerseqR::read.abif('sanger/raw/001_js_a16_515F.ab1')
temp@data$PCON.2

# specify filter path
filt_path <- 'sanger/fastq/trimmed'

# specify names of filtered files
filt_names <- file.path(filt_path, basename(files))
filt_names2 <- file.path(filt_path, basename(files2))

temp <- '~/Desktop/test2.fastq'

# plot quality profiles
pdf('sanger/figs/plot_qual_temp.pdf')
for(i in temp){
  print(plotQualityProfile(i))
}
dev.off()
pdf('sanger/figs/plot_qual_all.pdf')
for(i in files2){
  print(plotQualityProfile(i))
}
dev.off()

# filter and trim
out <- filterAndTrim(files2, filt_names2,
                     maxN = 0, 
                     rm.phix=TRUE,
                     compress=TRUE,
                     trimLeft = 40,
                     minLen = 50,
                     multithread=TRUE)

# these are pretty stringent criteria - cannot seem to get them to pass well
out

# dereplicate samples
filt_names <- list.files('sanger/fastq/trimmed', pattern = '.fastq', full.names = TRUE)
filt_names2 <- filt_names[filt_names == 'sanger/fastq/trimmed/96sampl.fastq']

errors <- learnErrors(filt_names2, multithread=TRUE)
plotErrors(errors)

derep_files <- derepFastq(filt_names2)

# sample inference
dadaFs <- dada(derep_files, err = NULL, selfConsist = TRUE)

# assign taxonomy
rdp_taxa <- assignTaxonomy(dadaFs, "sanger/ref_trainsets/rdp_train_set_16.fa", multithread = TRUE)
rdp_taxa <- addSpecies(rdp_taxa, "sanger/ref_trainsets/rdp_species_assignment_16.fa", allowMultiple = 2)
rdp_taxa <- compile_taxonomy(rdp_taxa, tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

temp <- sangerseqR::read.abif(files[1])

# extract sequence
seq <- temp@data$PBAS.2

# extract quality score
qual_scores <- temp@data$PCON.2

# make string set
seq2 <- BStringSet(seq)
names(seq2) <- tools::file_path_sans_ext(basename(files[1])) 

fq <- signature(sread = seq, quality = qual_scores, id = '')


fq <- signature(sread = c('ATGC'), quality = c('1234'), id = "missing")

# read in a single file
# check if the error profile changes
# change the error class
# check on one of the Illumina files from my analysis (Experiment 4)
# see here 
x <- readFastq(files[2])
qual <- as.character(quality(quality(x)))
errorProb <- decodeQuality(qual, type="Sanger")
qualSanger <- encodeQuality(errorProb, type="Sanger")
qualIllumina <- encodeQuality(errorProb, type = 'Illumina')
all.equal(qual, qualIllumina)

set <- ShortReadQ(sread(x), quality = BStringSet(qualIllumina), id = BStringSet(files[1]))
writeFastq(set, '~/Desktop/test2.fastq')

# get numeric quality
as(quality(x), "numeric")
as(quality(set), "numeric") # NOTHING HAS CHANGED HERE

x <- readFastq('~/Desktop/test2.fastq')
as.character(quality(quality(x))) == qualIllumina
