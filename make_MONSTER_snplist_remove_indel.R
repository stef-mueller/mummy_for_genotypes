#!/usr/bin/env Rscript


#
# Objective: Identify Indels (small insertions/deletions) in variant list and remove them
#

suppressMessages(require(tidyverse, quietly=T, warn.conflicts = F))
require(optparse, quietly = T, warn.conflicts = F)

### user input
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="variant IDs extracted from MONSTER genotype file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read data in from variant column of MONSTER genotype file
varsIN = suppressMessages(read_delim(opt$input,
 delim= "_", col_names = c("CHR", "POS", "REF", "ALT"),
  col_types ="iicc"))

# function to test alleles are only coded as ACTG and not for example "-" to indicate INDEL
isNucleotide = function(x){
  !grepl("[^ACTG]+$", x)
}

# remove Indels by counting character length
varsKEEP = varsIN %>% 
  filter(isNucleotide(REF) & isNucleotide(ALT)) %>% 
  filter(nchar(REF)<2) %>% 
  filter(nchar(ALT)<2)

# save output
varsKEEP %>% 
  mutate(OUT = paste(CHR,POS,REF,ALT, sep="_")) %>% 
  mutate(TABIX = paste0(CHR,":",POS,"-",POS)) %>% 
  select(OUT,TABIX, ALT, REF) %>% 
  write.table("monster_snps.txt", col.names = F, row.names = F, quote = F, sep ="\t")
