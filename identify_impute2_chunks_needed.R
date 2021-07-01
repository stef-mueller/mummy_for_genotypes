#!/usr/bin/env Rscript

suppressMessages(require(tidyverse, quietly=T, warn.conflicts = F))
require(data.table, quietly=T, warn.conflicts = F)
require(optparse, quietly = T, warn.conflicts = F)

### user input
option_list = list(
  make_option(c("-i", "--index"), type="character", default=NULL, 
              help="index file for impute2 chunks name", metavar="character"),
    make_option(c("-b", "--bed"), type="character", default=NULL,
              help="bed file with regions of interest", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load info to chromosomal positions of impute2 chunks
index = read.delim(opt$index)
# load info bed file
bed = read.table(opt$bed, header=F)
names(bed) = c("chr_bed", "start_bed", "end_bed", "name_bed")

chrom = index %>% 
  filter(chr == bed$chr_bed[1])

out = list()
for (i in 1:nrow(bed)){

  current = bed %>% 
    slice(i)
  
  # find standard overlap: region included in one impute2 chunk
  tmp1 = chrom %>% 
        filter(start <= current$start_bed & end >= current$end_bed)  
  
  # if succesful include in output
  if (nrow(tmp1)==1){
    
    out[[i]] = tmp1
    
  } else { # otherwise test for special case of region overlapping with split between chunks
    
    # first case: chunk which overlaps at start by calculating distance between region start and chunk start
    tmp2 = chrom %>% mutate(diffstart = start-current$start) %>% 
      filter(diffstart < 0) %>% 
      arrange(desc(diffstart)) %>% 
      slice(1) %>% 
      select(-diffstart)
    
    # second case: chunk which overlaps at end by calculating distance between region end and chunk end
    tmp3 = chrom %>% mutate(diffend = end-current$end) %>% 
      filter(diffend > 0) %>% 
      arrange(diffend) %>% 
      slice(1) %>% 
      select(-diffend)
    
    out[[i]] = rbind(tmp2, tmp3)
  }
  
  
}

outdf = data.table::rbindlist(out) %>% 
  distinct()

write.table(outdf, "impute2needed.txt", quote=F, col.names = F, row.names = F, sep ="\t")
