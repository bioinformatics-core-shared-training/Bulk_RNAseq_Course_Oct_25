# load libraries
library(tximport)
library(DESeq2)
library(tidyverse)


# load sample metadata
sampleinfo <- read_tsv( 'data/samplesheet.tsv' )

# read counts data

files <- file.path( "salmon", sampleinfo$SampleName, 'quant.sf' )
files <- set_names(files, sampleinfo$SampleName)
files

tx2gene <- read_tsv('references/tx2gene.tsv')

# load data
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)

# explore txi obeject

class(txi)
names(txi)
txi$abundance %>%  head()

txi$counts %>% head()
txi$length %>% head()

# save txi object
saveRDS( txi, file = "salmon_outputs/txi.rds")


# explore the counts data
rawCounts <- round(txi$counts, digits = 0)

head(rawCounts)

# filter raw counts
dim(rawCounts)

keep <- rowSums(rawCounts) > 5 

table(keep)

# create filtered counts matrix
filtCounts <- rawCounts[ keep, ]

dim(filtCounts)

# raw counts box-plot

boxplot(filtCounts, main="Raw counts", las=2 )

# raw counts: mean-variance relationship

plot( rowMeans(filtCounts), 
      rowSds(filtCounts),
      main = 'Raw counts: mean vs sd',
      xlim=c( 0, 10000),
      ylim = c( 0, 5000)
      )


filtCounts %>% 
  as.data.frame() %>% 
  View()

# log2 data transformation
logCounts <- log2( filtCounts + 1)


boxplot(logCounts, main = 'Log2 counts', las=2)

# log counts: mean-variance relationship

plot( rowMeans(logCounts), 
      rowSds(logCounts),
      main = 'Log counts: mean vs sd'
)



# rlog transformation
rlogcounts <- rlog( filtCounts )

head(rlogcounts)

boxplot(rlogcounts, main='rlog counts', las=2)


# log counts: mean-variance relationship

plot( rowMeans(rlogcounts), 
      rowSds(rlogcounts),
      main = 'rlog counts: mean vs sd'
)

# PCA analysis
library(ggfortify)

rlogcounts <- rlog(filtCounts)

# create PCA
pcDat <- prcomp( t(rlogcounts) )

# plot PCA
autoplot(pcDat,
         data=sampleinfo,
         colour = 'Status',
         shape = 'TimePoint',
         size= 5)


