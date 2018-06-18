library(ggplot2)
library(plotly)
library(reshape2)
library(phyloseq)

# Takes the relative abundances at an unranked level and
# upsamples them to higher taxonomical levels
upsample.ra <- function(ra.unranked, tax.table, newlev) {
  df = ra.unranked
  df$newlev = tax.table[[newlev]]
  df.melt = melt(df, id.vars = c("newlev"))
  df.melt.agg = aggregate(.~variable+newlev, data=df.melt, FUN=sum)
  df.ra = dcast(df.melt.agg, variable~newlev)
  rownames(df.ra) = df.ra$variable
  df.ra$variable = NULL 
  return(df.ra)
}

count.to.ra <- function(otu_count) {
  otu_relabu <- otu_count
  for (i in 1:ncol(otu_count))  {
    otu_relabu[,i] <- otu_relabu[,i]/sum(otu_relabu[,i])
  }
  return(otu_relabu)
}

# Trying to not use the pstat object directly
pstat.extraction <- function(p) {
  tables <- list()

  # Taxon Levels
  # unranked taxon x classifications
  tables$TAX <- as.data.frame(p@tax_table)

  # Relative Abundance
  # unranked taxon x samples
  tables$RLA <- as.data.frame(count.to.ra(p@otu_table))

  # Sample Data
  # samples x traits
  tables$SAM <- as.data.frame(p@sam_data)
  return(tables)
}

files.extraction <- function() {
  tables <- list()
  tables$TAX <- as.data.frame(read.table(file.path("data", "TAX_TABLE.txt"), header=TRUE, sep="\t", row.names = 1))
  OTU_TABLE = as.data.frame(read.table(file.path("data", "OTU_TABLE.txt"), header=TRUE, sep="\t", row.names = 1))
  tables$RLA <- count.to.ra(OTU_TABLE)
  tables$SAM <- as.data.frame(read.table(file.path("data", "SAM_TABLE.txt"), header=TRUE, sep="\t", row.names = 1))
  return(tables)
}