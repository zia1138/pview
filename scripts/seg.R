#data <- read.delim("Phospho1_BY-RM-Orbi_peptides.txt", stringsAsFactors=F)
data <- read.delim("QTOF_peptides.txt", stringsAsFactors=F)

# Remove columns that represent median, they're redundant.
data <- data[,-grep("median", names(data))]
#data <- data[ -grep("reverse_", data$protein.group, fixed=T),]


# Get schizo columns from table.
BY <- data[,grep("BY", names(data), perl=T)]
# Get the paired control samples.
RM <- data[,grep("RM", names(data), perl=T)]


pvalues <- c()
for( i in 1:dim(data)[1]) {
  by <- log2(as.numeric(BY[i,]))
  rm <- log2(as.numeric(RM[i,]))
  by <- by[is.na(by) == F]
  rm <- rm[is.na(rm) == F]

  test <- t.test(by,rm)
  pvalues <- c(pvalues, test$p.value)
}

library(qvalue)
g <- qvalue(pvalues, fdr.level=0.05)

pdf("boxplots.pdf")

for( i in 1:dim(data)[1]) {
  by <- log2(as.numeric(BY[i,]))
  rm <- log2(as.numeric(RM[i,]))
  by <- by[is.na(by) == F]
  rm <- rm[is.na(rm) == F]
  if(g$significant[i]) {
    boxplot(by, rm, paired=T, names=c("by", "rm"), mex=0.1)

    id <- data[i,]$group.id
    name <- data[i,]$protein.group
    seq <- data[i,]$ms2.seq
    mz <- data[i,]$group.mz
    rt <- data[i,]$group.retentionTime

    qvalue <- g$qvalue[i]
    ratio <- 2^(mean(rm) - mean(by))
    title(paste(round(ratio,3), " ", substr(name, 1,15),"\nseq=",seq,"\nid=", id,"mz=",round(mz,5),"rt=", round(rt,0), "qvalue=",round(qvalue,3), sep=""))
  }
}

dev.off()
