pgroup2gene <- function (group) {
  meta <- unlist(strsplit(group, " ||||| ", fixed=T))
  meta <- unlist(strsplit(gsub(".+Gene_Symbol=", "", meta), " ||||| ", fixed=T))

                                        # Get only the first gene symbol if there are many.
  first.symbol <- function(p)  { substr(p, 1, regexpr('\\s|;', p, perl=T)[1] - 1)  }
  
  gene.symbols <- unlist(lapply(meta, first.symbol))
  gene.symbols <- gene.symbols[gene.symbols != "" & gene.symbols != "-" ]

                                        # Remove some wierd gene symbols
  if(length(grep("LOC", gene.symbols)) != 0) { gene.symbols <- gene.symbols[-grep("LOC", gene.symbols)] }
  if(length(grep("OTT", gene.symbols)) != 0) { gene.symbols <- gene.symbols[-grep("OTT", gene.symbols)] }
  if(length(unique(gene.symbols)) == 1) { gene.symbols[1] }
  else {
                                        # Return nothing if no symbols or multiple symbols.
    "-"
  }
}



condition <- "HPvec_C1C2"
ratio.file <- "HPvec_C1C2.txt"
internal.corr.file <- "HPvec_C1C2_internal.txt"

ratio.data <- read.delim(ratio.file, stringsAsFactors=F)
ratio.data <- ratio.data[ratio.data$pair.c == condition,]

internal.data <- read.delim(internal.corr.file, stringsAsFactors=F)
internal.data <- internal.data[internal.data$description == condition,]

inliers <- abs(internal.data$group1.log2.ratio - internal.data$group2.log2.ratio) < log2(1.2)
changers <- abs(internal.data$group1.log2.ratio) >= log2(2) & abs(internal.data$group2.log2.ratio) >= log2(2)

inlier.data <- internal.data[inliers & changers == FALSE,]
outlier.data <- internal.data[inliers == FALSE, ]
cand.data <- internal.data[inliers & changers,]

pdf(paste(condition,"_cand.pdf", sep=""))

plot(internal.data$group1.log2.ratio, internal.data$group2.log2.ratio, xlim=c(-3,3),ylim=c(-3,3),
     type="n", xlab=("log2 ratio group 1"), ylab="log2 ratio group 2")
abline(0,1, col="gray")
points(outlier.data$group1.log2.ratio, outlier.data$group2.log2.ratio, pch=20,cex=0.3, col="gray")
points(inlier.data$group1.log2.ratio, inlier.data$group2.log2.ratio, pch=20,cex=0.3, col="red")

cand.data <-  merge(cand.data, ratio.data, by="group.id")
cand.data$gene <- unlist(lapply(cand.data$protein.group, pgroup2gene))

for(g in unique(cand.data$group.id)) {
  sub <- cand.data[cand.data$group.id == g,]
  x <- sub$group1.log2.ratio[1]
  y <- sub$group2.log2.ratio[1]
  text(x,y,sub$gene[1], adj=c(0,0), cex=0.6, offset=3)
  points(x,y, pch=20, cex=0.4, col="blue")
}


dev.off()

write.table(cand.data, paste(condition,"_cand.txt", sep=""), quote=FALSE, sep='\t', row.names=FALSE)
