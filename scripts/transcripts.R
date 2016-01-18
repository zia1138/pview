                                        # Load microRNA gene expression data
data <- read.delim("all-data.csv", stringsAsFactors=F)

                                        # Compute average of two cell types compared.
log2.HP.avg <- log2((data$HP + data$HP.2) / 2)
log2.C1C2.avg <- log2((data$C1C2 + data$C1C2.2) / 2)

                                        # Compute log2 of the ratio
data$log2.ratio.g <- log2.C1C2.avg - log2.HP.avg

                                        # Remove the single quote in the gene name.
data$Gene <- gsub("'", "", data$Gene, fixed=T)
data$Gene <- toupper(data$Gene)

                                        # Only keep the log2.ratio and the gene name.
g.data2 <- data.frame(gene = data$Gene, log2.ratio.g = data$log2.ratio.g, stringsAsFactors=F)
g.data2 <- g.data2[!duplicated(g.data2$gene),]
                                        # Free up the data.
rm(data)

                                        # Load protein abundance data.
p.data <- read.delim("Manav2009.csv", stringsAsFactors=F)

                                        # Save only protein group information and ratio.
p.data2 <- data.frame(gene = p.data$protein.group, log2.ratio.p = p.data$group.log2.HL.ratio, stringsAsFactors=F)
p.data2 <- unique(p.data2)


                                        # Map protein group to gene symbol.
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

                                        # Convert to gene symbol, remove duplicates.
p.data2$gene <- unlist(lapply(unique(p.data2$gene), pgroup2gene))
p.data2$gene <- toupper(p.data2$gene)
p.data2 <- p.data2[p.data2$gene != "-",]
p.data2 <- p.data2[!duplicated(p.data2$gene),]

cor <- merge(p.data2, g.data2)

print(dim(cor))

miR.200a <- read.delim("miR200a.csv", stringsAsFactors=F)
miR.200a$Gene <- toupper(miR.200a$Gene)

miR.200b <- read.delim("miR200b.csv", stringsAsFactors=F)
miR.200b$Gene <- toupper(miR.200b$Gene)

targets <- data.frame(gene = union(miR.200a$Gene, miR.200b$Gene))

cor.targets <- merge(targets, cor)
cor.targets$gene <- as.character(cor.targets$gene)

# Generate ecorrelation plot between protein abundance and mRNA abundance.
pdf("Manav2009/transcripts.pdf")

plot(c(-3.5,3.5), c(-3.5,3.5), type = "n", xlab="log2 mRNA ratio", ylab="log2 protein ratio" )
abline(0,1, col="gray", lty=1, lwd=2)
points(cor$log2.ratio.g, cor$log2.ratio.p, xlim=c(-4,4), ylim=c(-4,4), pch=20, cex=.3)

c <- cor.test(cor$log2.ratio.g, cor$log2.ratio.p, method="spearman")
title(paste("Spearman's rho=", round(c$estimate,2), ", N=", dim(cor)[1], sep=""))

dev.off()


# Genes that show a down regulation of 40% in both transcripts and proteins.
cand <- cor.targets[cor.targets$log2.ratio.g < log2(0.6) & cor.targets$log2.ratio.p < log2(0.6),]

# Plot only transcripts and proteins.
pdf("Manav2009/transcripts2.pdf")
plot(c(-3.5,3.5), c(-3.5,3.5), type = "n", xlab="log2 mRNA ratio", ylab="log2 protein ratio" )
abline(0,1, col="gray", lty=1, lwd=2)
points(cor.targets$log2.ratio.g, cor.targets$log2.ratio.p, xlim=c(-4,4), ylim=c(-4,4), pch=20, cex=0.3)

c <- cor.test(cor.targets$log2.ratio.g, cor.targets$log2.ratio.p, method="spearman")

# Show names of genes with 40% decrease in transcripts and proteins.
for (i in 1:dim(cand)[1]) {
  text(cand[i,]$log2.ratio.g, cand[i,]$log2.ratio.p, cand[i,]$gene, pos=4, offset=0.2, cex=0.55, col="red")
}

title(paste("Spearman's rho=", round(c$estimate,2), ", N=", dim(cor.targets)[1], sep=""))
dev.off()

cand2 <- cor[cor$log2.ratio.g < log2(0.6) & cor$log2.ratio.p < log2(0.6),]

print(paste("enrichment p-value", phyper(lower.tail = TRUE, dim(cand)[1], dim(targets)[1], dim(g.data2)[1], dim(cand2)[1])))


targets$has.target <- TRUE
cor.targets2 <- merge(cor, targets, all.x=TRUE)
cor.targets2$has.target[is.na(cor.targets2$has.target)] <- FALSE


cor.targets2 <- rbind(cor.targets2[cor.targets2$log2.ratio.g < log2(0.6) & cor.targets2$log2.ratio.p < log2(0.6),],
                      cor.targets2[!(cor.targets2$log2.ratio.g < log2(0.6) & cor.targets2$log2.ratio.p < log2(0.6)),])
      
write.table(cor.targets2, "TableS2.csv", sep="\t", row.names=FALSE)



pdf("Manav2009/transcripts3.pdf")
plot(c(-3.5,3.5), c(-3.5,3.5), type = "n", xlab="log2 mRNA ratio", ylab="log2 protein ratio" )
abline(0,1, col="gray", lty=1, lwd=2)
points(cor$log2.ratio.g, cor$log2.ratio.p, xlim=c(-4,4), ylim=c(-4,4), pch=20, cex=.3)
points(cor.targets$log2.ratio.g, cor.targets$log2.ratio.p, xlim=c(-4,4), ylim=c(-4,4), pch=20, cex=0.8, col="red")
dev.off()
