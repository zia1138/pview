data <- read.delim("DeplCys.txt", stringsAsFactors=F)

# Remove columns that represent median, they're redundant.
data <- data[,-grep("median", names(data))]

# Get schizo columns from table.
schizo <- data[,grep("115\\.01|117\\.01|119\\.02|120\\.01|122\\.01|124\\.02|127\\.02|203\\.01|204\\.01|398\\.01", names(data), perl=T)]
# Get the paired control samples.
control <- data[,grep("115\\.02|117\\.02|119\\.01|120\\.02|122\\.02|124\\.01|127\\.10|203\\.02|204\\.01|398\\.03", names(data), perl=T)]

pdf("boxplots.pdf")

pvalues <- c()
for( i in 1:dim(data)[1]) {
  s <- log2(as.numeric(schizo[i,]))
  c <- log2(as.numeric(control[i,]))
  s <- s[is.na(s) == F]
  c <- c[is.na(c) == F]
  if(length(s) > 1 & length(c) > 1) {
    test <- t.test(s,c)
    pvalues <- c(pvalues, test$p.value)
    if(test$p.value < 0.05) {
      boxplot(s, c, paired=T, names=c("schizo", "control"))

      id <- data[i,]$group.id
      name <- data[i,]$protein.group
      mz <- data[i,]$group.mz
      rt <- data[i,]$group.retentionTime
      
      title(paste(substr(name,1,20)," id= ", id,"mz = ",round(mz,5),"rt =", round(rt,0), "pvalue = ",round(test$p.value,3)))
    }
  }
}
dev.off()

library(qvalue)
g <- qvalue(pvalues, fdr.level=0.1)
