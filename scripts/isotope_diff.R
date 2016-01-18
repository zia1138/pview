data <- read.delim("HPvec_C1C2.txt", stringsAsFactors=F)


#data <- data[data$pair.c == "0h",]

# Get uniq set of protein group IDs.
groups <- unique(data$group.id)

# For each protein group get the log2 ratios.
p.values <- c()
group.id <- c()
for(g in groups) {
  log2ratios <- data[data$group.id == g,]$log2.HL.ratio

  if(length(log2ratios) > 1) {
    # Determine how significantly different these values are from zero.
    test <- t.test(log2ratios)
    p.values <- c(p.values, test$p.value)
    group.id <- c(group.id, g) # Save p-value and group ID.
  }
}
# Apply q-value FDR correctioin
library(qvalue)
d <- qvalue(p.values, fdr=0.05)

# Collect significant subset.
sig <- group.id[d$significant]

p <-  rep(FALSE, length(data$group.id))
for(s in sig) { p <- p | data$group.id == s; }

significant <- data[p,]

#down <- unique(significant[significant$group.log2.HL.ratio < 0,]$protein.group.avg.mass)
#up <- unique(significant[significant$group.log2.HL.ratio > 0,]$protein.group.avg.mass)

#pdf("PsuperSh61_updown.pdf")
#boxplot(down, up, ylab="mass in kDa", names=c("down","up"))
#dev.off()

#write.table(significant, "PsuperSh61_sig.txt", quote=FALSE, sep='\t', row.names=FALSE)

write.table(significant, "HPvec_C1C2_sig.txt", quote=FALSE, sep='\t', row.names=FALSE)

