file.from.zia='QTOF_peptides.txt'
file.cor.plots='phospho_parents.pdf'

pho=read.delim(file.from.zia, sep='\t', header=TRUE, stringsAsFactors=FALSE)
library(qvalue)

# Remove columns that represent median, they're redundant.
pho <- pho[,-grep("median", names(pho))]
#pho <- pho[ -grep("reverse_", pho$protein.group, fixed=T),]

# grab data
data.only=pho[,grep("rep.Rep1", names(pho), fixed=T)]
data.only=log2(data.only)
# clean up labels
names(data.only)=gsub('cond.', '', names(data.only))
#names(data.only)=gsub('rep.Rep1.scan.M10.', '', names(data.only))
names(data.only)=gsub('rep.Rep1.scan.E10_', '', names(data.only))

# input functions for pairs(see line 29)
panel.cor=function(x,y){
    test = cor.test(x,y, method='spearman', use='pairwise.complete.obs')
    #text(21,26, paste('Rho =', round(test$estimate,3)),cex=1.8, col='red')
    text(18,18, paste(round(test$estimate,3)),cex=1.8, col='red')
    # text(21,22, '95% conf int ->', cex=.8)
    # conf.left=round(test$conf.int,3)[1]
    # conf.right=round(test$conf.int,3)[2]
    # text(21,20, paste(conf.left, ' - ', conf.right, sep=''), cex=1.2, col='red')
}
plot.pairs = function(x,y){
    points(x,y, pch=20, cex=.4)
    # abline(0,1)
}

pdf(file=file.cor.plots, width=11, height=11)
pairs(data.only, upper.panel=plot.pairs, lower.panel=panel.cor)
dev.off()


# sig differences
ttest.results=apply(data.only, 1, function(x) {
            by=x[1:3]
            rm=x[4:6]
            by=na.omit(by)
            rm=na.omit(rm)
            if(length(by)>1 & length(rm)>1) {   return(t.test(by,rm)$p.value)  }
            else{return(NA)}
})

# multiple testing correction
qvalue.results=qvalue(ttest.results)
qs=qvalue.results$qvalue

#pdf(file='~/pvaluedist.pdf', width=11, height=11)
#hist(ttest.results, breaks=500, main='p-value distribution', xlab='p')
#dev.off()

#output sig pep count
length(qs[qs<.05])

#output for ludovic hand confirmation
#report.ind=which(qs<.05)
#prots=pho[report.ind,]
#out=data.frame(prots, ttest.results, qs)
#write.table(prots, file='~/Desktop/sig_prot_full.txt', sep='\t', quote=FALSE)
