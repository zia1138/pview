# R script contains different reporting plot generating functions.

group.plot <- function(dir, descr, subset, max.log2, m, dxlab = "log2 ratio group 1", dylab="log2 ratio group 2") {
  pdf(paste(dir, "/internal_", descr, ".pdf", sep=""))
  plot(c(-max.log2,max.log2), c(-max.log2,max.log2), type = "n", xlab=dxlab, ylab=dylab)
  abline(0,1, col="gray", lty=1, lwd=2)
  points(subset$group1.log2.ratio, subset$group2.log2.ratio, pch=20, cex=0.45)
  if(length(subset$group1.log2.ratio) > 1) {
    # only compute correlation if there's enough data!
    c <- cor.test(subset$group1.log2.ratio, subset$group2.log2.ratio, method="spearman")
    title(paste("Spearman's rho=", round(c$estimate,2), ", N=", dim(subset)[1], ", m>", m, sep=""))
  }
  dev.off()
}

mass.error.report <- function(dir, txt.file, outfile = "mass_error.pdf") {
                                        # Read isotope pairs.
  data <- read.delim(txt.file, stringsAsFactors=F)
                                        # Generate histgram of precursor mass error.
                                        # Computed mean and standard deviation of mass error.
  mu <- mean(data$ms2.ppm.error)
  stdev <- sd(data$ms2.ppm.error)
                                        # Generate the histogram.
  pdf(paste(dir, "/", outfile, sep=""))
  hist(data$ms2.ppm.error, breaks="Freedman-Diaconis", xlab="measured - calculated mass (p.p.m.)",
       main=paste("mean=",round(mu,3)," s.d.=", round(stdev,3),sep=""))
  dev.off()
}

# For a data set the save isotope pair should generate 3 TXT files.
# DataSetName.txt
# DataSetName_corr.txt - across data set correlation data
# DataSetName_mass_error.txt - mass error file
# DataSetName_internal - internal correlation data file
# Each of these three files is used to generate plots that report algorithm performance.
isotope.pair.report <- function(pair.txt.file, max.log2 = 2, m.internal = 1, m.across = 0) { 

  internal.corr.report <- function(dir, txt.file, group.N.cutoff, max.log2, check.mods = 0) {
    data <- read.delim(txt.file, stringsAsFactors=F)
                                        # For each unique description, generate an internal correlation.
                                        # only report internal correlation for modified peptides
    if(check.mods > 0) {
      if(check.mods == 1) { data <- data[data$Nmods == 0,] }
      else if(check.mods == 2) { data <- data[data$Nmods > 0,] }
    }
    
    descriptions <- unique(data$description)
    for(d in descriptions) {
      subset <- data[data$description == d, ]
                                        # Enforce minimum group size.
      subset <- subset[subset$group1.N + subset$group2.N > group.N.cutoff,]
                                        # Generate PDF plot.
      if(check.mods > 0) {
        if(check.mods == 2) { d <- paste(d, "_mod", sep="")  }
        group.plot(dir, d, subset, max.log2, group.N.cutoff, "log2 ratio charge = +2", "log2 ratio charge > +2")
      }
      else { group.plot(dir, d, subset, max.log2, group.N.cutoff) }
    }
  }
                                        # Generate plots correlating different conditions.
  cor.report <- function (dir, cor.file, group.N.cutoff, max.log2) {
    data <- read.delim(cor.file, stringsAsFactors=F)
                                        # Get all combined correlation IDs
    corr.ids <- unique(data$corr.id)

    for(c in corr.ids) {
      to.corr <- data[data$corr.id == c, ]
                                        # Read off replicate descriptions
      rep1 <- to.corr$rep1.id[1]
      rep2 <- to.corr$rep2.id[1]
                                        # Apply cutoff on replicates
      to.corr <- to.corr[to.corr$ratio1.N > group.N.cutoff & to.corr$ratio2.N > group.N.cutoff, ]

                                        # Generate replicate correlation plots.
      pdf(paste(dir, "/corr_", c, ".pdf", sep=""))
      plot(c(-max.log2,max.log2), c(-max.log2,max.log2), type = "n", xlab=paste(rep1, "log2 ratio"),  ylab=paste(rep2,"log2 ratio"))
      abline(0,1, col="gray", lty=1, lwd=2)
      points(to.corr$log2.ratio1, to.corr$log2.ratio2, pch=20, cex=0.45)
      c <- cor.test(to.corr$log2.ratio1, to.corr$log2.ratio2, method="spearman")
      title(paste("Spearman's rho=", round(c$estimate,2), " N=", dim(to.corr)[1], ", m>", group.N.cutoff, sep=""))
      dev.off()
    }
  }
                                        # Extract prefix.dir from pair TXT file.

  prefix.dir <- gsub(".txt", "", pair.txt.file, fixed=T)
                                        # Create directory for output
  if(file.exists(prefix.dir) == F) { dir.create(prefix.dir) }

                                        # Additional data files for output.
  internal.file <- paste(prefix.dir, "_internal.txt", sep="")
  internal.pep.file <- paste(prefix.dir, "_internal_pep.txt", sep="")
  cor.file <- paste(prefix.dir, "_corr.txt", sep="")
  mass.error.file <- paste(prefix.dir, "_mass_error_ms1.txt", sep="")
  mass.error2.file <- paste(prefix.dir, "_mass_error_ms2.txt", sep="")
  
                                        # If the corresponding data files exist, generate plots.
  if(file.exists(internal.file)) { internal.corr.report(prefix.dir, internal.file, m.internal, max.log2) }
  if(file.exists(internal.pep.file)) {
    internal.corr.report(prefix.dir, internal.pep.file, m.internal, max.log2, check.mods=1)
    internal.corr.report(prefix.dir, internal.pep.file, m.internal, max.log2, check.mods=2)
  }  
  if(file.exists(cor.file)) { cor.report(prefix.dir, cor.file, m.across, max.log2) }
  if(file.exists(mass.error.file)) { mass.error.report(prefix.dir, mass.error.file) }
  if(file.exists(mass.error2.file)) { mass.error.report(prefix.dir, mass.error2.file, "mass_error2.pdf") }

}

internal.align.report <- function(txt.file, group.N.cutoff = 1, max.log2=1.5) {
                                        # Extract dir from pair TXT file.
  dir <- gsub(".txt", "", txt.file, fixed=T)
                                          # Create directory for output
  if(file.exists(dir) == F) { dir.create(dir) }

                                        # Additional data files for output.
  internal.file <- paste(dir, "_internal.txt", sep="")

  data <- read.delim(internal.file, stringsAsFactors=F)
  descriptions <- unique(data$description)
  for(d in descriptions) {
    subset <- data[data$description == d, ]
                                        # Enforce minimum group size.
    subset <- subset[subset$group1.N + subset$group2.N > group.N.cutoff,]
                                        # Generate PDF plot.
    group.plot(dir, d, subset, max.log2, group.N.cutoff)
  }

  mass.error.file <- paste(dir, "_mass_error.txt", sep="")
  mass.error2.file <- paste(dir, "_mass_error_ms2.txt", sep="")
  
  mass.error.report(dir, mass.error.file)
  mass.error.report(dir, mass.error2.file, "mass_error2.pdf")
}


# Compute paired t-test
diff.express.cand <- function (pair.txt.file, cand.txt.file = "candidates.txt", q.value.cutoff = 0.01, log2.fold.change = log2(1.5)) {
                                        # Extract prefix.dir from pair TXT file.
  prefix.dir <- gsub(".txt", "", pair.txt.file, fixed=T)

  require(qvalue)
                                        # Read TXT data file.
  pairs <- read.delim(pair.txt.file, stringsAsFactors=F)

  pairs$group.log2.HL.ratio <- pairs$group.log2.HL.ratio - median(pairs$log2.HL.ratio)
  #pairs$log2.HL.ratio += median(pairs$log2.HL.ratio)
                                        # Elminate singleton groups as they are not valid for paired t-test.
  pairs <- pairs[pairs$group.N >= 2, ]

                                        # Compute paired t.test for log2 of XIC intensities
  cand <- NULL
  for(group.id in unique(pairs$group.id)) {
    # Get isotope pairs for a protein group
    group <- pairs[pairs$group.id == group.id,]

    # Compute paired t-test on the log2 intensities of each isotope pair.
    res <- t.test(group$log2.xicH, group$log2.xicL, paired=TRUE)

    p.value <- res$p.value
    group.log2.HL.ratio <- group$group.log2.HL.ratio[1]
    group.N <- group$group.N[1]
    protein.group <- group$protein.group[1]
    # Save as candidtes for FDR correction.
    cand <- rbind(cand, data.frame(group.id, protein.group, group.log2.HL.ratio, group.N, p.value))
  }

  # Compute FDR corrected q-value.
  cand$q.value <- qvalue(cand$p.value)$qvalues

  # Save data to a table.
  cand <- cand[cand$q.value < q.value.cutoff & abs(cand$group.log2.HL.ratio) > log2.fold.change ,]
  write.table(cand, paste(prefix.dir,cand.txt.file, sep="/"), quote=FALSE, sep='\t', row.names=FALSE)
}
