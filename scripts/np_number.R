data <- read.delim("HPvec_C1C2_table.txt", stringsAsFactors=F)

pgroup2np <-  function(group) {
  meta <- unlist(strsplit(group, " ||||| ", fixed=T))
  meta <- meta[grep("REFSEQ", meta, fixed=T)]
  q <- function (pg) { gsub(".*REFSEQ:","", pg, perl=T) }
  q2 <- function (pg) { gsub("\\|.*", "", pg, perl=T) }
  q3 <- function (pg) { gsub("\\s.*", "", pg, perl=T) }
  meta <- unlist(lapply(meta,q))
  meta <- unlist(lapply(meta,q2))
  meta <- unlist(lapply(meta,q3))
  if(length(meta) > 0){
    c <- paste(meta, collapse=" ")
    gsub(";", " ", c, fixed=T)
  }
  else {
    "-"
  }
}

data$np.numbers <- unlist(lapply(data$protein.group, pgroup2np))

write.table(data, "HPvec_C1C2_table_withnp.txt", quote=FALSE, sep='\t', row.names=FALSE)




