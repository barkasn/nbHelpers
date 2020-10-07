
counts_to_gct <- function(counts = NULL, file = NULL, species.org.eg.db = NULL) {
  # Calculate cpm from counts
  cpm <- sweep(x=counts * 1e6,MARGIN=2,STATS=apply(counts,2,sum),FUN='/') 
  
  # Map the ensembl identifiers to gene labels and create dataframe
  ens_id_noversion <- strpart(rownames(cpm),'.',1,fixed=T)
  ens_symbol_map <- select(species.org.eg.db, keys=ens_id_noversion, columns=c('ENSEMBL','SYMBOL'), keytype='ENSEMBL')
  rownames(cpm) <- ens_symbol_map$SYMBOL[match(ens_id_noversion, ens_symbol_map$ENSEMBL)]
  cpm <- cpm[!is.na(rownames(cpm)),]
  cpm.df <- data.frame(cpm)
  
  # Generate final dataframe
  
  df1 <- cbind(data.frame(NAME=rownames(cpm),Description=c("NA")), cpm.df)
  rownames(df1) <- NULL
  
  # remove duplicates
  df1 <- df1[!duplicated(df1$NAME),]
  
  output_dim <- dim(df1)
  output_dim[2] <- output_dim[2] - 2
  
  # write output file
  fh <- file(description = file, open="w" )
  writeLines(c("#1.2",paste0(output_dim[1],"\t",output_dim[2])),fh)
  write.table(df1,file = fh, sep='\t',quote=FALSE,row.names=FALSE,col.names = TRUE)
  close(fh)
}


writeClsFile <- function(x, file) {
  fh <- file(file, "w")
  writeLines(
    c(
      paste(length(x), nlevels(x), 1),
      paste(c('#',levels(x)), collapse=' '),
      paste(as.numeric(as.factor(as.character(x)))-1, collapse=' ')
    ),
    con = fh
  )
  close(fh)
}
