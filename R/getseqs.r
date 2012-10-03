getseqs <- function(taxon_name, gene, seqrange, getrelated)
{
#   taxon_name <- spvec[[1]]
  # Workflow
  message(paste("Working on ", taxon_name, "...", sep=""))
  ## Get GI numbers for a particular gene for a particular plant species
  message("...retrieving sequence IDs...")
  if(length(gene) > 1){ genes_ <- paste0(gene, collapse=" OR ") } else
    { genes_ <- paste0(gene, collapse=" ") }
  genes_ <- paste("(", genes_, ")")
  
  query <- list(db = "nuccore", term = paste(taxon_name, "[Organism] AND", genes_, "AND", seqrange, "[SLEN]", collapse=" "), RetMax=500)
  
#   query <- list(db = "nuccore", term = paste(taxon_name, "[Organism]", gene, "[Gene Name]", seqrange, "[SLEN]", collapse=" "), RetMax=500)
  
  out <- parsed_content(
    GET("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", query=query))$doc$children$eSearchResult
  if( as.numeric(xmlValue(xpathApply(out, "//Count")[[1]]))==0 ){
    message(paste("no sequences of ", gene, " for ", taxon_name, " - getting other sp.", sep=""))
    if(getrelated == FALSE){
      message(paste("no sequences of ", gene, " for ", taxon_name, sep=""))
      outt <- list(taxon_name, NA, NA, NA, NA, NA)
    } else
    {
      message("...retrieving sequence IDs for related species...")
      newname <- strsplit(taxon_name, " ")[[1]][[1]]
      query <- list(db = "nuccore", term = paste(newname, "[Organism] AND", genes_, "AND", seqrange, "[SLEN]", collapse=" "), RetMax=500)
#       query <- list(db = "nuccore", term = paste(newname, "[Organism]", gene, "[Gene Name]", seqrange, "[SLEN]", collapse=" "), RetMax=500)
      out <- parsed_content(
        GET("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", query=query))$doc$children$eSearchResult
      if( as.numeric(xmlValue(xpathApply(out, "//Count")[[1]]))==0 ){
        message(paste("no sequences of ", gene, " for ", taxon_name, " or ", newname, sep=""))
        outt <- list(taxon_name, NA, NA, NA, NA, NA)
      } else
      {
        ids <- xpathApply(out, "//IdList//Id")
        ids_ <- as.numeric(sapply(ids, xmlValue))
        
        ## For each species = get GI number with longest sequence
        message("...retrieving sequence ID with longest sequence length...")
        querysum <- list(db = "nucleotide", id = paste(ids_, collapse=" ")) # construct query for species
        outsum <- parsed_content( # API call
          GET("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", 
              query=querysum))$doc$children$eSummaryResult 
        names <- sapply(getNodeSet(outsum[[1]], "//Item"), xmlGetAttr, name="Name") # gets names of values in summary
        length_ <- as.numeric(sapply(getNodeSet(outsum, "//Item"), xmlValue)[str_detect(names, "Length")]) # gets seq lengths
        gis <- as.numeric(sapply(getNodeSet(outsum, "//Item"), xmlValue)[str_detect(names, "Gi")]) # gets GI numbers
        spnames <- sapply(getNodeSet(outsum, "//Item"), xmlValue)[str_detect(names, "Title")] # gets seq lengths # get spp names
        df <- data.frame(gis=gis, length=length_, spnames=laply(spnames, c)) # makes data frame
        gisuse <- df[which.max(x=df$length),] # picks longest sequnence length
        if(nrow(gisuse)>1){gisuse <- gisuse[sample(nrow(gisuse), 1), ]} else 
        {gisuse <- gisuse}
        
        ## Get sequence from previous
        message("...retrieving sequence...")
        queryseq <- list(db = "sequences", id = gisuse[,1], rettype = "fasta", retmode = "text")
        outseq <- parsed_content(
          GET("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", query = queryseq))  
        seq <- str_replace_all(str_split(str_replace(outseq, "\n", "<<<"), "<<<")[[1]][[2]], "\n", "")
        accessnum <- str_split(outseq, "\\|")[[1]][4]
        outt <- list(taxon_name, as.character(gisuse[,3]), gisuse[,1], accessnum, gisuse[,2], seq)
      }
    }
  } else
  {
    ids <- xpathApply(out, "//IdList//Id") # Get sequence IDs in list
    ids_ <- as.numeric(sapply(ids, xmlValue))  # Get sequence ID values
    
    ## For each species = get GI number with longest sequence
    message("...retrieving sequence ID with longest sequence length...")
    querysum <- list(db = "nucleotide", id = paste(ids_, collapse=" ")) # construct query for species
    outsum <- parsed_content( # API call
      GET("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", 
          query=querysum))$doc$children$eSummaryResult 
    names <- sapply(getNodeSet(outsum[[1]], "//Item"), xmlGetAttr, name="Name") # gets names of values in summary
    length_ <- as.numeric(sapply(getNodeSet(outsum, "//Item"), xmlValue)[str_detect(names, "Length")]) # gets seq lengths
    gis <- as.numeric(sapply(getNodeSet(outsum, "//Item"), xmlValue)[str_detect(names, "Gi")]) # gets GI numbers
    spnames <- sapply(getNodeSet(outsum, "//Item"), xmlValue)[str_detect(names, "Title")] # gets seq lengths # get spp names
    df <- data.frame(gis=gis, length=length_, spnames=laply(spnames, c)) # makes data frame
    gisuse <- df[which.max(x=df$length),] # picks longest sequnence length
    if(nrow(gisuse)>1){gisuse <- gisuse[sample(nrow(gisuse), 1), ]} else 
    {gisuse <- gisuse}
    
    ## Get sequence from previous
    message("...retrieving sequence...")
    queryseq <- list(db = "sequences", id = gisuse[,1], rettype = "fasta", retmode = "text")
    outseq <- parsed_content(
      GET("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", query = queryseq))	
    seq <- str_replace_all(str_split(str_replace(outseq, "\n", "<<<"), "<<<")[[1]][[2]], "\n", "")
    # 		accessnum <- str_split(outseq, "\\|")[[1]][
    # 			grep("[A-Za-z]{1}[0-9]{5}|[A-Za-z]{2}[0-9]{6}", str_split(outseq, "\\|")[[1]])
    # 		]
    accessnum <- str_split(outseq, "\\|")[[1]][4]
    outt <- list(taxon_name, as.character(gisuse[,3]), gisuse[,1], accessnum, gisuse[,2], seq)
  }
  message("...done.")
  #   return( outt )
  outoutout <- data.frame(outt)
  names(outoutout) <- NULL
  write.table(outoutout, file = "beeseqsoutput.txt", append=T, row.names=F)
}