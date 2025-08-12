infoExtract <- function(gtf, infoColName, key, type) {
  
  n <- length(key)
  if (n != length(type)) stop("key & type must be of same length")
  
  # parse info into data.frame keeping track of feature index
  info <- strsplit(gtf[[infoColName]], ";")
  infodt <- data.table::data.table(
    info = unlist(info),
    idx = unlist(mapply(function(v, i) rep(i, length(v)), info, 1:length(info)))
  )
  infodt$info <- gsub("\"", "", infodt$info)
  infodt$info <- trimws(infodt$info)
  
  # format keys for matching
  fkey <- ifelse(type == "tag", paste("tag", key), paste0("^", key, " "))
  
  # match
  for (i in 1:n) {
    
    if (type[i] == "tag") {
      
      r <- getTag(infodt, fkey[i], nrow(gtf))
      
    } else {
      
      r <- getMulti(infodt, fkey[i], nrow(gtf))
      
    }
    gtf[[ key[i] ]] <- r
    
  }
  
  return(gtf)
  
}

getMulti <- function(infodt, multi, nfeatures) {
  
  # get values for the queried multilevel key for each feature in gtf
  r <- rep(NA_character_, nfeatures)
  ismulti <- grepl(multi, infodt$info)
  r[infodt$idx[ismulti]] <- gsub(multi, "", infodt$info[ismulti])
  
  return(r)
  
}

getTag <- function(infodt, tag, nfeatures) {
  
  # get values for the queried binary key for each feature in gtf
  r <- rep(0L, nfeatures)
  r[infodt$idx[infodt$info == tag]] <- 1L
  
  return(r)
  
}

pvarCDSannotate <- function(cdsgtf, .genome) {
  
  # ensure order of this chromosome
  cdsgtf <- cdsgtf[order(Start_Position)]
  
  # get coding sites, their transcript of origin and strand
  gtfRanges <- GenomicRanges::GRanges(cdsgtf$Chromosome, IRanges::IRanges(cdsgtf$Start_Position, cdsgtf$End_Position))
  siteRanges <- unlist(GenomicRanges::slidingWindows(gtfRanges, 1))
  siteRanges$transcript_id <- rep(cdsgtf$transcript_id, width(gtfRanges))
  siteRanges$codingStrand <- rep(cdsgtf$Strand, GenomicRanges::width(gtfRanges))
  
  # split sites by transcript, if they are in negative strand reverse order so that everything is in transcription order
  siteRanges <- split(siteRanges, siteRanges$transcript_id)
  reverse <- sapply(siteRanges, function(gr) gr$codingStrand[1]) == "-"
  siteRanges[reverse] <- S4Vectors::endoapply(siteRanges[reverse], rev)
  
  # get relative position of each site within its transcript CDS
  relpos <- lapply(IRanges::elementNROWS(siteRanges), function(l) 1:l)
  siteRanges <- unlist(siteRanges)
  siteRanges$relpos <- unlist(relpos)
  
  # get nucleotides in the + strand
  siteRanges$ref <- .genome[siteRanges]
  
  # get nucleotides in the coding strand
  siteRanges$codingRef <- siteRanges$ref
  reverse <- siteRanges$codingStrand == "-"
  siteRanges$codingRef[reverse] <- Biostrings::reverseComplement(siteRanges$codingRef[reverse])
  
  # determine pyrimidine strand and get pyrimidine
  siteRanges$pyrimidineStrand <- ifelse(as.character(siteRanges$ref) %in% c("C", "T"), "+", "-")
  siteRanges$pyrimidine <- siteRanges$ref
  reverse <- siteRanges$pyrimidineStrand == "-"
  siteRanges$pyrimidine[reverse] <- Biostrings::reverseComplement(siteRanges$pyrimidine[reverse])
  
  # get the codon and aminoacid that each site influences
  siteRanges$frame <- ((siteRanges$relpos + 2L) %% 3L) + 1L # gives 1 if in frame, 2 for one base out of frame etc
  siteRanges$codon <- rep(
    xscat(
      siteRanges$codingRef[siteRanges$frame == 1],
      siteRanges$codingRef[siteRanges$frame == 2],
      siteRanges$codingRef[siteRanges$frame == 3]
    ),
    each = 3
  )
  siteRanges$wtAA <- Biostrings::translate(siteRanges$codon, no.init.codon = TRUE)
  
  # get all possible mutations in the examined sites oriented by the pyrimidine
  pmuts <- list(C = c("A", "G", "T"), T = c("A", "C", "G"))[as.character(siteRanges$pyrimidine)]
  pmutRanges <- siteRanges[rep(1:length(siteRanges), each = 3)]
  pmutRanges$pyrimidineMut <- Biostrings::DNAStringSet(unlist(pmuts))
  reverse <- rep(reverse, each = 3)
  pmutRanges$mut <- pmutRanges$pyrimidineMut
  pmutRanges$mut[reverse] <- Biostrings::reverseComplement(pmutRanges$mut[reverse])
  
  # get the mutation in the coding strand according to pyrimidine strand
  reverse <- pmutRanges$pyrimidineStrand != pmutRanges$codingStrand
  pmutRanges$codingMut <- pmutRanges$pyrimidineMut
  pmutRanges$codingMut[reverse] <- Biostrings::reverseComplement(pmutRanges$codingMut[reverse])
  
  # get the mutated codon and mutated aminoacid according to mutation in coding strand
  pmutRanges$mutCodon <- pmutRanges$codon
  Biostrings::subseq(pmutRanges$mutCodon, pmutRanges$frame, pmutRanges$frame) <- pmutRanges$codingMut
  pmutRanges$mutAA <- Biostrings::translate(pmutRanges$mutCodon, no.init.codon = TRUE)
  
  # determine the type of mutation
  pmutRanges$type <- rep("missense", length(siteRanges))
  pmutRanges$type[pmutRanges$wtAA == pmutRanges$mutAA] <- "syn"
  pmutRanges$type[pmutRanges$wtAA != "*" & pmutRanges$mutAA == "*"] <- "nonsense"
  pmutRanges$type[pmutRanges$wtAA == "*" & pmutRanges$mutAA != "*"] <- "nonstop"
  
  # convert to data.table
  names(pmutRanges) <- 1:length(pmutRanges)
  pmutdt <- data.table::data.table(as.data.frame(pmutRanges))
  pmutdt[, c("width", "strand", "end") := NULL]
  names(pmutdt)[2] <- "position"
  
  return(pmutdt)
  
}