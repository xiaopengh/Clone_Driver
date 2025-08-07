# ==============================================================================

get5utr <- function(mrnaGtf) {
  
  # split gtf by strand
  strandGtf <- split(mrnaGtf, mrnaGtf$strand)
  
  utr5Gtf <- lapply(
    names(strandGtf),
    function(s) {
      cdsGtf <- strandGtf[[s]][type == "CDS"]
      utrGtf <- strandGtf[[s]][type == "UTR"]
      if (s == "+") {
        firstcds <- tapply(cdsGtf$start, cdsGtf$transcript_id, min)
        return(utrGtf[utrGtf$end < firstcds[utrGtf$transcript_id]])
      }
      firstcds <- tapply(cdsGtf$end, cdsGtf$transcript_id, max)
      return(utrGtf[utrGtf$start > firstcds[utrGtf$transcript_id]])
    }
  )
  
  return(rbindlist(utr5Gtf))
  
}

get3utr <- function(mrnaGtf) {
  
  # split gtf by strand
  strandGtf <- split(mrnaGtf, mrnaGtf$strand)
  
  utr3Gtf <- lapply(
    names(strandGtf),
    function(s) {
      cdsGtf <- strandGtf[[s]][type == "CDS"]
      utrGtf <- strandGtf[[s]][type == "UTR"]
      if (s == "+") {
        lastcds <- tapply(cdsGtf$end, cdsGtf$transcript_id, max)
        return(utrGtf[utrGtf$start > lastcds[utrGtf$transcript_id]])
      }
      lastcds <- tapply(cdsGtf$start, cdsGtf$transcript_id, min)
      return(utrGtf[utrGtf$end < lastcds[utrGtf$transcript_id]])
    }
  )
  
  return(data.table::rbindlist(utr3Gtf))
  
}