clipSeq <- function(vdf = NULL) {

    # set up extra columns with null values
    vdf$softclip <- rep(0, nrow(vdf))
    vdf$leftclip <- rep(0, nrow(vdf))
    vdf$rightclip <- rep(0, nrow(vdf))
    vdf$clipseq   <- as.vector(vdf$seq)

    # we take a copy of the cigar string
    # and remove hard clipping information
    # for convenience
    vdf$cigar2 <- gsub("\\d+H", "", vdf$cigar)

    # Record the column indices for cigar columns
    cidx <- grep("cigar", colnames(vdf))[1]
    cidx2 <- grep("cigar", colnames(vdf))[2]

    # get the subset of rows with soft-clipping
    clipr <- grep("S", vdf$cigar)
    svdf <- vdf[clipr, ]

    # calculate the total number of soft clipped bases
    # by summing up numbers that are adjacent to an S
    # in the cigar string
    #svdf$softclip <- apply(svdf, 1, function(x) {
    #   sum(as.numeric(gsubfn::strapply(as.character(x[cidx]), "(\\d+)S",
    #       simplify=c)))
    #})

    # calculate the number of bases soft clipped from
    # the left of the read by extracting the number
    # adjacent to any S at the beginning of the string
    svdf$leftclip <- apply(svdf, 1, function(x) {
      sum(as.numeric(strapply(as.character(x[cidx2]), "^(\\d+)S", simplify = c)))
    })

    # calculate the number of bases soft clipped from
    # the right of the read by extracting the number
    # adjacent to any S at the end of the string
    svdf$rightclip <- apply(svdf, 1, function(x) {
      sum(as.numeric(strapply(as.character(x[cidx2]), "(\\d+)S$", simplify = c)))
    })

    # extract the sequence in the middle of leftclip and rightclip
    # this is the sequence that has actually aligned
    svdf$clipseq <- substr(svdf$seq,
                        as.numeric(svdf$leftclip) + 1,
                        nchar(as.vector(svdf$seq)) - as.numeric(svdf$rightclip))

    # substitute in the values
    # vdf$softclip[clipr] <- svdf$softclip
    vdf$leftclip[clipr] <- svdf$leftclip
    vdf$rightclip[clipr] <- svdf$rightclip
    vdf$clipseq[clipr] <- svdf$clipseq

    # calculate the length of the clipped sequence
    vdf$cliplen <- nchar(vdf$clipseq)

    # remove the second cigar string
    # we calculated above
    vdf <- vdf[, -cidx2]

    # return the data
    return(vdf)

}
