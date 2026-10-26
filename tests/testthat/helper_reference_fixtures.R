reference_fixture <- function(genome, mutation = FALSE) {
  chr <- GenomeInfoDb::seqlevels(genome)[1]
  start <- 10000001L
  # Keep this fixture's SNV in the pyrimidine orientation; the tests exercise
  # reference selection and persistence, not the existing SBS normalization.
  repeat {
    pos <- start + 20L
    ref <- as.character(BSgenome::getSeq(genome, chr, start = pos, end = pos))
    if (ref %in% c("C", "T")) break
    start <- start + 1L
  }
  alt <- if (ref == "C") "A" else "G"
  first <- as.character(BSgenome::getSeq(genome, chr, start = start, end = start + 49L))
  last <- as.character(BSgenome::getSeq(genome, chr, start = start + 10L, end = start + 59L))
  if (mutation) {
    substr(first, 21, 21) <- alt
    substr(last, 11, 11) <- alt
  }
  md1 <- if (mutation) paste0("20", ref, "29") else "50"
  md2 <- if (mutation) paste0("10", ref, "39") else "50"
  quality <- paste(rep("I", 50), collapse = "")
  stem <- tempfile("reference_fixture_")
  lengths <- GenomeInfoDb::seqlengths(genome)
  records <- c(
    paste("pair", 99, chr, start, 60, "50M", "=", start + 10L, 60, first, quality,
          paste0("NM:i:", as.integer(mutation)), paste0("MD:Z:", md1), sep = "\t"),
    paste("pair", 147, chr, start + 10L, 60, "50M", "=", start, -60, last, quality,
          paste0("NM:i:", as.integer(mutation)), paste0("MD:Z:", md2), sep = "\t"))
  writeLines(c("@HD\tVN:1.6\tSO:coordinate",
               paste0("@SQ\tSN:", names(lengths), "\tLN:", lengths), records), paste0(stem, ".sam"))
  Rsamtools::asBam(paste0(stem, ".sam"), destination = stem, overwrite = TRUE)
  mutation_file <- paste0(stem, ".tsv")
  utils::write.table(data.frame(chr = chr, pos = pos, ref = ref, alt = alt),
                     mutation_file, sep = "\t", quote = FALSE, row.names = FALSE)
  list(bam = paste0(stem, ".bam"), chr = chr, start = start, pos = pos, ref = ref,
       alt = alt, mutation_file = mutation_file)
}

without_reference_downloads <- function(expr) {
  forbidden <- function(...) stop("Network reference lookup is forbidden in this test.")
  testthat::local_mocked_bindings(getChromInfoFromUCSC = forbidden,
                                 getChromInfoFromNCBI = forbidden,
                                 .package = "GenomeInfoDb", .env = environment())
  testthat::expect_error(GenomeInfoDb::getChromInfoFromUCSC("hg38"), "forbidden")
  testthat::expect_error(GenomeInfoDb::getChromInfoFromNCBI("GRCh38"), "forbidden")
  force(expr)
}
