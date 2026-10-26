make_chromosome_test_bam <- function(build, mt_name = "chrM", mt_length = NULL, mate_offset = 100L, header_only = NULL) {
  reference <- switch(build,
    hg19 = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    hg38 = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  lengths <- GenomeInfoDb::seqlengths(reference)[c(paste0("chr", 1:22), "chrX", "chrY", "chrM")]
  # hg19 contains distinct chrM (16571) and chrMT (16569) sequences.
  if (mt_name %in% GenomeInfoDb::seqlevels(reference)) {
    lengths["chrM"] <- GenomeInfoDb::seqlengths(reference)[mt_name]
  }
  if (!is.null(mt_length)) lengths["chrM"] <- mt_length
  names(lengths)[names(lengths) == "chrM"] <- mt_name
  chromosomes <- c("chr1", "chr22", "chrX", "chrY", mt_name)
  if (!is.null(header_only)) {
    lengths <- lengths[header_only]
    chromosomes <- intersect(chromosomes, header_only)
  }
  bases <- paste(rep("A", 50), collapse = "")
  qualities <- paste(rep("I", 50), collapse = "")
  records <- unlist(lapply(seq_along(chromosomes), function(i) {
    chr <- chromosomes[i]
    c(paste(paste0("pair", i), 99, chr, 1000, 60, "50M", "=", 1000 + mate_offset, mate_offset + 50L,
            bases, qualities, "NM:i:0", "MD:Z:50", sep = "\t"),
      paste(paste0("pair", i), 147, chr, 1000 + mate_offset, 60, "50M", "=", 1000, -(mate_offset + 50L),
            bases, qualities, "NM:i:0", "MD:Z:50", sep = "\t"))
  }))
  stem <- tempfile("chromosomes_")
  sam <- paste0(stem, ".sam")
  writeLines(c("@HD\tVN:1.6\tSO:coordinate",
               paste0("@SQ\tSN:", names(lengths), "\tLN:", lengths), records), sam)
  Rsamtools::asBam(sam, destination = stem, overwrite = TRUE)
  paste0(stem, ".bam")
}

test_that("hg19 and hg38 chromosome filters retain exactly the selected pairs", {
  for (build in c("hg19", "hg38")) {
    bam <- make_chromosome_test_bam(build)
    for (curate in c(TRUE, FALSE)) {
      defaults <- readBam(bam, genome_label = build, curate_start_and_end = curate)
      expect_identical(as.character(GenomicRanges::seqnames(defaults)), c("chr1", "chr22"))
      for (selection in list("chrX", "chrY", "chrM", c("chrX", "chrY"), c("chr1", "chrX", "chrY", "chrM"))) {
        frag <- readBam(bam, genome_label = build, chromosome_to_keep = selection,
                        curate_start_and_end = curate)
        expect_setequal(as.character(GenomicRanges::seqnames(frag)), selection)
        expect_length(frag, length(selection))
        expect_true(all(BiocGenerics::width(frag) == 150L))
      }
      all_frag <- readBam(bam, genome_label = build, chromosome_to_keep = FALSE,
                         curate_start_and_end = curate)
      expect_length(all_frag, 5L)
      expect_setequal(as.character(GenomicRanges::seqnames(all_frag)), c("chr1", "chr22", "chrX", "chrY", "chrM"))
      expect_true(GenomeInfoDb::isCircular(all_frag)["chrM"])
      expect_identical(unname(GenomeInfoDb::seqlengths(all_frag)["chrM"]), if (build == "hg19") 16571L else 16569L)
    }
  }
})

test_that("mitochondrial aliases are resolved before BAM filtering and curation", {
  for (build in c("hg19", "hg38")) {
    for (bam_mt in c("chrM", "chrMT", "MT", "M")) {
      bam <- make_chromosome_test_bam(build, mt_name = bam_mt)
      for (requested_mt in c("chrM", "chrMT", "MT", "M")) {
        for (curate in c(TRUE, FALSE)) {
          frag <- readBam(bam, genome_label = build, chromosome_to_keep = requested_mt,
                          curate_start_and_end = curate)
          expect_length(frag, 1L)
          expected_name <- if (build == "hg19" && bam_mt == "chrMT") "chrMT" else "chrM"
          expect_identical(as.character(GenomicRanges::seqnames(frag)), expected_name)
          expect_identical(names(frag), "pair5")
          expect_true(GenomeInfoDb::isCircular(frag)[expected_name])
        }
      }
    }
  }
})

test_that("mitochondrial aliases do not bypass reference length validation", {
  for (build in c("hg19", "hg38")) {
    wrong_length <- if (build == "hg19") 16569L else 16571L
    for (bam_mt in c("chrM", "MT")) {
      bam <- make_chromosome_test_bam(build, mt_name = bam_mt, mt_length = wrong_length)
      expect_error(readBam(bam, genome_label = build, chromosome_to_keep = "chrM"), "incompatible seqlengths")
      expect_length(readBam(bam, genome_label = build), 2L)
    }
  }
})

test_that("invalid chromosome selections identify the requested missing chromosome", {
  bam <- make_chromosome_test_bam("hg38")
  expect_error(readBam(bam, genome_label = "hg38", chromosome_to_keep = "chrZ"), "not found in the BAM header: chrZ")
  expect_error(readBam(bam, genome_label = "hg38", chromosome_to_keep = character()), "non-empty character vector or FALSE")
  expect_error(readBam(bam, genome_label = "hg38", chromosome_to_keep = NA_character_), "non-empty character vector or FALSE")
})

test_that("default selection supports BAM headers containing only a regional subset", {
  bam <- make_chromosome_test_bam("hg38", header_only = "chr1")
  frag <- readBam(bam, genome_label = "hg38")
  expect_length(frag, 1L)
  expect_identical(as.character(GenomicRanges::seqnames(frag)), "chr1")
  expect_error(readBam(bam, genome_label = "hg38", chromosome_to_keep = c("chr1", "chrZ")), "not found in the BAM header: chrZ")
})

test_that("mutation-aware chromosome selection supports aliases and FALSE", {
  for (mt_name in c("chrM", "MT")) {
    bam <- make_chromosome_test_bam("hg38", mt_name = mt_name, mate_offset = 10L)
    mutations <- tempfile(fileext = ".tsv")
    utils::write.table(data.frame(chr = c("chr1", "chr22", "chrX", "chrY", mt_name),
                                 pos = 1020L, ref = "A", alt = "C"),
                       mutations, sep = "\t", row.names = FALSE, quote = FALSE)
    for (mut_only in c(TRUE, FALSE)) {
      for (selection in list(c("chrX", "chrM"), FALSE)) {
        frag <- readBam(bam, genome_label = "hg38", chromosome_to_keep = selection,
                        mutation_file = mutations, mut_fragments_only = mut_only,
                        galp_what = c("cigar", "mapq", "isize", "seq", "qual"))
        expected <- if (isFALSE(selection)) c("chr1", "chr22", "chrX", "chrY", "chrM") else selection
        expect_length(frag, length(expected))
        expect_setequal(as.character(GenomicRanges::seqnames(frag)), expected)
        expect_true(all(BiocGenerics::start(frag) == 1000L))
        expect_true(all(BiocGenerics::width(frag) == 60L))
        expect_true(all(frag$mutation_status == "REF:read_pair_concordant:A/A"))
      }
    }
  }
})
