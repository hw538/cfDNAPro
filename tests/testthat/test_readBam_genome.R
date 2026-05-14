test_that("get_genome_reference uses explicit genome argument", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
  g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  out <- cfDNAPro:::get_genome_reference(gr, genome = g)
  expect_identical(out, g)
})

test_that("get_genome_reference rejects non-BSgenome explicit genome", {
  gr <- GenomicRanges::GRanges()
  expect_error(
    cfDNAPro:::get_genome_reference(gr, genome = "not-a-BSgenome"),
    "BSgenome"
  )
})

test_that("get_genome_reference restores user_define from package metadata", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
  g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  gr <- GenomicRanges::GRanges()
  S4Vectors::metadata(gr) <- list(
    cfDNAPro_genome_label = "user_define",
    cfDNAPro_bsgenome_pkg = "BSgenome.Hsapiens.UCSC.hg19",
    cfDNAPro_bsgenome_symbol = "BSgenome.Hsapiens.UCSC.hg19"
  )
  out <- cfDNAPro:::get_genome_reference(gr)
  expect_identical(out, g)
})

test_that("get_genome_reference errors when user_define cannot be restored", {
  gr <- GenomicRanges::GRanges()
  S4Vectors::metadata(gr) <- list(cfDNAPro_genome_label = "user_define")
  expect_error(cfDNAPro:::get_genome_reference(gr), "user_define")
})
