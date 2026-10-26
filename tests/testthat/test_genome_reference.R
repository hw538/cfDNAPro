test_that("built-in and custom readers work without reference downloads", {
  for (label in c("hg19", "hg38", "hg38-NCBI")) {
    g <- cfDNAPro:::.cfDNAPro_builtin_genome(label)
    fixture <- reference_fixture(g)
    without_reference_downloads({
      regular <- readBam(fixture$bam, genome_label = label, chromosome_to_keep = fixture$chr)
      custom <- readBam(fixture$bam, genome_label = "user_define", genome = g,
                        chromosome_to_keep = fixture$chr)
      expect_identical(as.data.frame(custom), as.data.frame(regular))
      expect_identical(cfDNAPro:::get_genome_reference(custom), g)
      galp <- readGALP(fixture$bam, genome_label = "user_define", genome = g,
                       chromosome_to_keep = fixture$chr)
      expect_length(galp, 1L)
      expect_identical(cfDNAPro:::get_genome_reference(galp), g)
      counts <- read_bam_insert_metrics(bamfile = fixture$bam, genome_label = "user_define",
                                        genome = g, chromosome_to_keep = fixture$chr,
                                        isize_min = 60L, isize_max = 60L)
      expect_equal(counts$All_Reads.fr_count, 1)
      expect_equal(callMotif(custom, genome = g), callMotif(regular))
    })
  }
})

test_that("packaged custom references survive RDS and a fresh R process", {
  g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  fixture <- reference_fixture(g)
  frag <- readBam(fixture$bam, genome_label = "user_define", genome = g,
                  chromosome_to_keep = fixture$chr)
  path <- tempfile(fileext = ".rds")
  saveRDS(frag, path)
  restored <- readRDS(path)
  expect_identical(cfDNAPro:::get_genome_reference(restored), g)
  expect_equal(callMotif(restored), callMotif(frag, genome = g))
  expect_lt(file.info(path)$size, 100000)
  script <- tempfile(fileext = ".R")
  output <- tempfile(fileext = ".rds")
  writeLines(c("suppressPackageStartupMessages(library(cfDNAPro))",
               "args <- commandArgs(TRUE)",
               "saveRDS(callMotif(readRDS(args[1])), args[2])"), script)
  run <- system2(file.path(R.home("bin"), "Rscript"),
                 c("--vanilla", shQuote(script), shQuote(path), shQuote(output)),
                 stdout = TRUE, stderr = TRUE)
  expect_null(attr(run, "status"))
  expect_true(file.exists(output))
  if (file.exists(output)) expect_equal(readRDS(output), callMotif(frag, genome = g))
})

test_that("ad-hoc references keep their identity and require an explicit object after RDS", {
  g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  custom <- g
  GenomeInfoDb::genome(custom) <- "local_custom_assembly"
  names <- GenomeInfoDb::seqnames(custom)
  names[1] <- "local_contig"
  GenomeInfoDb::seqnames(custom) <- names
  fixture <- reference_fixture(custom)
  without_reference_downloads({
    frag <- readBam(fixture$bam, genome_label = "user_define", genome = custom,
                    chromosome_to_keep = fixture$chr)
    expect_identical(as.character(GenomicRanges::seqnames(frag)), "local_contig")
    expect_identical(unique(GenomeInfoDb::genome(frag)), "local_custom_assembly")
    expect_null(S4Vectors::metadata(frag)$cfDNAPro_reference$package)
    path <- tempfile(fileext = ".rds")
    saveRDS(frag, path)
    restored <- readRDS(path)
    expect_error(callMotif(restored), "pass `genome` explicitly")
    motif <- callMotif(restored, genome = custom)
    expected <- as.character(BSgenome::getSeq(custom, fixture$chr,
                                             start = fixture$start, end = fixture$start + 2L))
    expect_equal(as.character(motif$motif[motif$n == 1]), expected)
  })
})

test_that("explicit references and restored provenance cannot silently change assemblies", {
  hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  fixture <- reference_fixture(hg19)
  frag <- readBam(fixture$bam, chromosome_to_keep = fixture$chr)
  expect_error(callMotif(frag, genome = hg38), "incompatible sequence lengths")
  expect_error(callMotif(frag, genome = "hg19"), "BSgenome")
  changed <- frag
  md <- S4Vectors::metadata(changed)
  md$cfDNAPro_reference$package_version <- "0.0.0"
  S4Vectors::metadata(changed) <- md
  expect_error(callMotif(changed), "version has changed")
  expect_equal(callMotif(changed, genome = hg19), callMotif(frag))
  md$cfDNAPro_reference$seqinfo$seqlengths[1] <- 1L
  S4Vectors::metadata(changed) <- md
  expect_error(callMotif(changed, genome = hg19), "recorded reference metadata")
  legacy <- frag
  S4Vectors::metadata(legacy) <- list()
  expect_identical(cfDNAPro:::get_genome_reference(legacy), hg19)
  GenomeInfoDb::genome(legacy) <- "unknown_assembly"
  expect_error(cfDNAPro:::get_genome_reference(legacy), "Cannot infer")
  fixture38 <- reference_fixture(hg38)
  legacy38 <- readBam(fixture38$bam, genome_label = "hg38", chromosome_to_keep = fixture38$chr)
  S4Vectors::metadata(legacy38) <- list()
  GenomeInfoDb::genome(legacy38) <- "GRCh38"
  expect_identical(cfDNAPro:::get_genome_reference(legacy38, genome = hg38), hg38)
})

test_that("mutation TSV chromosome names and T alleles stay character data", {
  path <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(chr = "1", pos = 123L, ref = "T", alt = "C"),
                     path, sep = "\t", quote = FALSE, row.names = FALSE)
  mutations <- cfDNAPro:::read_mutation_file(path)
  expect_identical(mutations$chr, "1")
  expect_identical(mutations$ref, "T")
  expect_identical(mutations$pos, 123L)
})

test_that("reader validation and existing argument positions are preserved", {
  g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  fixture <- reference_fixture(g)
  expect_identical(names(formals(readBam))[7], "outdir")
  expect_identical(names(formals(readGALP))[6], "outdir")
  expect_identical(names(formals(read_bam_insert_metrics))[6], "outdir")
  expect_length(readBam(fixture$bam, TRUE, TRUE, fixture$chr, 1, "hg38", FALSE), 1L)
  expect_length(readBam(fixture$bam, genome = "hg38", chromosome_to_keep = fixture$chr), 1L)
  expect_error(readBam(fixture$bam, genome_label = "user_define"), "BSgenome")
  expect_error(readBam(fixture$bam, genome_label = "user_define", genome = "hg38"), "BSgenome")
  for (bad in list(character(), c("hg19", "hg38"), NA_character_, 1, list("hg19"))) {
    expect_error(readBam(fixture$bam, genome_label = bad), "genome_label must be")
  }
  expect_warning(readBam(fixture$bam, genome_label = "hg38", genome = g,
                          chromosome_to_keep = fixture$chr), "ignored")
})

test_that("custom references propagate through mutation context and table export", {
  g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  fixture <- reference_fixture(g, mutation = TRUE)
  without_reference_downloads({
    frag <- readBam(fixture$bam, genome_label = "user_define", genome = g,
                    chromosome_to_keep = fixture$chr, mutation_file = fixture$mutation_file,
                    mut_fragments_only = TRUE, galp_what = c("cigar", "mapq", "isize", "seq", "qual"))
    summary <- callTrinucleotide(frag)
    expect_equal(summary, callTrinucleotide(frag, genome = g))
    context <- as.character(BSgenome::getSeq(g, fixture$chr, start = fixture$pos - 1L, end = fixture$pos + 1L))
    expected <- paste0(substr(context, 1, 1), "[", fixture$ref, ">", fixture$alt, "]", substr(context, 3, 3))
    expect_identical(summary$SBS96, expected)
    output <- tempfile(fileext = ".tsv")
    expect_identical(writeMutTable(frag, output, genome = g), output)
    expect_true(file.exists(output))
    expect_equal(utils::read.delim(output)$SBS96, expected)
  })
})
