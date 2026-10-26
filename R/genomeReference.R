# Helpers shared by BAM readers and sequence-based analyses. Reference objects
# stay outside fragment RDS files; only their small provenance manifest is saved.

.cfDNAPro_builtin_genome <- function(label) {
  switch(label,
    hg19 = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    hg38 = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    `hg38-NCBI` = BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38,
    stop("Unknown built-in genome label.", call. = FALSE))
}

.cfDNAPro_resolve_genome <- function(genome_label, genome = NULL,
                                    label_missing = FALSE) {
  builtins <- c("hg19", "hg38", "hg38-NCBI")
  # Before genome was a formal argument, genome="hg38" partially matched
  # genome_label. Keep that existing calling convention working.
  if (label_missing && is.character(genome) && length(genome) == 1L &&
      !is.na(genome) && genome %in% builtins) {
    genome_label <- genome
    genome <- NULL
  }
  if (!is.character(genome_label) || length(genome_label) != 1L ||
      is.na(genome_label) || !genome_label %in% c(builtins, "user_define")) {
    stop("genome_label must be one of 'hg19', 'hg38', 'hg38-NCBI', or 'user_define'.",
         call. = FALSE)
  }
  if (genome_label == "user_define") {
    if (!methods::is(genome, "BSgenome")) {
      stop("When genome_label is 'user_define', supply a BSgenome object via `genome`.",
           call. = FALSE)
    }
    return(list(label = genome_label, genome = genome))
  }
  if (!is.null(genome)) {
    if (!methods::is(genome, "BSgenome")) {
      stop("`genome` must be a BSgenome object; use genome_label for a built-in label.",
           call. = FALSE)
    }
    warning("Argument `genome` is ignored unless genome_label is 'user_define'.",
            call. = FALSE)
  }
  list(label = genome_label, genome = .cfDNAPro_builtin_genome(genome_label))
}

.cfDNAPro_reference_manifest <- function(genome) {
  info <- GenomeInfoDb::seqinfo(genome)
  result <- data.frame(
    seqnames = GenomeInfoDb::seqlevels(info),
    seqlengths = unname(GenomeInfoDb::seqlengths(info)),
    isCircular = unname(GenomeInfoDb::isCircular(info)),
    genome = unname(GenomeInfoDb::genome(info)), stringsAsFactors = FALSE)
  result <- result[order(result$seqnames, method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  result
}

.cfDNAPro_attach_reference <- function(object, label, genome) {
  package <- NULL
  version <- NULL
  if ("pkgname" %in% methods::slotNames(genome)) {
    candidate <- methods::slot(genome, "pkgname")
    if (is.character(candidate) && length(candidate) == 1L &&
        !is.na(candidate) && startsWith(candidate, "BSgenome.")) {
      installed <- tryCatch(BSgenome::getBSgenome(candidate), error = function(e) NULL)
      # Modified/ad-hoc objects must be supplied again after an RDS round trip.
      # A package name alone must never restore a different reference object.
      if (!is.null(installed) && identical(installed, genome)) {
        package <- candidate
        version <- as.character(utils::packageVersion(package))
      }
    }
  }
  provenance <- list(schema_version = 1L, label = label, package = package,
                     package_version = version,
                     seqinfo = .cfDNAPro_reference_manifest(genome))
  md <- S4Vectors::metadata(object)
  md$cfDNAPro_reference <- provenance
  S4Vectors::metadata(object) <- md
  object
}

.cfDNAPro_validate_reference <- function(object, genome) {
  if (!methods::is(genome, "BSgenome")) {
    stop("`genome` must be a BSgenome object when provided.", call. = FALSE)
  }
  used <- as.character(unique(GenomeInfoDb::seqnames(object)))
  used <- used[!is.na(used)]
  missing <- setdiff(used, GenomeInfoDb::seqlevels(genome))
  if (length(missing)) {
    stop("The supplied genome is missing fragment chromosomes: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  observed <- GenomeInfoDb::seqlengths(object)[used]
  expected <- GenomeInfoDb::seqlengths(genome)[used]
  bad <- used[!is.na(observed) & !is.na(expected) & observed != expected]
  if (length(bad)) {
    stop("The supplied genome has incompatible sequence lengths for: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  observed_labels <- unique(GenomeInfoDb::genome(object)[used])
  expected_labels <- unique(GenomeInfoDb::genome(genome)[used])
  observed_labels <- observed_labels[!is.na(observed_labels)]
  expected_labels <- expected_labels[!is.na(expected_labels)]
  assembly_aliases <- c(hg19 = "GRCh37", GRCh37 = "GRCh37", `GRCh37.p13` = "GRCh37",
                        hg38 = "GRCh38", GRCh38 = "GRCh38", `hg38-NCBI` = "GRCh38")
  normalize_labels <- function(labels) {
    known <- labels %in% names(assembly_aliases)
    labels[known] <- unname(assembly_aliases[labels[known]])
    labels
  }
  if (length(observed_labels) && length(expected_labels) &&
      !setequal(normalize_labels(observed_labels), normalize_labels(expected_labels))) {
    stop("The supplied genome has incompatible assembly labels.", call. = FALSE)
  }
  recorded <- S4Vectors::metadata(object)$cfDNAPro_reference
  if (!is.null(recorded) &&
      !identical(recorded$seqinfo, .cfDNAPro_reference_manifest(genome))) {
    stop("The supplied genome does not match the recorded reference metadata.",
         call. = FALSE)
  }
  genome
}

get_genome_reference <- function(frag_obj_mut, genome = NULL) {
  if (!is.null(genome)) {
    return(.cfDNAPro_validate_reference(frag_obj_mut, genome))
  }
  md <- S4Vectors::metadata(frag_obj_mut)
  recorded <- md$cfDNAPro_reference
  if (!is.null(recorded)) {
    if (!is.list(recorded) || !identical(recorded$schema_version, 1L)) {
      stop("Unsupported reference metadata; pass `genome` explicitly.", call. = FALSE)
    }
    package <- recorded$package
    if (is.null(package) || !is.character(package) || length(package) != 1L ||
        is.na(package) || !startsWith(package, "BSgenome.")) {
      stop("This custom reference cannot be restored from a package; pass `genome` explicitly.",
           call. = FALSE)
    }
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Reference package '", package, "' is unavailable; install it or pass `genome` explicitly.",
           call. = FALSE)
    }
    if (!identical(as.character(utils::packageVersion(package)), recorded$package_version)) {
      stop("Reference package version has changed; pass `genome` explicitly to confirm the reference.",
           call. = FALSE)
    }
    genome <- BSgenome::getBSgenome(package)
    return(.cfDNAPro_validate_reference(frag_obj_mut, genome))
  }
  # Older fragment RDS files have no package provenance. Infer only a known,
  # unambiguous assembly and still check its chromosome names and lengths.
  if (identical(md$cfDNAPro_genome_label, "user_define")) {
    stop("Custom reference metadata is incomplete; pass `genome` explicitly.", call. = FALSE)
  }
  labels <- unique(GenomeInfoDb::genome(frag_obj_mut))
  labels <- labels[!is.na(labels)]
  aliases <- c(hg19 = "hg19", hg38 = "hg38", GRCh38 = "hg38-NCBI",
               `hg38-NCBI` = "hg38-NCBI")
  if (length(labels) != 1L || !labels %in% names(aliases)) {
    stop("Cannot infer an unambiguous reference genome; pass `genome` explicitly.",
         call. = FALSE)
  }
  .cfDNAPro_validate_reference(frag_obj_mut,
                              .cfDNAPro_builtin_genome(unname(aliases[labels])))
}
