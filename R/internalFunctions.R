#' Internal functions
#'
#' @noRd
#' @importFrom stats median na.omit
#' @importFrom utils read.table
#' @importFrom rlang .data
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select rename group_by mutate summarise summarize
#' @importFrom dplyr lag bind_rows
#' @importFrom quantmod findValleys
#' @importFrom magrittr %>%



get_group_path <- function(path, groups) {
    # Check whether path input ends with file.seq.
    if (stringr::str_sub(path, start = -1) != .Platform$file.sep) {
        path <- paste0(path, .Platform$file.sep)
    }
    # Check whether the root path exists.
    if (dir.exists(path)) {
        group_path <- paste0(path, as.vector(groups))
        # Check whether the group path exists.
        if (dir.exists(group_path)) {
            return(group_path)
        } else {
            cat("Error: Group does not exist! ")
        }
    } else {
        cat("Error:Path does not exist!")
    }
}

get_full_sample_name <- function(x, input_type) {
    if (input_type == "bam") {
        full_sample_name <-
            list.files(x, full.names = TRUE, pattern = "\\.bam$", )
    } else if (input_type == "picard") {
        full_sample_name <-
            list.files(x, full.names = TRUE, pattern = "\\.txt$")
    } else if (input_type == "cfdnapro") {
        full_sample_name <-
            list.files(x, full.names = TRUE, pattern = "\\.txt$")
    } 

    # Change vector to a named list.
    named_vector <- stats::setNames(full_sample_name, full_sample_name)
    named_list <- as.list(named_vector)

    return(named_list)
}

read_picard_insert_metrics <- function(x) {
    insert_metrics <- read.table(x,
        header = TRUE,
        skip = 10,
        sep = "\t"
    )
    return(insert_metrics)
}

read_cfdnapro_insert_metrics <- function(x) {
  insert_metrics <- read.table(x,
                               header = TRUE,
                               skip = 0,
                               sep = "\t"
  )
  return(insert_metrics)
}

#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBam
#' @importFrom rlang .data

read_bam_insert_metrics_from_tlen <- function(file, ...) {
    flags <- Rsamtools::scanBamFlag(
        isPaired = TRUE,
        isProperPair = TRUE,
        isSecondMateRead = FALSE,
        isUnmappedQuery = FALSE,
        hasUnmappedMate = FALSE,
        isSecondaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE,
        isDuplicate = FALSE,
        isSupplementaryAlignment = FALSE,
    )

    # tags <- c("NM", "RG")

    # Match distance was stored in NM tag, and only select those <= 2.
    tag_filters <- list("NM" = c(0, 1, 2))

    what <- c("isize", "cigar")

    # We could not use simpleCigar to filter out those with INDELs,
    # Because we have soft clipping in the alignments. Try to use filter.

    params <- Rsamtools::ScanBamParam(
        what = what,
        flag = flags,
        # simpleCigar = TRUE,
        # tag = tags,
        tagFilter = tag_filters,
        mapqFilter = 30,
    )

    r1 <- Rsamtools::scanBam(
        file = file,
        param = params,
    )

    # Convert r1 to data.frame.

    r1 <- do.call(data.frame, r1[[1]])
    r1$isize <- abs(r1$isize)

    r2 <- r1 %>%
        dplyr::filter(.data$isize <= 1000 & !str_detect(
            string = .data$cigar,
            pattern = "[ID]"
        )) %>%
        dplyr::select(.data$isize) %>%
        dplyr::rename("insert_size" = .data$isize) %>%
        dplyr::group_by(.data$insert_size) %>%
        dplyr::summarise("All_Reads.fr_count" = dplyr::n())

    return(r2)
}

loop_read_insert_metrics_in_list <- function(x, input_type, ...) {

    # x should be a named list.
    if (missing(input_type)) {
        input_type <- "picard"
    }

    if (input_type == "picard") {
        d <- lapply(x, read_picard_insert_metrics, ...)
        
    } else if (input_type == "bam") {
        d <- lapply(x, read_bam_insert_metrics, ...)
        
    } else if (input_type == "cfdnapro") {
        d <- lapply(x, read_cfdnapro_insert_metrics, ...)
    }
    return(d)
}


calculate_sample_max_prop <- function(x) {
    # Transform to dataframe format.
    result <-
        as_tibble(x) %>%
        # Only keep those >=30 bp inserts.
        dplyr::filter(.data$insert_size >= 30) %>%
        # Calculate proportion and add new column called "prop".
    mutate(prop = .data$All_Reads.fr_count / sum(.data$All_Reads.fr_count)) %>%
        # Only keep the "insert_size" and "prop" column.
        select(c("insert_size", "prop")) %>%
        # Only keep the max prop.
        dplyr::filter(.data$prop == max(.data$prop))
    return(result)
}

calculate_peaks <- function(x) {
    peaks <-
        x %>%
        dplyr::mutate(id = row_number()) %>%
        dplyr::filter(.data$id %in% as.vector(quantmod::findPeaks(x$prop) - 1))
    return(peaks)
}

calculate_valleys <- function(x) {
    valleys <-
        x %>%
        dplyr::mutate(id = row_number()) %>%
        dplyr::filter(.data$id %in% as.vector(findValleys(x$prop) - 1))
    return(valleys)
}

filter_insert_size <- function(x, limit) {
    filtered <-
        x %>%
        as_tibble() %>%
        dplyr::filter(.data$insert_size >= as.vector(limit)[1] &
            .data$insert_size <= as.vector(limit)[2])
    return(filtered)
}

calcualte_interpeak_dist <- function(x) {
    dist <-
        x %>%
        mutate(interpeak_dist = .data$insert_size - lag(.data$insert_size))
    return(dist)
}

calcualte_intervalley_dist <- function(x) {
    dist <-
        x %>%
        mutate(intervalley_dist = .data$insert_size - lag(.data$insert_size))
    return(dist)
}


calculate_prop <- function(x) {
    result <-
        # Tranform to dataframe format.
        as_tibble(x) %>%
        # Only keep those >=30 bp inserts.
        dplyr::filter(.data$insert_size >= 30) %>%
        # Calcuate proportion and add new column called "prop".
    mutate(prop = .data$All_Reads.fr_count / sum(.data$All_Reads.fr_count)) %>%
        # Only keep the "insert_size" and "prop" column.
        dplyr::select(c("insert_size", "prop"))
    return(result)
}

calculate_prop_cdf <- function(x) {
    result <-
        # Tranform to dataframe format.
        as_tibble(x) %>%
        # Only keep those >=30 bp inserts.
        dplyr::filter(.data$insert_size >= 30) %>%
        # Calcuate proportion and add new column called "prop".
    mutate(prop = .data$All_Reads.fr_count / sum(.data$All_Reads.fr_count)) %>%
        # Add cdf column.
        dplyr::mutate(cdf = cumsum(.data$prop)) %>%
        # Add 1-cdf column.
        dplyr::mutate(one_minus_cdf = 1 - .data$cdf) %>%
        # Only keep those columns indicated.
        dplyr::select(c("insert_size", "prop", "cdf", "one_minus_cdf"))
    return(result)
}

# read the isize and freq data from bam/txt file ------------------------------

read_raw_data <- function(groups, path, input_type, ...) {
    # Convert groups into a named list.
    groups_list <- as.list(stats::setNames(groups, groups))

    result <-
        groups_list %>%
        # Set the path to  the folder containing .txt files.
        lapply(get_group_path, path = path) %>%
        # Get the full file name of each sample.
        lapply(get_full_sample_name, input_type = input_type) %>%
        # Read the insert size and frequency.
        lapply(loop_read_insert_metrics_in_list, input_type = input_type, ...)

    return(result)
}


# Calculate group prop, cdf, 1-cdf --------------------------------------------



AllGroupProp_List <- function(path, groups, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result1 <- read_raw_data(
        path = path,
        groups = groups,
        input_type = input_type,
        ...
    )
    result2 <- lapply(result1, function(x) {
        lapply(x, calculate_prop)
    })

    return(result2)
}

AllGroupProp_Dataframe <- function(path, groups, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result1 <- read_raw_data(
        path = path,
        groups = groups,
        input_type = input_type,
        ...
    )

    result2 <- lapply(result1, function(x) {
        lapply(x, calculate_prop)
    }) %>%
        lapply(bind_rows, .id = "file_name") %>%
        bind_rows(.id = "group")

    return(result2)
}

AllGroupPropCdf_List <- function(path, groups, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result1 <- read_raw_data(
        path = path,
        groups = groups,
        input_type = input_type,
        ...
    )

    result2 <- lapply(result1, function(x) {
        lapply(x, calculate_prop_cdf)
    })

    return(result2)
}

AllGroupPropCdf_Dataframe <- function(path,
                                      groups,
                                      input_type,
                                      ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result1 <- read_raw_data(
        path = path,
        groups = groups,
        input_type = input_type,
        ...
    )

    result2 <- lapply(
        result1,
        function(x) {
            lapply(x, calculate_prop_cdf)
        }
    ) %>%
        lapply(bind_rows, .id = "file_name") %>%
        bind_rows(.id = "group")

    return(result2)
}



# calculate mean and median -------------------------------------------------



calculate_mean_for_list <- function(x) {
    insert_size <- prop <- NULL

    r <-
        x %>%
        group_by(insert_size) %>%
        summarise(mean = mean(prop, na.rm = TRUE))
    return(r)
}

calculate_mean_for_df <- function(x) {
    insert_size <- prop <- group <- NULL
    r <-
        x %>%
        group_by(group, insert_size) %>%
        summarise(mean = mean(prop, na.rm = TRUE))
    return(r)
}

calculate_median_for_list <- function(x) {
    insert_size <- prop <- NULL
    r <-
        x %>%
        group_by(insert_size) %>%
        summarise(median = median(prop, na.rm = TRUE))
    return(r)
}

calculate_median_for_df <- function(x) {
    insert_size <- prop <- group <- NULL
    r <-
        x %>%
        group_by(group, insert_size) %>%
        summarise(median = median(prop, na.rm = TRUE))
    return(r)
}

calculate_propcdf_median_for_list <- function(x) {
    insert_size <- prop <- cdf <- one_minus_cdf <- NULL

    result <-
        x %>%
        group_by(insert_size) %>%
        summarise(
            "prop_median" = median(prop, na.rm = TRUE),
            "cdf_median" = median(cdf, na.rm = TRUE),
            "one_minus_cdf_median" = median(one_minus_cdf, na.rm = TRUE)
        )
    return(result)
}

calculate_propcdf_mean_for_list <- function(x) {
    insert_size <- prop <- cdf <- one_minus_cdf <- NULL

    result <-
        x %>%
        group_by(insert_size) %>%
        summarise(
            "prop_mean" = mean(prop, na.rm = TRUE),
            "cdf_mean" = mean(cdf, na.rm = TRUE),
            "one_minus_cdf_mean" = mean(one_minus_cdf, na.rm = TRUE)
        )
    return(result)
}

calculate_propcdf_meanmedian_for_list <- function(x) {
    insert_size <- prop <- cdf <- one_minus_cdf <- NULL

    result <-
        x %>%
        group_by(insert_size) %>%
        summarise(
            "prop_mean" = mean(prop, na.rm = TRUE),
            "prop_median" = median(prop, na.rm = TRUE),
            "cdf_mean" = mean(cdf, na.rm = TRUE),
            "cdf_median" = median(cdf, na.rm = TRUE),
            "one_minus_cdf_mean" = mean(one_minus_cdf, na.rm = TRUE),
            "one_minus_cdf_median" = median(one_minus_cdf, na.rm = TRUE)
        )
    return(result)
}


# calculate_breaks --------------------------------------------------------

calculate_breaks_for_mode_plot <- function(df, order) {
    # Define a vector "breaks".
    breaks <- vector()
    # Set the first value in breaks.
    breaks[1] <- nrow(dplyr::filter(df, .data$group == order[1])) / 2

    if (length(order) == 1) {
        return(breaks)
    } else if (length(order) > 1) {
        for (i in c(2:length(order))) {
            breaks[i] <-
                nrow(dplyr::filter(df, .data$group %in%
                    as.vector(order[seq_len(i - 1)]))) +
                nrow(dplyr::filter(df, .data$group == order[i])) / 2
        }
        return(breaks)
    }
}



# calculate_labels_for_mode_plot ------------------------------------------


calculate_labels_for_mode_plot <- function(df, order) {
    # Define a vector "labels".
    labels <- vector()

    for (i in c(seq_len(length(order)))) {
        labels[i] <- paste0(
            stringr::str_to_title(order[i]),
            "\n",
            "(",
            "n=",
            dplyr::filter(df, .data$group == order[i]) %>%
                group_by(.data$file_name) %>%
                count() %>%
                nrow(),
            ")"
        )
    }
    return(labels)
}



# check_dup_mode ----------------------------------------------------------

check_dup_mode <- function(x) {
    if ("TRUE" %in% x$duplication) {
        cat("ATTENTION! These samples have multiple modes: ")
        dup <- as.data.frame(dplyr::filter(x, x$duplication == "TRUE"))
        print(dup)
    }
}

# Fragment curation functions ----------------------------------------------

#' @importFrom GenomeInfoDb keepSeqlevels

bam_to_galp <- function(bamfile,
                        chromosome_to_keep = paste0("chr", 1:22),
                        strand_mode = 1,
                        genome_name = NULL,
                        ...) {
    # Check parameters
    
    if (!file.exists(bamfile)) {
        m <- paste0(bamfile, " doesn't exist, please check your input.")
        stop(m)
    }
    
    
    # Read bam into galp
    message("Reading bam into galp...")
    
    flag <- Rsamtools::scanBamFlag(
        isPaired = TRUE,
        isDuplicate = FALSE,
        isSecondaryAlignment = FALSE,
        isUnmappedQuery = FALSE,
        isSupplementaryAlignment = FALSE
    )
    
    what <- c("cigar", "mapq", "isize")
    
    
    param <- Rsamtools::ScanBamParam(flag = flag,
                                     what = what)
    
    # read bam file as galp object
    galp <- GenomicAlignments::readGAlignmentPairs(
        bamfile,
        index = bamfile,
        param = param,
        use.names = TRUE,
        strandMode = strand_mode
    )
    
    if(!is.null(genome_name)) {
      genome(galp) <- genome_name
    }
    
    # strandMode should be one for downstream operations
    stopifnot(GenomicAlignments::strandMode(galp) == 1)
    
    # remove read pairs without correct seqnames and strand information
    
    message("Curating seqnames and strand information...")
    # only keep chr1:22
    galp <-
        keepSeqlevels(galp, chromosome_to_keep, pruning.mode = "coarse")
    
    # remove read pairs without seqname information
    galp2 <- galp[!is.na(GenomicAlignments::seqnames(galp))]
    
    # remove read pairs without strand information
    galp3 <- galp2[GenomicAlignments::strand(galp2) != '*']
    
    
    return(galp3)
    
}







#' @importFrom BiocGenerics start end strand
curate_start_and_end <- function(galp)  {
    message("Correcting start and end coordinates of fragments...")
    get_start <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "+",
            start(GenomicAlignments::first(galp_object)),
            start(GenomicAlignments::second(galp_object))
        )
    }
    
    get_end <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "-",
            end(GenomicAlignments::first(galp_object)),
            end(GenomicAlignments::second(galp_object))
        )
    }
    
    # Add the strand and mutation locus information back to fragmentwise object
    fragmentwise <-  GRanges(
        seqnames = seqnames(GenomicAlignments::first(galp)),
        ranges = IRanges(start = get_start(galp),
                         end = get_end(galp)),
        strand = strand(galp),
        which_label.first = mcols(GenomicAlignments::first(galp))$which_label
    )
    
    return(fragmentwise)
}


#' @importFrom GenomeInfoDb isCircular
#' @importFrom BiocGenerics start end strand
get_out_of_bound_index <- function (x) {
    if (length(x) == 0L)
        return(integer(0))
    
    x_seqnames_id <- as.integer(seqnames(x))
    x_seqlengths <- unname(seqlengths(x))
    seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
    seqlength_is_na <- is.na(x_seqlengths)
    seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
    
    which(seqlevel_has_bounds[x_seqnames_id] &
              (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))
}

remove_out_of_bound_reads <- function(granges_object){
    message("Removing out-of-bound reads...")
    idx <-get_out_of_bound_index(granges_object)
    if(length(idx) != 0L)  granges_object <- granges_object[-idx]
    return(granges_object)
}


# motif functions

# Motif annotation

#' @importFrom GenomeInfoDb seqinfo seqnames seqlevels genome seqlengths
#' @importFrom BiocGenerics start end strand unstrand
#' @importFrom BSgenome getSeq

get_motif <- function(obj,
                      genome,
                      motif_type, 
                      motif_length){
  
  obj <- BiocGenerics::unstrand(obj)
  
  if(motif_type == "s") {
    coord <- GenomicRanges::GRanges(seqnames = seqnames(obj),
                     ranges = IRanges::IRanges(start = GenomicAlignments::start(obj) ,
                                      end = GenomicAlignments::start(obj) + 
                                        motif_length - 1))
  }
  
  
  if(motif_type == "e") {
    coord <- GenomicRanges::GRanges(seqnames = seqnames(obj),
                     ranges = IRanges::IRanges(start = GenomicAlignments::end(obj) - 
                                        motif_length + 1 ,
                                      end = GenomicAlignments::end(obj)))
  }
  
  if(motif_type == "u") {
    coord <- GenomicRanges::GRanges(seqnames = seqnames(obj),
                     ranges = IRanges::IRanges(start = GenomicAlignments::start(obj) - 
                                        motif_length,
                                      end = GenomicAlignments::start(obj) - 1))
  }
  
  if(motif_type == "d") {
    coord <- GenomicRanges::GRanges(seqnames = seqnames(obj),
                     ranges = IRanges::IRanges(start = GenomicAlignments::end(obj) + 1,
                                      end = GenomicAlignments::end(obj) + 
                                        motif_length))
  }
  
  seqlengths(coord) <- seqlengths(genome)[levels(seqnames(coord))]
  genome(coord) <- seqinfo(genome)@genome %>% unique()
  
  id <- get_out_of_bound_index(coord)
  
  if(length(id) != 0)  {
    coord <- plyranges::mutate(coord, index = 1:length(coord)) 
    coord_tidy <- coord[-id]
    coord_tidy <- plyranges::mutate(coord_tidy, 
                                    anno  = getSeq(genome, 
                                                   coord_tidy, 
                                                   as.character = TRUE))
    
    df <- tibble::tibble(index = coord_tidy$index, anno = coord_tidy$anno)
    result <- dplyr::bind_rows(df, tibble(index = id)) %>% 
      dplyr::arrange(index) %>%
      dplyr::pull(anno)
    
  } else {
    result <- BSgenome::getSeq(genome, coord, as.character = TRUE)
  }
  
  return(result)
  
}


# rle id 
rleid_cfdnapro <- function(x) {
  1L + cumsum(dplyr::coalesce(x != dplyr::lag(x), FALSE))
}





remove_outward_facing_readpairs <- function(galp) {
    message("Removing outward facing fragments...")
    # Get FR read pairs(i.e. first read aligned to positive strand)
    
    FR_id <-
        BiocGenerics::which(strand(GenomicAlignments::first(galp)) == "+")
    FR_galp <- galp[FR_id]
    
    # Get RF read pairs(i.e. second read aligned to positive strand)
    RF_id <-
        BiocGenerics::which(strand(GenomicAlignments::first(galp)) == "-")
    RF_galp <- galp[RF_id]
    
    # Get FR inward facing pairs
    FR_galp_inward <-
        FR_galp[GenomicAlignments::start(GenomicAlignments::first(FR_galp)) <
                    GenomicAlignments::end(GenomicAlignments::second(FR_galp))]
    
    # Get RF inward facing pairs
    RF_galp_inward <-
        RF_galp[GenomicAlignments::start(GenomicAlignments::second(RF_galp)) <
                    GenomicAlignments::end(GenomicAlignments::first(RF_galp))]
    
    # Combine RF and FR inward read pairs
    galp_corrected <- c(FR_galp_inward, RF_galp_inward)
    
    return(galp_corrected)
    
}

