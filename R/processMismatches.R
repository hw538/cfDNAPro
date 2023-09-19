processMismatches <- function(which_processed,
                              bed_full,
                              bamfile,
                              galp_flag,
                              galp_mapqFilter) {
  ############################################################################
  # Read the bam read by read to obtain mismatch information
  ############################################################################
  what <- c("seq", "mapq")

  param <- ScanBamParam(which = which_processed,
           flag = galp_flag,
           what = what,
           mapqFilter = galp_mapqFilter,
           tag = c("NM", "MD"))

  galp_modified_raw <- readGAlignments(file = bamfile,
                              use.names = TRUE,
                              with.which_label = TRUE,
                              param = param)

  message("BAM has been read ")

  ############################################################################
  # Convert GRanges to DF and remove indels
  ############################################################################
  gr_df <- as.data.frame(galp_modified_raw, use.outer.mcols = TRUE)

  gr_df <- gr_df[!grepl("I", gr_df$cigar), ]

  gr_df <- gr_df[!grepl("D", gr_df$cigar), ]

  ############################################################################
  # Clip the read sequences to remove soft clipped bases to match the MD tags
  ############################################################################
  gr_df_clipped <- clipSeq(gr_df)

  message("ClipSeq Finished")

  gr_df_clipped$read_id <- rownames(gr_df_clipped)

  gr_df_clipped$read_id  <- gsub("\\.\\d+", "", gr_df_clipped$read_id)

  gr_df_clipped <- gr_df_clipped %>%
    mutate(read_id = paste(read_id, strand, sep = ""))

  gr_df_clipped$read_id <- make.unique(gr_df_clipped$read_id, sep = ".")

  gr_df_clipped$read_id <- gsub("[+-]", "", gr_df_clipped$read_id)

  ############################################################################
  # Tidy up the clipped read dataframe and remove duplicate .1 annotation
  ############################################################################
  gr_df_clipped_reads <-  gr_df_clipped %>% select(
                          read_id,
                          seqnames,
                          start,
                          strand,
                          NM,
                          MD,
                          cigar,
                          clipseq,
                          which_label)

  # when fragments overlap multiple mutations, there will >= 2 duplicates
  gr_df_clipped_reads$read_id <- gsub("\\.1",
                                     ".1",
                                     as.character(gr_df_clipped_reads$read_id)
                                     )

  gr_df_clipped_reads$read_id <- gsub("\\.[2-3]",
                                     ".2",
                                     as.character(gr_df_clipped_reads$read_id)
                                     )

  gr_df_clipped_reads$read_id <- gsub("\\.[4-5]",
                                     ".3",
                                     as.character(gr_df_clipped_reads$read_id)
                                     )

  gr_df_clipped_reads$read_id <- gsub("\\.[6-7]",
                                     ".4",
                                     as.character(gr_df_clipped_reads$read_id)
                                     )

  colnames(gr_df_clipped_reads) <- c("read_id",
                                     "chr_start",
                                     "start_pos",
                                     "read_strand",
                                     "nm_tag",
                                     "md_tag",
                                     "cigar_string",
                                     "clip_seq",
                                     "which_locus"
  )

  ############################################################################
  # Annotated each read with mismatch information and derive a dataframe
  ############################################################################
  suppressWarnings(
    output_f <- purrr::pmap(gr_df_clipped_reads, mismatchPosition)
  )
  message("mismatchPosition finished")

  output_f_df <- as.data.frame(stri_list2matrix(output_f, byrow = TRUE))

  ############################################################################
  # Remove read name. They will be added back
  ############################################################################
  output_f_df_clean <- as.data.frame(lapply(
                        output_f_df, function(x) gsub(".*\\_", "\\1", x)
                        ))

  ############################################################################
  # Remove unnecessary locus positions from mismatches
  ############################################################################
  output_f_df_clean <- as.data.frame(lapply(
                        output_f_df_clean, function(x) gsub("\\:c.*", "\\1", x)
                          ))

  ############################################################################
  # Add back read's locus of interest as a separate column
  ############################################################################
  output_f_df_clean$which_locus <- gsub(".*\\:c", "\\1c", output_f_df$V1)

  ############################################################################
  # Add back read ID with strand information (+/-) to a separate column
  ############################################################################
  output_f_df_clean$read_id <- gsub("\\_.*", "\\1", output_f_df$V1)

  # TODO: parallel processing and chunking
  # MAY need to change to mclapply

  # get just the locus
  output_f_df_clean$locus <- gsub("-\\d+", "", output_f_df_clean$which_locus)

  # Define the V column names
  v_columns <- grep("^V\\d+$", names(output_f_df_clean), value = TRUE)

  # Define the number of rows per chunk
  chunk_size <- 10000

  # Calculate the number of chunks needed
  num_chunks <- ceiling(nrow(output_f_df_clean) / chunk_size)

  # Split the dataframe into chunks
  df_chunks <- split(
    output_f_df_clean, rep(
      1:num_chunks, each = chunk_size, length.out = nrow(output_f_df_clean)))

  # Create a function to process each chunk
  process_chunk <- function(chunk) {
    # Iterate over each row of the chunk
    for (i in 1:nrow(chunk)) {
      # Get the 'which_locus' value for the current row
      target_locus <- chunk$locus[i]

      # Check if V1 already contains the string "read.position=0"
      if (!grepl("read.position=0", chunk$V1[i], fixed = TRUE)) {
        # Check which V columns partially match with the 'which_locus' value
        matching_columns <- v_columns[
          sapply(v_columns, function(column) {
            grepl(target_locus, chunk[i, column], fixed = TRUE)
          })
        ]
        # If matching columns are found
        if (length(matching_columns) > 0) {
          # Move the first matching V column to V1
          chunk$V1[i] <- chunk[i, matching_columns[1]]

          # Set all other V columns (excluding V1) to NA
          chunk[i, v_columns[-1]] <- NA
        } else {
          # If no matching column is found,
          # set V1 as the 'which_locus' value plus "::read.position=0"
          chunk$V1[i] <- paste0(target_locus, "::read.position=0")

          # Set all V columns (excluding V1) to NA
          chunk[i, v_columns[-1]] <- NA
        }
      }
    }

    # Return the processed chunk
    return(chunk)
  }

  # Apply the function to each chunk in parallel
  # TODO: core number should be provided by the USER
  cl <- makeCluster(4)  # Create a cluster with 4 cores

  clusterExport(cl, c("v_columns"), envir = environment())  # Export the required objects to cluster

  processed_chunks <- parLapply(cl, df_chunks, process_chunk)

  stopCluster(cl)  # Stop the cluster

  # Merge the processed chunks back into a single dataframe
  output_f_df_clean <- do.call(rbind, processed_chunks)

  ############################################################################
  # Convert wide format into long format. This gets rid of the NAs.
  ############################################################################
  output_f_df_clean_gathered <- output_f_df_clean %>%
                                  gather(key = "key",
                                         value = "mismatch",
                                         starts_with("V"),
                                         na.rm = TRUE)

  ############################################################################
  # Tidy up the dataframe by selecting columns and removing unnecessary info
  ############################################################################
  output_f_df_clean_gathered <- output_f_df_clean_gathered %>%
                                    select(which_locus,
                                           read_id,
                                           mismatch)

  output_f_df_clean_gathered$which_locus <- gsub(
                                    "\\-.*",
                                    "\\1",
                                    output_f_df_clean_gathered$which_locus)

  ############################################################################
  # read.position=0 indicates that there there were no mismatches in a read
  ############################################################################

  ############################################################################
  # Thus, reads with read.position=0 read get assigned REF as their base
  ############################################################################

  ############################################################################
  # It indicates the read overlapped the locus of interest but had REF as base
  ############################################################################
  output_f_df_clean_gathered$allele <- ifelse(
                                        grepl(
                                        "read.position=0",
                                        output_f_df_clean_gathered$mismatch),
                                        "REF",
                                        "ALT")

  ############################################################################
  # Tidy up the dataframe by splitting the mismatch column into extra columns
  ############################################################################

  ############################################################################
  # Mismatch position within the read placed in its own column
  # Mismatch positions are in the 5'-> 3' direction regardless of read strand
  ############################################################################
  output_f_df_clean_gathered$read_position <- gsub(
                                        ".*\\:r",
                                        "\\1r",
                                        output_f_df_clean_gathered$mismatch)

  ############################################################################
  # Removing mismatch position information from the mismatch column
  ############################################################################
  output_f_df_clean_gathered$mismatch <- gsub(
                                        "\\:read.position.*",
                                        "\\1",
                                        output_f_df_clean_gathered$mismatch)

  ############################################################################
  # Listing SNV as REF if there were no mismatches found in the read
  ############################################################################
  output_f_df_clean_gathered$snv <- ifelse(grepl(
            "read.position=0",
            output_f_df_clean_gathered$read_position),
            paste0(output_f_df_clean_gathered$which_locus, sep = ":", "REF"),
            output_f_df_clean_gathered$mismatch)

  ############################################################################
  # Using the Mutation file to query the mismatch dataframe
  ############################################################################
  bed_full_df <- bed_full %>% select(chr, start, alt)

  ############################################################################
  # Obtain a vector of targeted mutations
  ############################################################################
  bed_vector <- apply(bed_full_df, 1, function(x) paste(x, collapse = ":"))

  bed_vector <- gsub(" ", "", bed_vector)

  ############################################################################
  # Query whether each mismatch in every read matches the targeted mutation
  ############################################################################
  output_f_df_clean_gathered$match_mutation_bed <- ifelse(
                        output_f_df_clean_gathered$snv %in% bed_vector,
                        TRUE,
                        FALSE)

  ############################################################################
  # Assign fragment ID to each mismatch
  ############################################################################
  output_f_df_clean_gathered$fragment_id <- gsub(
                                    ":([^:]+)$",
                                    "",
                                    output_f_df_clean_gathered$read_id)

  ############################################################################
  # Place mismatch genomic positions into a separate column
  ############################################################################
  output_f_df_clean_gathered$snv_pos <- gsub(
                ":([^:]+)$",
                "",
                output_f_df_clean_gathered$snv)

  ############################################################################
  # Place read strand information into a separate column
  ############################################################################
  output_f_df_clean_gathered$strand <- gsub(
                    ".*\\:",
                    "",
                    output_f_df_clean_gathered$read_id)

  ############################################################################
  # Add SNV+Strand column to indicate which paired-end read had the mismatch
  ############################################################################
  output_f_df_clean_gathered$snv_read_specific <- paste0(
                    output_f_df_clean_gathered$snv,
                    sep = ":", output_f_df_clean_gathered$strand)

  ############################################################################
  # Obtain mutation of interest positions as a vector
  ############################################################################
  bed_full_df_pos <- bed_full %>% select(chr, start)

  bed_vector_pos <- apply(
                        bed_full_df_pos,
                        1,
                        function(x) paste(x, collapse = ":"))

  bed_vector_pos <- gsub(" ", "", bed_vector_pos)

  ############################################################################
  # Focus on the mutation file positions to extract relevant information
  ############################################################################
  `%nin%` <- Negate(`%in%`)

  ############################################################################
  # Reads with mismatch positions absent in mutation file, get assigned REF
  ############################################################################
  output_f_df_clean_gathered$snv_read_specific <- ifelse(
    output_f_df_clean_gathered$snv_pos %nin% bed_vector_pos,
    paste0(output_f_df_clean_gathered$which_locus,
           sep = ":REF:",
           output_f_df_clean_gathered$strand),
    paste0(output_f_df_clean_gathered$snv,
           sep = ":",
           output_f_df_clean_gathered$strand))

  ############################################################################
  # New column to keep reads that match the positions of mutation file
  ############################################################################
  output_f_df_clean_gathered$match_mutation_bed_pos <- ifelse(
        output_f_df_clean_gathered$snv_pos %in% bed_vector_pos,
        TRUE,
        FALSE)

  ############################################################################
  # Dataframe tidy up
  ############################################################################
  final_output <- output_f_df_clean_gathered %>% select(
                      fragment_id,
                      snv_read_specific,
                      match_mutation_bed,
                      which_locus)

  final_output$snv_read_specific <- paste(final_output$snv_read_specific,
                                          sep = "|",
                                          final_output$match_mutation_bed)

  ############################################################################
  # Wide format to convert DF from read-level to fragment-level
  ############################################################################
  # Add a row number column to handle duplicates
  df <- final_output[, -3]

  # Convert fragment_id to character type
  df <- df %>%
    mutate(fragment_id = as.character(fragment_id))

  df$new_id <- paste(
    sub("\\..*$", "", df$fragment_id), ":", df$which_locus, sep = "")

  df_wide <- df %>%
          group_by(new_id) %>%
          mutate(row_num = row_number()) %>%
          ungroup() %>%
          pivot_wider(
            id_cols = new_id,
            names_from = row_num,
            values_from = snv_read_specific,
            names_prefix = "snv_read_specific_"
          ) %>%
          unnest(c(starts_with("snv_read_specific_")))

  df_wide$fragment_id <- gsub(":chr[0-9]+:[0-9]+", "", df_wide$new_id)

  df_wide$fragment_id <- make.unique(df_wide$fragment_id, sep = ".")

  final_output_wide <- df_wide %>% select(
    fragment_id,
    snv_read_specific_1,
    snv_read_specific_2
    )

  ############################################################################
  # Add fragment level information the mismatch and paired end read support
  ############################################################################
  final_output_wide$concordance <- ifelse(
    is.na(final_output_wide$snv_read_specific_1) |
    is.na(final_output_wide$snv_read_specific_2),
    NA,
    str_extract(
        final_output_wide$snv_read_specific_1, ".*:") == str_extract(
            final_output_wide$snv_read_specific_2, ".*:")
  )

  final_output_wide$bed_locus <- ifelse(
    !is.na(final_output_wide$snv_read_specific_1),
    gsub("(:[^:]+){2}$", "", final_output_wide$snv_read_specific_1),
    ifelse(
        !is.na(final_output_wide$snv_read_specific_2),
        gsub("(:[^:]+){2}$", "", final_output_wide$snv_read_specific_2),
    NA)
  )

  ############################################################################
  # Discern between: base = MUT file, base = REF, discordant R1/R2, other base
  ############################################################################
  final_output_wide$match  <- ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_1) &
        is.na(final_output_wide$snv_read_specific_2),
        gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_1),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_2) &
        is.na(final_output_wide$snv_read_specific_1),
        gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_2),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_1) &
        is.na(final_output_wide$snv_read_specific_2),
        paste0("REF"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_2) &
        is.na(final_output_wide$snv_read_specific_1),
        paste0("REF"),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_1) &
        grepl("REF", final_output_wide$snv_read_specific_2),
        paste0("discordant"),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_2) &
        grepl("REF", final_output_wide$snv_read_specific_1),
        paste0("discordant"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_1) &
        grepl("^(?!.*REF).*FALSE.*$",
              final_output_wide$snv_read_specific_2, perl = TRUE),
        paste0("discordant"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_2) &
        grepl("^(?!.*REF).*FALSE.*$",
              final_output_wide$snv_read_specific_1, perl = TRUE),
        paste0("discordant"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_1) &
        grepl("REF", final_output_wide$snv_read_specific_2),
        paste0("REF"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_2) &
        grepl("REF", final_output_wide$snv_read_specific_1),
        paste0("REF"),
    ifelse(
        grepl("^(?!.*REF).*FALSE.*$",
            final_output_wide$snv_read_specific_1, perl = TRUE) &
        is.na(final_output_wide$snv_read_specific_2),
        paste0(gsub(":([^:]+)$", "",
               final_output_wide$snv_read_specific_1), sep = ":not_mut"),
    ifelse(
        grepl("^(?!.*REF).*FALSE.*$",
              final_output_wide$snv_read_specific_2, perl = TRUE) &
        is.na(final_output_wide$snv_read_specific_1),
        paste0(gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_2),
               sep = ":not_mut"),
    ifelse(
        grepl("^(?!.*REF).*FALSE.*$",
            final_output_wide$snv_read_specific_1, perl = TRUE) &
        grepl("^(?!.*REF).*FALSE.*$",
            final_output_wide$snv_read_specific_2, perl = TRUE),
        paste0(gsub(":([^:]+)$", "",
               final_output_wide$snv_read_specific_1), sep = ":not_mut"),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_1) &
        grepl("TRUE", final_output_wide$snv_read_specific_2),
        gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_1),
  NA))))))))))))))

  ############################################################################
  # Can wrap into separate function for each long if else - check notes (Phil)
  ############################################################################

  ############################################################################
  # Summarise fragment by R1/R2 mismatch and concordance
  ############################################################################
  final_output_wide$mismatch_status  <- ifelse(
        grepl("TRUE", final_output_wide$concordance) &
        !grepl("REF", final_output_wide$match),
        paste0(final_output_wide$match, sep = ":concordant"),
    ifelse(
        is.na(final_output_wide$concordance) &
        !grepl("REF", final_output_wide$match),
        paste0(final_output_wide$match, sep = ":single_read"),
    ifelse(
        is.na(final_output_wide$concordance) &
        grepl("REF", final_output_wide$match),
        paste0(final_output_wide$bed_locus, sep = ":REF:single_read"),
    ifelse(
        !is.na(final_output_wide$concordance) &
        grepl("REF", final_output_wide$match),
        paste0(final_output_wide$bed_locus, sep = ":REF:concordant"),
    paste0(final_output_wide$bed_locus, sep = ":discordant")))))

  ############################################################################
  # If fragment has only one read mutated from the pair, assign strand info
  ############################################################################
  final_output_wide$single_read_strand <- ifelse(
    grepl("single_read", final_output_wide$mismatch_status) &
    !is.na(final_output_wide$snv_read_specific_1),
    str_extract(final_output_wide$snv_read_specific_1,
                "(?<=:)[^:|]*(?=\\|[^|]*$)"),
    ifelse(
    grepl("single_read", final_output_wide$mismatch_status) &
    !is.na(final_output_wide$snv_read_specific_2),
    str_extract(final_output_wide$snv_read_specific_2,
                "(?<=:)[^:|]*(?=\\|[^|]*$)"),
    ""))

  ############################################################################
  # Tidy up Dataframe
  ############################################################################
  final_output_wide$fragment_mismatch_bed <- paste0(
    final_output_wide$mismatch_status,
    final_output_wide$single_read_strand)

  final_df <- final_output_wide  %>% select(fragment_id, fragment_mismatch_bed)

  colnames(final_df) <- c("fragment_id", "locus_info")

  final_df$locus <- ifelse(!grepl("not_mut", final_df$locus_info),
                          gsub("(:[^:]+){2}$", "", final_df$locus_info),
                          ifelse(grepl("not_mut", final_df$locus_info),
                          gsub("(:[^:]+){3}$", "", final_df$locus_info), NA))

  ############################################################################
  # Final tidy up and slight renaming
  ############################################################################
  final_df$locus_status <- ifelse(
        !grepl("REF", final_df$locus_info) &
        !grepl("not_mut", final_df$locus_info),
        paste0("MUT", sep = ":", gsub(".*:([^:]+)$", "\\1",
               final_df$locus_info)),
    ifelse(
        grepl("REF", final_df$locus_info),
        paste0("REF", sep = ":",
               gsub(".*:([^:]+)$", "\\1", final_df$locus_info)),
    ifelse(
        grepl("not_mut", final_df$locus_info),
        paste0("other_base", sep = ":",
               gsub(".*:([^:]+)$", "\\1", final_df$locus_info)),
    paste0("discordant"))))

  ############################################################################
  # Remove redundant strand information
  ############################################################################
  final_df$locus_status <- gsub("[+|-]", "", final_df$locus_status)

  #############################################################################
  # process fragment mismatch dataframe for merging
  #############################################################################
  mismatch_df <- final_df  %>% select(locus_info,
                                      locus_status,
                                      fragment_id)

  rownames(mismatch_df) <- final_df$fragment_id

  return(mismatch_df)
}
