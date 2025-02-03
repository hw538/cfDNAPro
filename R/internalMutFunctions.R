
#' Internal Functions for mutational processing
#' @importFrom magrittr '%>%'
#' @noRd

# Helper function to process a set of columns with a given suffix
#' @importFrom gsubfn strapply
#' @importFrom dplyr mutate filter select
#' @noRd
process_columns <- function(df, suffix) {
  # Setup extra columns with null values
  df[[paste0("softclip", suffix)]] <- rep(0, nrow(df))
  df[[paste0("left_softclip", suffix)]] <- rep(0, nrow(df))
  df[[paste0("right_softclip", suffix)]] <- rep(0, nrow(df))
  df[[paste0("clipseq", suffix)]] <- as.vector(df[[paste0("seq", suffix)]])
  df[[paste0("clipqual", suffix)]] <- as.vector(df[[paste0("qual", suffix)]])

  # Remove hard clipping information from the cigar string
  df[[paste0("cigar_processed", suffix)]] <- gsub(
    "\\d+H", "", df[[paste0("cigar", suffix)]])

  # Get the subset of rows with soft-clipping
  softclip_rows <- grep("S", df[[paste0("cigar", suffix)]])
  softclip_df <- df[softclip_rows, ]

  # Calculate the number of bases soft clipped from the left and right
  softclip_df[[paste0("left_softclip", suffix)]] <- apply(
    softclip_df, 1, function(x) {
      sum(as.numeric(gsubfn::strapply(
        as.character(x[[paste0("cigar", suffix)]]), "^(\\d+)S",
        simplify = c)))
  })
  softclip_df[[paste0("right_softclip", suffix)]] <- apply(
    softclip_df, 1, function(x) {
      sum(as.numeric(gsubfn::strapply(
        as.character(x[[paste0("cigar", suffix)]]), "(\\d+)S$",
        simplify = c)))
  })

  # Extract the sequence and quality scores
  softclip_df[[paste0("clipseq", suffix)]] <- substr(
    softclip_df[[paste0("seq", suffix)]],
    softclip_df[[paste0("left_softclip", suffix)]] + 1,
    nchar(softclip_df[[paste0("seq", suffix)]]) -
    softclip_df[[paste0("right_softclip", suffix)]])
  softclip_df[[paste0("clipqual", suffix)]] <- substr(
    softclip_df[[paste0("qual", suffix)]],
    softclip_df[[paste0("left_softclip", suffix)]] + 1,
    nchar(softclip_df[[paste0("qual", suffix)]]) -
    softclip_df[[paste0("right_softclip", suffix)]])

  # Substitute in the values
  df[[paste0("left_softclip", suffix)]][softclip_rows] <-
    softclip_df[[paste0("left_softclip", suffix)]]
  df[[paste0("right_softclip", suffix)]][softclip_rows] <-
    softclip_df[[paste0("right_softclip", suffix)]]
  df[[paste0("seq", suffix)]][softclip_rows] <-
    softclip_df[[paste0("clipseq", suffix)]]
  df[[paste0("qual", suffix)]][softclip_rows] <-
    softclip_df[[paste0("clipqual", suffix)]]

  # Calculate the length of the clipped sequence
  df[[paste0("clipped_len", suffix)]] <-
    nchar(df[[paste0("clipseq", suffix)]])

  return(df)
}

#' Function to remove soft clipped bases to match the MD tags
#' #' @param reads_df A dataframe containing read pairs with 'cigar.first' and
#' 'cigar.last' columns that may include soft clipping indicators.
#' @noRd
clip_read_pair <- function(reads_df = NULL) {

  # Check and process columns with soft clipping in 'cigar.first'
  if (any(grepl("S", reads_df$cigar.first, fixed = TRUE))) {
    reads_df <- process_columns(reads_df, ".first")
  }

  # Check and process columns with soft clipping in 'cigar.last'
  if (any(grepl("S", reads_df$cigar.last, fixed = TRUE))) {
    reads_df <- process_columns(reads_df, ".last")
  }

  return(reads_df)
}


#' Function to check Mutation file columns
#' @param df A dataframe representing the mutation file.
#' @noRd
check_mutfile_columns <- function(df) {
  # Define expected columns
  expected_cols <- c("chr", "pos", "ref", "alt")
  expected_types <- c("character", "integer", "character", "character")

  # Check if all expected columns are present
  if (!all(expected_cols %in% colnames(df))) {
      return(paste("Mutation file column names are incorrect or missing.",
                   "Expected columns:", toString(expected_cols), "."))
  }

  # Check if the column types are correct
  actual_types <- sapply(df[expected_cols], class)
  if (!all(actual_types == expected_types)) {
    return(paste("Mutation File column types are incorrect, found: ",
                 toString(actual_types), " expected: ",
                 toString(expected_types)))
  }

  # If all checks pass, return confirmation
  return("Mutation file checks passed successfully!")
}

#' Remove closely clustered mutations
#'
#' This function sorts a mutation data frame by chromosome and position,
#' and removes mutations that are less than 1000 base pairs apart to ensure
#' no mutations are closer than this threshold. It avoids duplicate read pairs
#' in the GRanges object, crucial for accurate genomic analysis. The threshold
#' was chosen based on the findings on kataegis, described as having an average
#' intermutation distance of ≤1000 bp in Alexandrov LB et al., Nature, 2013.
#'
#' @param mutation_data A data frame with mutation information.
#' @return A data frame with mutations spaced more than 1000 base pairs apart.
#' @importFrom dplyr arrange
#' @examples
#' # Assuming df has columns 'chr' and 'pos'
#' cleaned_data <- remove_clustered_mutations(df)
#' @noRd
remove_clustered_mutations <- function(mutation_data) {
  # Sort by chromosome and position
  sorted_data <- mutation_data %>%
    arrange(chr, pos)

  # Initialize a vector to keep track of indices to keep
  keep_indices <- c()

  # Loop through the sorted data frame
  last_pos <- -Inf
  last_chr <- ""

  for (i in seq_len(nrow(sorted_data))) {
    current_chr <- sorted_data$chr[i]
    current_pos <- sorted_data$pos[i]

    # Check if current mut is on the same chr and close to the last one
    if (current_chr != last_chr || current_pos - last_pos >= 1000) {
      keep_indices <- c(keep_indices, i)
      last_pos <- current_pos
      last_chr <- current_chr
    }
  }

  # Subset the data frame to keep only the non-clustered mutations
  non_clustered_data <- sorted_data[keep_indices, ]
  return(non_clustered_data)
}

#' Function to read in the .tsv mutation file
#' @param mutation_file A string path to the TSV file containing mutation data.
#'        The file should contain columns named 'chr', 'pos', 'ref', and 'alt'.
#' @importFrom utils read.table
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @noRd
read_mutation_file <- function(mutation_file) {
  # Check if the file exists
  if (!file.exists(mutation_file)) {
    stop("Mutation file does not exist.")
  }

  # Read in the TSV mutation file
  message("Reading in the provided mutation file...")
  bed <- read.table(mutation_file,
                    header = TRUE, sep="\t",
                    quote="", stringsAsFactors = FALSE)

  # Check whether the file format is correct
  result <- check_mutfile_columns(bed)

  if (result == "Mutation file checks passed successfully!") {
    message("Mutation file is formatted correctly.")
  
    # Remove mutations that are less than 1000 bp apart
    bed <- remove_clustered_mutations(bed)

    return(head(bed))
  } else {
    stop(result)
  }
}





#' Function to parse CIGAR strings and extract insertions
#' @noRd
insertion_pos <- function(
    unique_read_pair_id,
    chr,
    start_pos,
    end_pos,
    read_strand,
    nm_tag,
    md_tag,
    cigar_string,
    clip_seq,
    clip_qual,
    mapq,
    which_locus
) {
    # Initialize position trackers
    position_in_read <- 1
    position_in_ref <- start_pos

    # Array to store insertions
    insertions <- list()

    # Parse CIGAR string for operations
    elements <- strsplit(cigar_string,
                         "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",
                         perl = TRUE)[[1]]
    numbers <- as.integer(grep("\\d+", elements, value = TRUE))
    letters <- grep("[MIDNSHPX=]", elements, value = TRUE)

    # Loop through elements of the CIGAR string
    for (i in seq_along(letters)) {
        len <- numbers[i]
        op <- letters[i]

        if (op %in% c("M", "=", "X")) {
            # Advance position counters for match or sequence match/mismatch
            position_in_read <- position_in_read + len
            position_in_ref <- position_in_ref + len
        } else if (op == "I") {
            # Handle insertion
            bases <- substring(clip_seq,
                               position_in_read,
                               position_in_read + len - 1)
            insertions[[length(insertions) + 1]] <- list(
                PositionInRef = position_in_ref,
                PositionInRead = position_in_read,
                Length = len,
                Bases = bases
            )
            position_in_read <- position_in_read + len
        } else if (op == "D") {
            # Advance reference and read position for deletion
            # Deletion characters (D) will be added later to the seq columns
            position_in_read <- position_in_read + len
            position_in_ref <- position_in_ref + len
        } else if (op == "S" || op == "H") {
            # Handle soft clipping and hard clipping
            position_in_read <- position_in_read + len
        }
    }

    # Construct genomic coordinates for insertions
    insertion_coordinates <- sapply(insertions, function(ins) {
        paste(
              paste0("r_str:", read_strand),
              paste0("g_pos:", chr, ":", ins$PositionInRef, "-",
                     ins$PositionInRef - 1 + ins$Length, ":", "*", ins$Bases),
              paste0("r_pos:", ins$PositionInRead, "-",
                     ins$PositionInRead - 1 + ins$Length),
              "ins",
              sep="|")
    })

    # Return a list containing the insertion coordinates
    return(list(Insertions = insertion_coordinates))
}

#' Function to extract insertion information from paired reads
#' @noRd
apply_to_paired_reads_ins <- function(data_frame) {
    # Initialize a list to hold results
    results <- list()

    # Iterate over each row of the DataFrame
    for (i in seq_len(nrow(data_frame))) {
        row <- data_frame[i, ]

        # Extract information for the first part of the read pair
        first_results <- insertion_pos(
            unique_read_pair_id = row$unique_read_pair_id_first,
            chr = row$chr_first,
            start_pos = row$start_pos_first,
            end_pos = row$end_pos_first,
            read_strand = row$read_strand_first,
            nm_tag = row$nm_tag_first,
            md_tag = row$md_tag_first,
            cigar_string = row$cigar_string_first,
            clip_seq = row$clip_seq_first,
            clip_qual = row$clip_qual_first,
            mapq = row$mapq_first,
            which_locus = row$which_locus
        )

        # Extract information for the last part of the read pair
        last_results <- insertion_pos(
            unique_read_pair_id = row$unique_read_pair_id_last,
            chr = row$chr_last,
            start_pos = row$start_pos_last,
            end_pos = row$end_pos_last,
            read_strand = row$read_strand_last,
            nm_tag = row$nm_tag_last,
            md_tag = row$md_tag_last,
            cigar_string = row$cigar_string_last,
            clip_seq = row$clip_seq_last,
            clip_qual = row$clip_qual_last,
            mapq = row$mapq_last,
            which_locus = row$which_locus
        )

        # Combine the results from first and last into one list for each row
        results[[i]] <- list(First = first_results, Last = last_results)
    }

    # Return the list of results
    return(results)
}

#' Function to remove insertion bases from the seq columns
#' @importFrom stringr str_split str_detect
#' @noRd
remove_insertion_bases <- function(cigar, clip_seq) {
  elements <- strsplit(cigar, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",
                       perl=TRUE)[[1]]
  numbers <- as.integer(grep("\\d+", elements, value = TRUE))
  letters <- grep("[MDISHPX=]", elements, value = TRUE)

  adjusted_seq <- clip_seq
  position_in_read <- 1

  # Iterate over CIGAR elements
  for (i in seq_along(letters)) {
    len <- numbers[i]
    op <- letters[i]

    if (op == "I") {
      # Calculate the actual position of the insertion in the adjusted sequence
      insertion_start <- position_in_read
      insertion_end <- position_in_read + len - 1

      # Remove insertion bases from the sequence
      if (insertion_start <= nchar(adjusted_seq)) {
        part_before_insertion <- substring(adjusted_seq, 1,
                                           insertion_start - 1)
        part_after_insertion <- substring(adjusted_seq,
                                          insertion_end + 1)
        adjusted_seq <- paste0(part_before_insertion,
                               part_after_insertion)
      }
      position_in_read <- position_in_read + len
    } else if (op %in% c("M", "X", "=")) {
      # Advance the position in the read by the length of M, DEL 
      position_in_read <- position_in_read + len
    } else if (op == "S" || op == "H" || op == "D") {
      # Skip soft clips and hard clips in sequence position calculations
      if (i == 1 || i == length(letters)) {
        next
      }
    }
  }

  return(adjusted_seq)
}

#' Function to remove insertion quality scores from the base quality columns
#' @noRd
remove_insertion_qual_bases <- function(cigar, clip_qual) {
  elements <- strsplit(cigar, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",
                       perl=TRUE)[[1]]
  numbers <- as.integer(grep("\\d+", elements, value = TRUE))
  letters <- grep("[MDISHPX=]", elements, value = TRUE)

  adjusted_qual <- clip_qual
  position_in_read <- 1

  # Iterate over CIGAR elements
  for (i in seq_along(letters)) {
    len <- numbers[i]
    op <- letters[i]

    if (op == "I") {
      # Calculate the actual position of the insertion in the quality string
      insertion_start <- position_in_read
      insertion_end <- position_in_read + len - 1

      # Remove insertion quality bases from the string
      if (insertion_start <= nchar(adjusted_qual)) {
        part_before_insertion <- substring(adjusted_qual, 1,
                                           insertion_start - 1)
        part_after_insertion <- substring(adjusted_qual,
                                          insertion_end + 1)
        adjusted_qual <- paste0(part_before_insertion,
                                part_after_insertion)
      }
      position_in_read <- position_in_read + len
    } else if (op %in% c("M", "X", "=")) {
      position_in_read <- position_in_read + len
    } else if (op == "S" || op == "H" || op == "D") {
      if (i == 1 || i == length(letters)) {
        next
      }
    }
  }

  return(adjusted_qual)
}

#' Function to process insertion information
#' @importFrom dplyr mutate filter
#' @importFrom magrittr %>%
#' @noRd
process_insertions <- function(reads_df) {
  if (any(grepl("I", reads_df$cigar_string_first)) ||
      any(grepl("I", reads_df$cigar_string_last))) {

    # Filter for rows with insertions in the CIGAR string
    ins_reads_df <- reads_df[grepl("I", reads_df$cigar_string_first) |
                             grepl("I", reads_df$cigar_string_last), ]

    # Process insertions
    insertion_results <- apply_to_paired_reads_ins(ins_reads_df)

    # Alternative way to display mismatches, deletions, insertions
    ins_reads_df$MDI_first <- ""
    ins_reads_df$MDI_last <- ""

    # Iterate over the results to populate the new columns
    for (i in seq_along(insertion_results)) {
        # Check and assign insertions from First part of the pair
        if (!is.null(insertion_results[[i]]$First) &&
            length(insertion_results[[i]]$First$Insertions) > 0) {
            # Convert the list of insertions to a character string, if necessary
            ins_reads_df$MDI_first[i] <- paste(
              unlist(insertion_results[[i]]$First$Insertions), collapse = "; ")
        } else {
            # Handle cases where there are no insertions
            ins_reads_df$MDI_first[i] <- ""
        }

        # Check and assign insertions from Last part of the pair
        if (!is.null(insertion_results[[i]]$Last) &&
            length(insertion_results[[i]]$Last$Insertions) > 0) {
            ins_reads_df$MDI_last[i] <- paste(
              unlist(insertion_results[[i]]$Last$Insertions), collapse = "; ")
        } else {
            ins_reads_df$MDI_last[i] <- ""
        }
    }

    # Process each row of the dataframe and remove insertion bases read-pairs
    ins_reads_df$clip_seq_first <- mapply(remove_insertion_bases,
                                          ins_reads_df$cigar_string_first,
                                          ins_reads_df$clip_seq_first)
    ins_reads_df$clip_seq_last <- mapply(remove_insertion_bases,
                                         ins_reads_df$cigar_string_last,
                                         ins_reads_df$clip_seq_last)

    # Remove insertion base qualities
    ins_reads_df$clip_qual_first <- mapply(remove_insertion_qual_bases,
                                           ins_reads_df$cigar_string_first,
                                           ins_reads_df$clip_qual_first)
    ins_reads_df$clip_qual_last <- mapply(remove_insertion_qual_bases,
                                          ins_reads_df$cigar_string_last,
                                          ins_reads_df$clip_qual_last)

    # Combine the ins and mm_del information
    # Function to check if there are no insertions in the CIGAR strings
    has_no_insertions <- function(cigar_str) {
      !grepl("I", cigar_str)
    }

    # Filter dataframe to include only rows with no insertions in CIGAR strings
    mm_del_reads_df <- reads_df[
      has_no_insertions(reads_df$cigar_string_first) &
      has_no_insertions(reads_df$cigar_string_last),
    ]

    # Identifying missing columns in mm_del dataframe
    missing_in_mm_del <- base::setdiff(
      names(ins_reads_df), names(mm_del_reads_df))

    # Adding missing columns as NA to mm_del dataframe
    mm_del_reads_df[missing_in_mm_del] <- lapply(missing_in_mm_del,
                                                 function(x) "")

    # Ensure columns are in the same order
    ins_reads_df <- ins_reads_df[names(mm_del_reads_df)]

    # Combine mm_del and ins dataframes
    reads_df <- rbind(ins_reads_df, mm_del_reads_df)

  } else {

  # If no insertions, add MDI_first and MDI_last columns with empty values
  reads_df <- reads_df %>%
    dplyr::mutate(MDI_first = "", MDI_last = "")
  }

  return(reads_df)

}




#' Function to add 'D' for bases that were deleted in the seq columns
#' @noRd
adjust_sequence_for_deletions <- function(md_tag, clip_seq) {
    elements <- strsplit(md_tag, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",
                         perl=TRUE)[[1]]
    adjusted_seq <- clip_seq
    position_in_read <- 1  # Tracks the current position in the read

    for (element in elements) {
        if (grepl("^\\d+$", element)) {
            # Advance position by the number of matches
            position_in_read <- position_in_read + as.numeric(element)
        } else if (grepl("^\\^[ACGTNacgtn]+$", element)) {
            # Deletion from the reference; extract the bases
            deletion_length <- nchar(substring(element, 2, nchar(element)))
            # Insert 'D' at current pos for the length of the deletion
            adjusted_seq <- paste0(substring(adjusted_seq,
                                             1,
                                             position_in_read - 1),
                                   paste(rep("D", deletion_length),
                                         collapse = ""),
                                   substring(adjusted_seq,
                                             position_in_read))
        } else if (grepl("^[ACGTNacgtn]+$", element)) {
            # Handle mismatches: advance by the length of the mismatch element
            position_in_read <- position_in_read + nchar(element)
        }
    }

    return(adjusted_seq)
}

#' Function to add '!' for bases that were deleted in the qual columns
#' @noRd
adjust_quality_for_deletions <- function(md_tag, clip_qual) {
    elements <- strsplit(md_tag, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",
                         perl=TRUE)[[1]]
    adjusted_qual <- clip_qual
    position_in_read <- 1

    for (element in elements) {
        if (grepl("^\\d+$", element)) {
            # Advance position by the number of matches
            position_in_read <- position_in_read + as.numeric(element)
        } else if (grepl("^\\^[ACGTNacgtn]+$", element)) {
            # Deletion from the reference; extract the bases
            deletion_length <- nchar(substring(element, 2, nchar(element)))
            # Insert '!' at current pos in quality string for the length of del
            adjusted_qual <- paste0(substring(adjusted_qual,
                                              1,
                                              position_in_read - 1),
                                    paste(rep("!", deletion_length),
                                          collapse = ""),
                                    substring(adjusted_qual,
                                              position_in_read))
        } else if (grepl("^[ACGTNacgtn]+$", element)) {
            # Handle mismatches: advance by the length of the mismatch element
            position_in_read <- position_in_read + nchar(element)
        }
    }

    return(adjusted_qual)
}




#' Function to add mismatch and deletion info as separate columns for readpairs
#' @noRd
mm_del_pos <- function(
    unique_read_pair_id,
    chr,
    start_pos,
    end_pos,
    read_strand,
    nm_tag,
    md_tag,
    cigar_string,
    clip_seq,
    clip_qual,
    mapq,
    which_locus
    ) {

    # Split the MD tag into components (numbers and letters)
    md_elements <- regmatches(md_tag, gregexpr("([0-9]+|\\^[ACGT]+|[ACGT]+)",
                              md_tag, perl = TRUE))[[1]]

    # Initialize position trackers
    position_in_read <- 1
    position_in_ref <- start_pos + 1

    # Arrays to store mismatches and deletions
    mismatches <- list()
    deletions <- list()

    for (element in md_elements) {
        if (grepl("^[0-9]+$", element)) {
            # Number of bases matched (no discrepancy)
            num_matches <- as.integer(element)
            position_in_read <- position_in_read + num_matches
            position_in_ref <- position_in_ref + num_matches
        } else if (grepl("^\\^[ACGT]+", element)) {
            # Deletion from the reference
            deletion_length <- nchar(substr(element, 2, nchar(element)))
            deletions[[length(deletions) + 1]] <- list(
                PositionInRef = position_in_ref,
                PositionInRead = position_in_read,
                Length = deletion_length,
                Bases = substring(element, 2, nchar(element))
            )
            position_in_read <- position_in_read + deletion_length
            position_in_ref <- position_in_ref + deletion_length
        } else if (grepl("^[ACGT]+$", element)) {
            # Mismatches
            for (base_index in seq(nchar(element))) {
                base_from_clip_seq <- substring(clip_seq, position_in_read,
                                                position_in_read)
                mismatches[[length(mismatches) + 1]] <- list(
                    PositionInRef = position_in_ref,
                    PositionInRead = position_in_read,
                    Base = substring(element, base_index, base_index),
                    ClipSeqBase = base_from_clip_seq
                )
                position_in_read <- position_in_read + 1
                position_in_ref <- position_in_ref + 1
            }
        }
    }

    # Construct genomic coordinates for mismatches
    mismatch_coordinates <- sapply(mismatches, function(m) {
        paste(
              paste0("r_str:", read_strand),
              paste0("g_pos:", chr, ":", m$PositionInRef - 1, "-",
                     m$PositionInRef - 1, ":", m$Base, "-", m$ClipSeqBase),
              paste0("r_pos:", m$PositionInRead),
              "mm",
              sep="|")
    })

    # Construct genomic coordinates for deletions
    deletion_coordinates <- sapply(deletions, function(d) {
        paste(
              paste0("r_str:", read_strand),
              paste0("g_pos:", chr, ":", d$PositionInRef - 1, "-",
                     d$PositionInRef - 2 + d$Length, ":", "^",d$Bases),
              paste0("r_pos:", d$PositionInRead, "-",
                     d$PositionInRead - 1 + d$Length),
              "del",
              sep="|")
    })

    # Return list containing both mismatches and del with detailed coordinates
    list(
        Mismatches = mismatch_coordinates,
        Deletions = deletion_coordinates
    )
}

#' Wrapper function applying mm_del_pos to each row with _first/_last suffix
#' @noRd
apply_to_paired_reads <- function(data_frame) {
    results <- lapply(seq_len(nrow(data_frame)), function(i) {
        row <- data_frame[i, ]
        first_results <- mm_del_pos(row$unique_read_pair_id,
                                    row$chr_first,
                                    row$start_pos_first,
                                    row$end_pos_first,
                                    row$read_strand_first,
                                    row$nm_tag_first,
                                    row$md_tag_first,
                                    row$cigar_string_first,
                                    row$clip_seq_first,
                                    row$clip_qual_first,
                                    row$mapq_first,
                                    row$which_locus
                                    )

        last_results <- mm_del_pos(row$unique_read_pair_id,
                                   row$chr_last,
                                   row$start_pos_last,
                                   row$end_pos_last,
                                   row$read_strand_last,
                                   row$nm_tag_last,
                                   row$md_tag_last,
                                   row$cigar_string_last,
                                   row$clip_seq_last,
                                   row$clip_qual_last,
                                   row$mapq_last,
                                   row$which_locus
                                   )

        list(First = first_results, Last = last_results)
    })
    return(results)
}

#' Function to obtain updated mismatch and deletion information
#' @noRd
update_mdi_columns <- function(reads_df) {
  results <- apply_to_paired_reads(reads_df)

  # Iterate over the results to update the MDI_first and MDI_last columns
  for (i in seq_along(results)) {
    # Update MDI_first column with mismatches and deletions
      if (!is.null(results[[i]]$First) &&
          !is.null(results[[i]]$First$Mismatches)) {
      # Create a string from the list of mismatches
          new_first_mismatches <- paste(unlist(results[[i]]$First$Mismatches),
                                        collapse = "; ")
          # Append this new data to the existing data in the dataframe
          if (nzchar(reads_df$MDI_first[i])) {
              reads_df$MDI_first[i] <- paste(reads_df$MDI_first[i],
                                             new_first_mismatches,
                                             sep = "; ")
          } else {
            reads_df$MDI_first[i] <- new_first_mismatches
          }
      }

      if (!is.null(results[[i]]$First) &&
          !is.null(results[[i]]$First$Deletions)) {
          # Create a string from the list of deletions
          new_first_deletions <- paste(unlist(results[[i]]$First$Deletions),
                                       collapse = "; ")
          # Append or set data as done with mismatches
          if (nzchar(reads_df$MDI_first[i])) {
              reads_df$MDI_first[i] <- paste(reads_df$MDI_first[i],
                                             new_first_deletions, sep = "; ")
          } else {
              reads_df$MDI_first[i] <- new_first_deletions
          }
      }

      # Update MDI_last column with mismatches and deletions
      if (!is.null(results[[i]]$Last) &&
          !is.null(results[[i]]$Last$Mismatches)) {
          # Create a string from the list of mismatches
          new_last_mismatches <- paste(unlist(results[[i]]$Last$Mismatches),
                                       collapse = "; ")
          # Append this new data to the existing data in the dataframe
          if (nzchar(reads_df$MDI_last[i])) {
              reads_df$MDI_last[i] <- paste(reads_df$MDI_last[i],
                                            new_last_mismatches,
                                            sep = "; ")
          } else {
              reads_df$MDI_last[i] <- new_last_mismatches
          }
      }

    if (!is.null(results[[i]]$Last) && !is.null(results[[i]]$Last$Deletions)) {
      # Create a string from the list of deletions
          new_last_deletions <- paste(unlist(results[[i]]$Last$Deletions),
                                      collapse = "; ")
          # Append or set data as done with mismatches
          if (nzchar(reads_df$MDI_last[i])) {
              reads_df$MDI_last[i] <- paste(reads_df$MDI_last[i],
                                            new_last_deletions,
                                            sep = "; ")
      } else {
        reads_df$MDI_last[i] <- new_last_deletions
      }
    }
  }
  return(reads_df)
}




#' Function to add the base information to mutation_status
#' @importFrom dplyr mutate case_when select
#' @importFrom stringr str_extract str_detect
#' @noRd
add_base_info <- function(reads_df) {
  reads_df <- reads_df %>%
    mutate(
      # Extract substrings from MDI columns based on mutation_location
      MDI_first_details = str_extract(
        MDI_first, paste0(mutation_location, "\\:[^|]+")),
      MDI_last_details = str_extract(
        MDI_last, paste0(mutation_location, "\\:[^|]+")),

      # Create new column with reference base(s) before '-' in which_mutation
      ref_base = gsub("-.*", "", gsub(".*\\:", "", which_mutation)),

      # For discordant cases, concatenate details from MDI first/last columns
      mutation_status = case_when(
        # When discordant
        str_detect(mutation_status, "discordant") ~ paste0(
          mutation_status,
          ":",
          gsub(".*\\:", "", MDI_first_details),
          "/",
          gsub(".*\\:", "", MDI_last_details)),

        # When concordant readpair and mutations are detected as described
        str_detect(mutation_status,
                   "MUT:read_pair_concordant|OTHER:read_pair_concordant") ~
                   paste0(mutation_status, ":",
                          gsub("^.*\\:", "", MDI_first_details), "/",
                          gsub("^.*\\:", "", MDI_last_details)),

        # When there's a single read positive or negative with a mutation
        (str_detect(mutation_status, "MUT:single_read") &
                    mutation_location_present_last) ~
                    paste0(mutation_status, ":",
                           gsub("^.*\\:", "", MDI_last_details)),
        (str_detect(mutation_status, "MUT:single_read") &
                    mutation_location_present_first) ~
                    paste0(mutation_status, ":",
                           gsub("^.*\\:", "", MDI_first_details)),

        # For single-read REF, add the extracted base to the mutation_status
        str_detect(mutation_status, "REF:single_read") ~
                   paste0(mutation_status, ":", ref_base),

        # For read-pair REF, add the extracted base twice separated by '/'
        str_detect(mutation_status, "REF:read_pair_concordant") ~
                   paste0(mutation_status, ":", ref_base, "/", ref_base),

        # Default to existing mutation_status if none of the conditions apply
        TRUE ~ mutation_status
      )
    ) %>%
    select(-c(mutation_present_first,
              mutation_present_last,
              mutation_location_present_last,
              mutation_location_present_first,
              MDI_first_details, MDI_last_details))
}

#' Function to add mutation overlap status to a dataframe
#' @importFrom dplyr mutate case_when select
#' @importFrom stringr str_extract
#' @noRd
add_mutation_overlap_status <- function(df) {
  # Update the dataframe with new mutation overlap status
  df_updated <- df %>%
    mutate(
      # Extract mutation positions from the which_mutation column
      start_mutation = as.numeric(str_extract(which_mutation,
                                  "(?<=:)[0-9]+(?=-)")),
      end_mutation = as.numeric(str_extract(which_mutation,
                                "(?<=-)[0-9]+")),

      # Determine if there is any overlap with the mutation positions
      mutation_overlap = case_when(
        # Check for overlap between both reads and the mutation range
        (start_pos_first <= end_mutation & end_pos_first >= start_mutation) &
        (start_pos_last <= end_mutation & end_pos_last >= start_mutation) ~
        "read_pair",

        # Check if only the first read overlaps the mutation range
        (start_pos_first <= end_mutation & end_pos_first >= start_mutation) ~
        paste0("single_read", tolower(read_strand_first)),

        # Check if only the last read overlaps the mutation range
        (start_pos_last <= end_mutation & end_pos_last >= start_mutation) ~
        paste0("single_read", tolower(read_strand_last)),

        # Default case where there is no overlap
        TRUE ~ "no_overlap"
      )
    ) %>%
    select(-start_mutation, -end_mutation)

  # Return the updated dataframe
  return(df_updated)
}

#' Function to derive information about the mutation locus from each read-pair
#' @importFrom dplyr mutate case_when
#' @importFrom stringr str_extract str_detect
#' @noRd
process_mutation_status <- function(reads_df) {
  reads_df <- reads_df %>%
    mutate(
      # Extract just the location part of the mutation (before the "|")
      mutation_location = str_extract(which_mutation, "^[^:]+:[^:]+"),

      # Check presence of mutation location in MDI columns
      mutation_location_present_first = str_detect(
        MDI_first, mutation_location),
      mutation_location_present_last = str_detect(
        MDI_last, mutation_location),

      # Check presence of full mutation (including bases) in MDI columns
      mutation_present_first = str_detect(MDI_first, which_mutation),
      mutation_present_last = str_detect(MDI_last, which_mutation),

      # Create mutation_status based on conditions described
      mutation_status = case_when(
        # When there's a read pair overlap and no detections
        mutation_overlap == "read_pair" &
        !mutation_location_present_first & !mutation_location_present_last ~
        "REF:read_pair_concordant",

        # When read-pair overlaps and mutations are detected as described
        mutation_overlap == "read_pair" &
        mutation_present_first & mutation_present_last ~
        "MUT:read_pair_concordant",

        # When read-pair overlaps and only one read has the mutation
        mutation_overlap == "read_pair" &
        mutation_present_first != mutation_present_last ~
        "MUT:read_pair_discordant",

        # When both MDI columns match the location but not the full mutation
        mutation_overlap == "read_pair" &
        mutation_location_present_first & mutation_location_present_last &
        !mutation_present_first & !mutation_present_last ~
        "OTHER:read_pair_concordant",

        # When read-pair overlaps and only one of MDI column matches location
        mutation_overlap == "read_pair" &
        mutation_location_present_first != mutation_location_present_last &
        !mutation_present_first & !mutation_present_last ~
        "OTHER:read_pair_discordant",

        # For single_read+/- overlap where the read displays another mutation
        mutation_overlap %in% c("single_read+", "single_read-") &
        (mutation_location_present_first != mutation_location_present_last) &
        !mutation_present_first & !mutation_present_last ~
        paste("OTHER:", mutation_overlap, sep=""),

        # For single_read+/- where the mutation matches exactly as described
        mutation_overlap %in% c("single_read+", "single_read-") &
        mutation_present_first != mutation_present_last ~
        paste("MUT:", mutation_overlap, sep=""),

        # For single_read+/- where there are no alterations at the mut location
        mutation_overlap %in% c("single_read+", "single_read-") &
        mutation_present_first == mutation_present_last ~
        paste("REF:", mutation_overlap, sep=""),

        # Default case if none of the above conditions are met
        TRUE ~ "not_applicable"
      )
    )
  return(reads_df)
}

#' Funcion to process base qualities
#' @importFrom dplyr mutate select
#' @importFrom stringr str_extract
#' @noRd
process_base_qualities <- function(df) {
  df %>%
    mutate(
      # Extract the numerical position part of the mutation_location
      mutation_position = str_extract(
        mutation_location, "(?<=:)[0-9]+-[0-9]+"),

      # Extract start and end positions from mutation_position
      mut_start = as.integer(str_extract(mutation_position, "^[0-9]+")),
      mut_end = as.integer(str_extract(mutation_position, "[0-9]+$")),

      # Determine overlap with the mutation location
      overlap_first = (start_pos_first <= mut_end &
                       end_pos_first >= mut_start),
      overlap_last = (start_pos_last <= mut_end &
                      end_pos_last >= mut_start),

      # Calculate read position for the mutation within the reads
      read_pos_first = mut_start - start_pos_first + 1,
      read_pos_last = mut_start - start_pos_last + 1,

      # Extract the base quality at the mutation position from quality strings
      bq_first = ifelse(overlap_first, substr(clip_qual_first,
                                              read_pos_first,
                                              read_pos_first), NA),
      bq_last = ifelse(overlap_last, substr(clip_qual_last,
                                            read_pos_last,
                                            read_pos_last), NA),

      # Convert Phred+33 quality scores to integer scores
      bq_first_int = ifelse(!is.na(bq_first) & nchar(bq_first) == 1,
                            sapply(
                              bq_first, function(x) as.integer(
                                charToRaw(x)) - 33),
                            NA),
      bq_last_int = ifelse(!is.na(bq_last) & nchar(bq_last) == 1,
                          sapply(
                            bq_last, function(x) as.integer(
                              charToRaw(x)) - 33),
                          NA),

      # Calculate mean quality score when both reads overlap the mut location
      bq_first_int = as.numeric(bq_first_int),
      bq_last_int = as.numeric(bq_last_int),
      mutation_locus_bq = ifelse(
        overlap_first & overlap_last,
        rowMeans(cbind(bq_first_int, bq_last_int), na.rm = TRUE),
        ifelse(overlap_first, bq_first_int,
              ifelse(overlap_last, bq_last_int, NA)))
      ) %>%
      select(-c(bq_first, bq_last,
                read_pos_first, read_pos_last,
                bq_first_int, bq_last_int,
                overlap_first, overlap_last,
                mut_start, mut_end))

}




# Function to calculate read length from CIGAR string, considering only 'M'
#' @importFrom stringr str_extract_all
#' @noRd
calculate_read_length <- function(cigar) {
  # Return 0 immediately if cigar is NA
  if (is.na(cigar)) {
    return(0)
  }

  # Extract numbers and their corresponding operations
  numbers <- as.integer(unlist(str_extract_all(cigar, "\\d+")))
  letters <- unlist(str_extract_all(cigar, "[A-Z]"))

  # Identify indices where the operation is 'M'
  m_indices <- letters == "M"

  # Calculate the total length contributing to the read's footprint
  # Sum only lengths where the operation is 'M'
  # Return 0 if there are no 'M' operations
  if (any(m_indices)) {
    length <- sum(numbers[m_indices])
  } else {
    length <- 0
  }
  
  return(length)
}




#' Function to process duplicate read-pair IDs in galp/granges objects
#' and generate a dataframe for later processing
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicAlignments GAlignmentPairs
#' @importFrom dplyr mutate
#' @importFrom stringr str_detect
#' @noRd
process_duplicate_read_pairs <- function(genomic_data, mut_fragments_only) {
  # Initialize the reads_df variable
  reads_df <- NULL
  
  # Check if input is a GAlignments or GRanges object and adjust accordingly
  if (is(genomic_data, "GAlignmentPairs")) {
    reads_df <- as.data.frame(genomic_data)
    # Update which_label.first by removing '.' and characters that follow it
    reads_df <- reads_df %>%
      mutate(which_label.first = gsub("\\..*$", "", which_label.first))
    # Create unique identifier for GAlignmentPairs as before
    reads_df$unique_read_pair_id <- paste(gsub("\\.1$", "",
                                          rownames(reads_df)),
                                          reads_df$which_label.first,
                                          sep = "_")
  } else if (is(genomic_data, "GRanges")) {
    # Convert GRanges to dataframe without row names to avoid duplicates
    reads_df <- as.data.frame(genomic_data, row.names = NULL)
    # Use names of GRanges to create id
    reads_df$unique_read_pair_id <- paste(names(genomic_data))
    if(mut_fragments_only) {
      # Update which_label.first by removing '.' and following characters
      reads_df <- reads_df %>%
        mutate(which_label.first = gsub("\\..*$", "", which_label.first))
      # Use names of GRanges and which_label.first to create a unique id
      reads_df$unique_read_pair_id <- paste(names(genomic_data),
                                            reads_df$which_label.first,
                                            sep = "_")
    } else if (!mut_fragments_only) {
      # Update which_label.first by removing '.' and following characters
      reads_df <- reads_df %>%
        mutate(which_label.first = gsub("\\..*$", "", which_label.first))
      # Use names of GRanges and which_label.first to create a unique id
      reads_df$unique_read_pair_id <- paste(names(genomic_data),
                                            reads_df$which_label.first,
                                            sep = "_")
      reads_df <- reads_df %>%
        mutate(unique_read_pair_id = gsub("_NA$", "", unique_read_pair_id))
    }
  } else {
    stop("The input object is neither a GAlignmentPairs nor a GRanges object.")
  }

  # Remove any old row names
  rownames(reads_df) <- NULL

  # Return the modified dataframe
  return(reads_df)
}




#' Function to subset alignments to mutational positions
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom GenomicAlignments readGAlignmentPairs strandMode seqnames strand
#' @importFrom GenomeInfoDb merge seqinfo keepSeqlevels
#' @importFrom S4Vectors isSingleStringOrNA
#' @importFrom Rsamtools ScanBamParam
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @import Rsamtools
#' @noRd
bam_to_galp_mut <- function(bamfile,
                            use_names = TRUE,
                            param = galp_param,
                            chromosome_to_keep = FALSE,
                            strand_mode = 1,
                            genome = NA_character_) {
  # Check parameters
  stopifnot(file.exists(bamfile))
  stopifnot(isSingleStringOrNA(genome) || is(genome, "Seqinfo"))

  # Read bam into galp
  message("Reading bam into galp...")

  galp <- GenomicAlignments::readGAlignmentPairs(file = bamfile,
                            use.names = use_names,
                            strandMode = strand_mode,
                            param = param,
                            with.which_label = TRUE)

  # Message indicating successful read and possibly some details about 'galp'
  message("BAM file has been read into galp.")

  # add genome information
  if (isSingleStringOrNA(genome)) {
    genome <- Seqinfo(genome = genome)
  }
  seqinfo(galp) <- merge(GenomeInfoDb::seqinfo(galp), genome)

  # strandMode should be one for downstream operations
  stopifnot(GenomicAlignments::strandMode(galp) == 1)

  # only keep needed seqnames
  if (!isFALSE(chromosome_to_keep)) {
    galp <- keepSeqlevels(galp, chromosome_to_keep, pruning.mode = "coarse")

  }

  message("Curating seqnames and strand information...")
  # remove read pairs without correct seqnames and strand information
  galp2 <- galp[!is.na(GenomicAlignments::seqnames(galp))]
  galp3 <- galp2[GenomicAlignments::strand(galp2) != "*"]

  return(galp3)

}

#' Function to generate GRanges for specific positions
#' @importFrom IRanges IRanges
#' @noRd
make_granges <- function(loci) {
  loci_gr <- GenomicRanges::GRanges(seqnames = loci[, 1],
                     ranges = IRanges::IRanges(start = loci[, 2],
                                      end = loci[, 2]))
  return(loci_gr)
}

#' Function to process and curate mutational fragment-level information
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @import Rsamtools
#' @importFrom methods is
#' @importFrom stringr str_detect str_extract
#' @importFrom dplyr mutate select filter left_join rename
#' @importFrom tibble as_tibble
#' @importFrom tidyr replace_na
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqinfo
#' @noRd
process_mutation_fragments <- function(
            bamfile = bamfile,
            genome = genome,
            galp = galp,
            mutation_file = mutation_file,
            frag = frag,
            galp_bqFilter = galp_bqFilter,
            chromosome_to_keep = chromosome_to_keep,
            mut_fragments_only = mut_fragments_only) {

  # check if works
  if (!is.null(mutation_file) && !mut_fragments_only) {
      # Extract the count from the metadata
      metadata_string <- galp@metadata[[1]]
      count <- as.numeric(sub(".*count:", "", metadata_string))

      # Check if the count is valid and within the appropriate range
      if (!is.na(count) && count > 0 && count <= length(galp)) {
          # Subset galp to only include the range from 1 to the extracted count
          galp <- galp[1:count]
      } else {
          message("Count is out of range or not specified correctly.")
      }
  }

  # Process mutation list data
  loci_df <- read_mutation_file(mutation_file = mutation_file)
  loci_df <- base::subset(loci_df, chr %in% chromosome_to_keep)

  # Process duplicate read-pair names
  reads_df <- process_duplicate_read_pairs(
    genomic_data = galp, mut_fragments_only = mut_fragments_only)

  # Identify columns that end with ".1"
  columns_to_remove <- grep("\\.1$", names(reads_df), value = TRUE)

  # Exclude these columns from the data frame
  reads_df <- reads_df[, !names(reads_df) %in% columns_to_remove]

  # Check if there is soft clipping in 'cigar.first' or 'cigar.last'
  if (any(grepl("S", reads_df$cigar.first)) || 
      any(grepl("S", reads_df$cigar.last))) {
        
        # Call the function if soft clipping is detected
        reads_df <- clip_read_pair(reads_df)
  }

  # Select relevant columns from the 'first' and 'last' alignment pairs
  reads_df <- reads_df %>% dplyr::select(unique_read_pair_id,
                                        seqnames.first,
                                        start.first,
                                        end.first,
                                        strand.first,
                                        NM.first,
                                        MD.first,
                                        cigar.first,
                                        seq.first,
                                        qual.first,
                                        mapq.first,
                                        which_label.first,
                                        seqnames.last,
                                        start.last,
                                        end.last,
                                        strand.last,
                                        NM.last,
                                        MD.last,
                                        cigar.last,
                                        seq.last,
                                        qual.last,
                                        mapq.last,
                                        which_label.last)

  # Rename columns to a more uniform and readable format
  colnames(reads_df) <- c("unique_read_pair_id",
                          "chr_first",
                          "start_pos_first",
                          "end_pos_first",
                          "read_strand_first",
                          "nm_tag_first",
                          "md_tag_first",
                          "cigar_string_first",
                          "clip_seq_first",
                          "clip_qual_first",
                          "mapq_first",
                          "which_locus_first",
                          "chr_last",
                          "start_pos_last",
                          "end_pos_last",
                          "read_strand_last",
                          "nm_tag_last",
                          "md_tag_last",
                          "cigar_string_last",
                          "clip_seq_last",
                          "clip_qual_last",
                          "mapq_last",
                          "which_locus_last")

  # Use loci_df to annotate df with mutational information
  loci_df$which_locus <- with(loci_df, paste0(chr, ":", pos, "-", pos))

  loci_df$which_mutation <- with(loci_df, {
    locus <- paste0(chr, ":", pos, "-", pos)
    mutation <- paste0(ref, "-", alt)
    paste0(locus, ":", mutation)
  })

  # Left join to add which_mutation
  reads_df <- reads_df %>%
    left_join(loci_df %>% select(which_locus, which_mutation),
              by = c("which_locus_first" = "which_locus")) %>%
    dplyr::rename("which_locus" = "which_locus_first") %>%
    select(-which_locus_last)

  # Process insertions
  reads_df <- process_insertions(reads_df)

  # Process deletions
  # First, check if any rows contain 'D' in the CIGAR strings
  if (any(grepl("D", reads_df$cigar_string_first)) ||
      any(grepl("D", reads_df$cigar_string_last))) {
    # Applying adjustments only to rows that have deletions in their MD tags
    reads_df <- reads_df %>%
      mutate(
        clip_seq_first = ifelse(grepl("\\^", md_tag_first),
                                        mapply(adjust_sequence_for_deletions,
                                               md_tag_first, clip_seq_first),
                                        clip_seq_first),
        clip_seq_last = ifelse(grepl("\\^", md_tag_last),
                                        mapply(adjust_sequence_for_deletions,
                                               md_tag_last, clip_seq_last),
                                        clip_seq_last),
        clip_qual_first = ifelse(grepl("\\^", md_tag_first),
                                          mapply(adjust_quality_for_deletions,
                                                 md_tag_first, clip_qual_first),
                                          clip_qual_first),
        clip_qual_last = ifelse(grepl("\\^", md_tag_last),
                                        mapply(adjust_quality_for_deletions,
                                               md_tag_last, clip_qual_last),
                                        clip_qual_last)
      )
  }

  # Update MDI columns
  reads_df <- update_mdi_columns(reads_df)

  # Add read overlap information to read-pairs
  reads_df <- add_mutation_overlap_status(reads_df)

  # Derive information about the mutation locus from each read-pair
  reads_df <- process_mutation_status(reads_df)

  # Add base information
  reads_df <- add_base_info(reads_df)

  # Add the mutation locus basequalities
  reads_df <- process_base_qualities(reads_df)

  # Filter the dataframe
  reads_df <- reads_df %>%
    filter(mutation_locus_bq >= galp_bqFilter)

  # Select only the specified columns
  reads_df <- reads_df %>%
    select(unique_read_pair_id, which_mutation,
           mutation_status, mutation_locus_bq)

  # Extract Seqinfo from the existing frag GRanges object to add it back later
  seqinfo_frag <- GenomeInfoDb::seqinfo(frag)

  # Generate the dataframe containing fragment length data
  frag <- process_duplicate_read_pairs(
    genomic_data = frag, mut_fragments_only = mut_fragments_only)

  # Subset the fragment length dataframe
  frag <- frag %>%
    dplyr::select(unique_read_pair_id, seqnames, start, end, width, strand)

  # Merge dataframes
  frag <- suppressWarnings(merge(
    frag,
    reads_df,
    by = "unique_read_pair_id",
    all.x = TRUE,
    all.y = FALSE
    ))

  # Remove rows with NA in the 'chr' column
  frag <- frag %>%
    dplyr::filter(!is.na(seqnames))

  # Rename target column to target_mutation
  names(frag)[names(frag) == "which_mutation"] <- "target_mutation"

  # Convert df to GRanges:
  frag <- makeGRangesFromDataFrame(frag,
                                   keep.extra.columns = TRUE,
                                   seqinfo = seqinfo_frag)

  gr_meta <- suppressWarnings(mcols(frag))

  gr_meta$target_mutation[is.na(gr_meta$target_mutation)] <- "outer_fragment"
  gr_meta$mutation_status[is.na(gr_meta$mutation_status)] <- "outer_fragment"
  gr_meta$mutation_locus_bq[is.na(
    gr_meta$mutation_locus_bq)] <- "outer_fragment"
  mcols(frag) <- gr_meta

  frag@metadata <- list("GRanges object with fragment and mutational data")

  return(frag)

}

#' Function to process BAM files based on mutational or general alignment data
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @import Rsamtools
#' @importFrom GenomicAlignments last cigar
#' @importFrom dplyr filter
#' @importFrom GenomicRanges GRanges
#' @noRd
process_bam_file <- function(bamfile, 
                             mutation_file = NULL,
                             mut_fragments_only = FALSE,
                             use_names = FALSE,
                             chromosome_to_keep = NULL,
                             strand_mode = 1,
                             genome_name,
                             galp_flag,
                             galp_what,
                             galp_tag,
                             galp_mapqFilter) {

  # Define general and ScanBam parameters
  galp_param <- Rsamtools::ScanBamParam(
    flag = galp_flag,
    what = galp_what,
    tag = galp_tag,
    mapqFilter = galp_mapqFilter
  )

  # Check the presence of a mutation file and the mode of operation
  if (!is.null(mutation_file)) {
    loci_df <- suppressMessages(
      read_mutation_file(mutation_file = mutation_file))
    loci_df <- base::subset(loci_df, chr %in% chromosome_to_keep)

    if (nrow(loci_df) == 0) {
      message("No mutations for selected chromosomes in the mutation file.")
      return(create_empty_galp())
    }

    # Create genomic ranges from loci
    which_loci <- make_granges(loci = loci_df)
    # Define galp params for mutational processing
    galp_param_mut <- Rsamtools::ScanBamParam(
      flag = galp_flag,
      what = galp_what,
      tag = galp_tag,
      mapqFilter = galp_mapqFilter,
      which = which_loci
    )

    if (mut_fragments_only) {
      galp <- bam_to_galp_mut(bamfile = bamfile, 
                              use_names = use_names,
                              chromosome_to_keep = chromosome_to_keep,
                              strand_mode = strand_mode,
                              genome = genome_name,
                              param = galp_param_mut)
      # Update metadata with content specifics
      galp@metadata <- list(
        "GRanges with mutation-specific read-pairs only")

      if (length(galp) == 0) {
        message("No reads found in bam file for specified mutation loci.")
        return(create_empty_galp())
      }
    } else {
      # Fetch both general and mutation-specific alignments
      galp <- bam_to_galp2(bamfile = bamfile,
                           use_names = use_names,
                           chromosome_to_keep = chromosome_to_keep,
                           strand_mode = strand_mode,
                           genome = genome_name,
                           param = galp_param)

      galp_mm <- suppressMessages(bam_to_galp_mut(
                                  bamfile = bamfile,
                                  use_names = use_names,
                                  chromosome_to_keep = chromosome_to_keep,
                                  strand_mode = strand_mode,
                                  genome = genome_name,
                                  param = galp_param_mut))

      if (length(galp_mm) == 0) {
        message("No mutation-specific reads. Processing general alignments.")
      } else {
        galp <- c(galp_mm, galp)
        galp <- galp[!duplicated(galp@NAMES)]
        # Update metadata with content specifics
        msg1 <- "GRanges with general and mutation-specific read-pairs;"
        msg2 <- paste("Mutation-specific count:", length(galp_mm))
        galp@metadata <- list(paste(msg1, msg2))
      }
    }
  } else {
    # Process all alignments when no mutation file is provided
    galp <- bam_to_galp2(bamfile = bamfile, 
                         use_names = use_names,
                         chromosome_to_keep = chromosome_to_keep,
                         strand_mode = strand_mode,
                         genome = genome_name,
                         param = galp_param)
  }

  return(galp)
}

#' Create an Empty GAlignmentPairs Object
#'
#' Initializes an empty GAlignmentPairs object with two GAlignments objects,
#' one for each pair member, and combines them.
#'
#' @return An empty GAlignmentPairs object.
#' @import GenomicAlignments
#' @import GenomicRanges
#' @importFrom S4Vectors Rle
#' @noRd
create_empty_galp <- function() {
  # Create empty GAlignments for pair members
  empty_first <- GenomicAlignments::GAlignments(
    seqnames=character(), cigar=character()
  )
  empty_last <- GenomicAlignments::GAlignments(
    seqnames=character(), cigar=character()
  )
  
  # Combine into a GAlignmentPairs object
  empty_galp <- GenomicAlignments::GAlignmentPairs(
    first=empty_first, last=empty_last
  )
  
  return(empty_galp)
}




#' Function to determine the correct BSgenome based on the seqinfo result
#' @importFrom GenomeInfoDb seqinfo
#' @noRd
get_genome_reference <- function(frag_obj_mut) {
  # Extract genome sequence from GRanges object's seqinfo
  genome_seq <- unique(
    as.character(GenomeInfoDb::seqinfo(frag_obj_mut)@genome))
  
  # Define genome versions and corresponding BSgenome data packages
  genome_versions <- c("GRCh38", "hg38-NCBI", "hg38", "hg19")
  genome_packages <- c(
    "BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38", 
    "BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38", 
    "BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38", 
    "BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19")
  
  # Find the index of the genome sequence in the genome_versions vector
  index <- match(genome_seq, genome_versions)
  
  # Return the corresponding genome package; default to hg19 if not found
  if (!is.na(index)) {
    return(eval(parse(text = genome_packages[index])))
  } else {
    return(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  }
}

#' Prepare and filter GRanges data
#' @noRd
prepare_data <- function(frag_obj_mut) {
  frag_obj_mut <- frag_obj_mut[!grepl("outer_fragment",
                               frag_obj_mut$target_mutation), ]
  as.data.frame(frag_obj_mut)
}

#' Process mutation status and split into components
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract
#' @noRd
stratify_mutation_status <- function(gr_df) {
  gr_df %>%
    mutate(
      sbs_mut_status = str_extract(mutation_status, "^[^:]+"),
      sbs_read_overlap = str_extract(mutation_status, "(?<=:)(single_read|read_pair)"),
      sbs_read_type = str_extract(mutation_status, "(?<=_)(concordant|discordant)|[+-]"),
      sbs_read_bases = str_extract(mutation_status, "(?<=:)[ACGT][-ACGT/]+$")
    )
}

#' Function to check for mutational metadata in the GRanges object
#' @noRd
check_mutation_in_metadata <- function(frag_obj) {
  # Check if the input is a GRanges object
  if (!inherits(frag_obj, "GRanges")) {
    stop("The provided object is not a GRanges object.")
  }

  # Extract the metadata
  metadata_list <- frag_obj@metadata

  # Check for the specific metadata entry
  required_info <- "GRanges object with fragment and mutational data"

  if (!any(sapply(metadata_list,
                  function(x) grepl(required_info, x, ignore.case = TRUE)))) {
    stop("The GRanges object does not contain fragment data ",
         "with annotated mutational information.\n ",
         " Re-run readBam() with ",
         "parameters mutation_file = '/path/to/mutation_file.tsv'.\n ")
  }

  message("The GRanges object contains mutational information. ",
          "However, not all fragments may have the target mutant base; ",
          "some or all may contain only the reference base.")

}

#' Function to check for mutant base fragments in a dataframe
#' @noRd
check_mutation_status <- function(reads_df) {
  # Check if 'locus_status' contains the string "MUT"
  if (!any(grepl("MUT", reads_df$mutation_status, ignore.case = TRUE))) {
    stop("The GRanges object provided does not have any ",
         "fragments with mutant bases.\n  Please check your data.")
  }

  # Proceed with the rest of function if "MUT" is found
  message("Fragments with mutant bases found. Proceeding with the analysis.")
}

#' Helper function to update locus status
#' @importFrom dplyr mutate case_when select
#' @importFrom stringr str_sub str_extract str_replace
#' @noRd
update_locus_status <- function(gr_df) {
  gr_df %>%
    mutate(
      # Extract the initial part of the mutation before the "-"
      locus_status = sub("-.*$", "", target_mutation),

      # Determine the relevant base from the sbs_read_bases based on read type
      base = case_when(
        sbs_read_type %in% c("concordant", "+", "-") ~
          str_sub(sbs_read_bases, -1),  # Last character for these types
        sbs_read_type == "discordant" ~
          str_extract(sbs_read_bases, "(?<=-)[ACGT](?=/?)"),
        TRUE ~ NA_character_  # Default case for unexpected or missing data
      ),

      # Formulate the complete locus_status value
      locus_status = paste0(locus_status, ":", base, ":", sbs_mut_status)
    ) %>%
    mutate(
      # Adjust 'locus_status' for 'discordant' and 'OTHER' categories
      locus_status = case_when(
        sbs_read_type == "discordant" ~
          str_replace(locus_status, "MUT", "discordant"),
        sbs_mut_status == "OTHER" ~
          str_replace(locus_status, "OTHER", paste0(
            "other_base_", sbs_read_overlap)),
        TRUE ~ locus_status  # Keep existing values where no changes are needed
      )
    ) %>%
    select(-base)
}

#' Helper function to summarize mutational data
#' @noRd
summarize_mutational_data <- function(gr_df) {
  gr_df %>%
    group_by(target_mutation) %>%
    summarize(
      CO_MUT = sum(sbs_mut_status == "MUT" &
                   sbs_read_type == "concordant"),
      SO_MUT = sum(sbs_mut_status == "MUT" &
                   sbs_read_overlap == "single_read"),
      CO_REF = sum(sbs_mut_status == "REF" &
                   sbs_read_type == "concordant"),
      SO_REF = sum(sbs_mut_status == "REF" &
                   sbs_read_overlap == "single_read"),
      DO = sum((sbs_mut_status == "MUT" |
                sbs_mut_status == "OTHER") &
               sbs_read_type == "discordant"),
      SO_OTHER = sum(sbs_mut_status == "OTHER" &
                     sbs_read_overlap == "single_read"),
      CO_OTHER = sum(sbs_mut_status == "OTHER" &
                     sbs_read_type == "concordant")
    )
}

#' Helper function to summarize fragment lengths
#' @importFrom dplyr group_by summarize
#' @importFrom stats median
#' @noRd
summarize_fragment_lengths <- function(gr_df) {
  gr_df %>%
    group_by(target_mutation) %>%
    summarize(
      CO_MUT_flength = median(width[sbs_mut_status == "MUT" &
                                    sbs_read_type == "concordant"],
                              na.rm = TRUE),
      SO_MUT_flength = median(width[sbs_mut_status == "MUT" &
                                    sbs_read_overlap == "single_read"],
                              na.rm = TRUE),
      CO_REF_flength = median(width[sbs_mut_status == "REF" &
                                    sbs_read_type == "concordant"],
                              na.rm = TRUE),
      SO_REF_flength = median(width[sbs_mut_status == "REF" &
                                    sbs_read_overlap == "single_read"],
                              na.rm = TRUE),
      DO_flength = median(width[(sbs_mut_status == "MUT" |
                                 sbs_mut_status == "OTHER") &
                                sbs_read_type == "discordant"],
                          na.rm = TRUE),
      SO_OTHER_flength = median(width[sbs_mut_status == "OTHER" &
                                      sbs_read_overlap == "single_read"],
                                na.rm = TRUE),
      CO_OTHER_flength = median(width[sbs_mut_status == "OTHER" &
                                      sbs_read_type == "concordant"],
                                na.rm = TRUE)
    )
}

#' Function to update consensus_mismatch based on the highest priority match
#' @importFrom dplyr filter
#' @noRd
update_consensus_mismatch <- function(merged_table_df, gr_df, mapping) {
  set.seed(123)  # For reproducibility in random selection
  
  # Iterate over each row in merged_table_df to update consensus_mismatch
  for (i in 1:nrow(merged_table_df)) {
    target_mutation <- merged_table_df$target_mutation[i]
    selected_column <- get_highest_column(merged_table_df[i, ])

    # Define the required status based on the selected column
    required_status <- mapping[[selected_column]]

    # Find matching rows in gr_df by target_mutation and filter
    filtered_rows <- gr_df %>%
      filter(target_mutation == !!target_mutation) %>%
      filter(
        (required_status == "MUT" & sbs_mut_status == "MUT") |
        (required_status == "discordant" & sbs_read_type == "discordant") |
        (required_status == "other_base_single_read" &
         sbs_mut_status == "OTHER" & sbs_read_overlap == "single_read") |
        (required_status == "other_base_read_pair" &
         sbs_mut_status == "OTHER" & sbs_read_type == "concordant")
      )

    # Randomly select one row if multiple matches exist
    if (nrow(filtered_rows) > 1) {
      selected_row <- filtered_rows[sample(1:nrow(filtered_rows), 1), ]
    } else {
      selected_row <- filtered_rows
    }

    # Update consensus_mismatch column if a match is found
    if (!is.null(selected_row) && nrow(selected_row) > 0) {
      merged_table_df$consensus_mismatch[i] <- selected_row$locus_status
    } else {
      merged_table_df$consensus_mismatch[i] <- NA  # Ensure NAs for unmatched
    }
  }

  # Filter out rows without a valid consensus_mismatch
  merged_table_df <- merged_table_df[complete.cases(
    merged_table_df$consensus_mismatch), ]
  return(merged_table_df)
}

#' Helper function which selects consensus mutation type based on counts
#' @noRd
get_highest_column <- function(row) {
  values <- row[c("CO_MUT", "SO_MUT", "DO", "SO_OTHER", "CO_OTHER")]
  max_value <- max(values)
  highest_columns <- names(values)[values == max_value]

  # If the choice is between SO_OTHER, CO_OTHER,
  # and DO with equal values, randomly choose one
  if (length(highest_columns) > 1) {
    selected_column <- sample(highest_columns, 1)
  } else {
    selected_column <- highest_columns
  }

  # Always prioritize CO_MUT or SO_MUT
  if ("CO_MUT" %in% highest_columns) {
    selected_columns <- "CO_MUT"
  }
  if ("SO_MUT" %in% highest_columns) {
    selected_columns <- "SO_MUT"
  }
  return(selected_column)
}

#' Helper function which generates SBS 96 trinucleotide channel information
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr filter select
#' @importFrom stats setNames
#' @noRd
get_trinucleotide <- function(merged_table_df, genome) {
  # Split consensus_mismatch into three fields: chr, start, ALT Base
  split_values <- strsplit(merged_table_df$consensus_mismatch, ":")

  # Form a dataframe out of chr,start,end required for GRanges to get ref bases
  tri_df <- data.frame(
    chromosome = sapply(split_values, "[", 1),
    start = as.integer(sapply(split_values, "[", 2)) - 1,
    end = as.integer(sapply(split_values, "[", 2)) + 1
  )

  # Form a GRanges to obtain reference trinucleotide
  tri_loci_gr <- GRanges(seqnames = tri_df$chromosome,
                         ranges = IRanges(start = tri_df$start,
                                          end = tri_df$end))

  # Get Vector of trinucleotide
  ref_tri_bases <- as.vector(BSgenome::getSeq(genome, names = tri_loci_gr))

  # Add reference trinucleotide to table
  merged_table_df$ref_tri <- ref_tri_bases
  merged_table_df$mut_base <- sapply(split_values, "[", 3)
  merged_table_df$ref_base <- substr(merged_table_df$ref_tri, 2, 2)

  # Split df into two parts: A/G and C/T
  merged_table_df_1 <- merged_table_df %>% filter(ref_base %in% c("A", "G"))
  merged_table_df_2 <- merged_table_df %>% filter(ref_base %in% c("C", "T"))

  # Reverse complement for A/G
  merged_table_df_1$ref_tri <- chartr("ATGC", "TACG",
                                      merged_table_df_1$ref_tri)
  merged_table_df_1$ref_base <- chartr("ATGC", "TACG",
                                       merged_table_df_1$ref_base)
  merged_table_df_1$mut_base <- chartr("ATGC", "TACG",
                                       merged_table_df_1$mut_base)

  merged_table_df <- rbind(merged_table_df_1, merged_table_df_2)

  # Construct SBS96 notation
  merged_table_df$SBS96 <- paste0(
    merged_table_df$left <- substr(merged_table_df$ref_tri, 1, 1),
    "[",
    merged_table_df$ref_base,
    ">",
    merged_table_df$mut_base,
    "]",
    merged_table_df$right <- substr(merged_table_df$ref_tri, 3, 3)
  )

  # Add chr and pos columns
  merged_table_df$chr <- sapply(split_values, "[", 1)
  merged_table_df$pos <- as.integer(sapply(split_values, "[", 2))

  # Clean up the DataFrame
  merged_table_df <- merged_table_df %>%
    select(-mut_base, -ref_base, -ref_tri, -left, -right, -chr, -pos)

  return(merged_table_df)
}

#' Function which generates complete SBS96 profiles with readpair overlap info
#' Function which generates complete SBS96 profiles with readpair overlap info
#' @importFrom dplyr mutate across arrange select
#' @importFrom tidyr pivot_longer
#' @importFrom stats setNames
#' @importFrom stringr str_extract
#' @importFrom rlang .data
#' @noRd
processTrinucleotideData <- function(trinuc_df,
                                     exclude_if_type_present = NULL,
                                     retain_if_type_present = NULL,
                                     remove_type = NULL,
                                     normalize_counts = TRUE) {

  # Exclude rows based on exclude_if_type_present
  if (!is.null(exclude_if_type_present) &&
        length(exclude_if_type_present) > 0) {
    trinuc_df <- trinuc_df[!apply(trinuc_df[exclude_if_type_present] > 0,
                                  1, any), ]
  }

  # Retain rows based on retain_if_type_present
  if (!is.null(retain_if_type_present) && length(retain_if_type_present) > 0) {
    trinuc_df <- trinuc_df[apply(trinuc_df[retain_if_type_present] > 0,
                                 1,
                                 any), ]
  }

  # Convert all counts of specified types to 0 in remove_type columns
  if (!is.null(remove_type) && length(remove_type) > 0) {
    trinuc_df[remove_type] <- lapply(trinuc_df[remove_type],
                                     function(x) ifelse(x > 0, 0, x))
  }

  # Define the relevant columns to check for the maximum value
  relevant_columns <- c("SO_MUT", "CO_MUT", "DO", "SO_OTHER", "CO_OTHER")

  # Define a function to find the column name with the highest value
  get_max_column <- function(x) {
    # Find the maximum value among the specified columns
    max_val <- max(x[relevant_columns])

    # Get all column names that have this maximum value
    max_cols <- names(x[relevant_columns])[x[relevant_columns] == max_val]

    # If more than one column shares the max value, randomly choose one
    if (length(max_cols) > 1) {
      return(sample(max_cols, 1))
    } else {
      return(max_cols)
    }
  }

  # Apply this function across all rows of the dataframe
  trinuc_df$consensus_mismatch_type <- apply(trinuc_df[, relevant_columns],
                                             1,
                                             get_max_column)

  # Filter out rows where all values in relevant columns are zero
  trinuc_df <- trinuc_df[rowSums(trinuc_df[, relevant_columns] > 0) > 0, ]

  # Create a 96 channel trinucleotide matrix from the 6 SBS types
  mutation.types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  bases <- c("A", "T", "C", "G")
  combinations <- expand.grid(start = bases, mut = mutation.types, end = bases)
  sub.types.96 <- with(combinations, paste0(start, "[", mut, "]", end))

  # Count occurrences of each subtype in trinuc_df$SBS96
  subtype_counts <- table(factor(trinuc_df$SBS96, levels = sub.types.96))
  count_df <- data.frame(sample = as.integer(subtype_counts))
  rownames(count_df) <- sub.types.96

  # Count the consensus_mismatch_type occurrences for each SBS96
  for (subtype in sub.types.96) {
    subtype_rows <- trinuc_df[trinuc_df$SBS96 == subtype, ]
    count_df[subtype, "CO_MUT"] <- sum(
      subtype_rows$consensus_mismatch_type == "CO_MUT")
    count_df[subtype, "SO_MUT"] <- sum(
      subtype_rows$consensus_mismatch_type == "SO_MUT")
    count_df[subtype, "DO"] <- sum(
      subtype_rows$consensus_mismatch_type == "DO")
    count_df[subtype, "SO_OTHER"] <- sum(
      subtype_rows$consensus_mismatch_type == "SO_OTHER")
    count_df[subtype, "CO_OTHER"] <- sum(
      subtype_rows$consensus_mismatch_type == "CO_OTHER")
  }

  #Add the counts of SO_OTHER_sample to SO_MUT_sample
  count_df$SO_MUT <- count_df$SO_MUT + count_df$SO_OTHER

  # Add the counts of CO_OTHER_sample to CO_MUT_sample
  count_df$CO_MUT <- count_df$CO_MUT + count_df$CO_OTHER

  # Remove the SO_OTHER_sample and CO_OTHER_sample columns
  count_df <- count_df[, !colnames(count_df) %in% c("SO_OTHER", "CO_OTHER")]

  # Normalize counts if normalize_counts is TRUE
  if (normalize_counts) {
    total_sample_sum <- sum(count_df$sample)
    count_df <- count_df %>%
      mutate(across(c(CO_MUT, SO_MUT, DO), ~ .x / total_sample_sum))
  }

  # Remove the sample column
  count_df <- count_df[, !colnames(count_df) %in% c("sample")]

  # Add the SBS names as a column and convert to long format
  count_df$SBS <- rownames(count_df)
  count_df <- pivot_longer(count_df, cols = c(CO_MUT, SO_MUT, DO),
                           names_to = "overlap_type",
                           values_to = "value")

  # Mutation types in the desired order
  mutation.types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

  # Function to extract the mutation type from the SBS column
  extract_mutation_type <- function(s) {
    sub(".*\\[([A-Z]>[A-Z])\\].*", "\\1", s)
  }

  # Add a column for the extracted mutation type
  count_df <- count_df %>%
    mutate(mutation_type = extract_mutation_type(SBS))

  # Convert the mutation_type column to a factor with the specified levels
  count_df <- count_df %>%
    mutate(mutation_type = factor(mutation_type, levels = mutation.types))

  # Arrange the dataframe based on the mutation_type factor
  count_df <- count_df %>%
    arrange(mutation_type) %>%
    select(-mutation_type)  # Remove the helper column

  count_df <- as.data.frame(count_df)

  return(count_df)
}

#' Function which runs the the plotting operations for SBS96 profile
#' @import ggplot2
#' @importFrom ggpattern geom_bar_pattern
#' @importFrom patchwork plot_layout wrap_elements
#' @importFrom grid unit
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select filter
#' @importFrom stringr str_extract str_replace
#' @importFrom scales percent
#' @noRd
plotTrinucData <- function(
    count_df,
    ylim = ylim,
    show_overlap_type = TRUE,
    plot_title = "Trinucleotide Profile",
    y_axis_title = "Proportion of Single Base Substitutions",
    draw_x_axis_labels = TRUE,
    draw_y_axis_labels = TRUE,
    draw_y_axis_title = TRUE,
    output_file = output_file,
    ggsave_params = ggsave_params) {

  # Define the mutation types
  mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

  # Define Plot colours
  plot.colors <- c("#232f7c", "#000000",
                   "#612370", "#474242",
                   "#436e2d", "#ad5c74")

  strip.text.x.colors <- c("white", "white", "white", "white", "white",
                           "white")

  # Initialize a list to hold data frames for each mutation type
  df.plot.data <- list()

  # Populate the list with segmented data frames for each mutation type
  for (i in 1:6) {
    start.idx <- ((i - 1) * 48) + 1
    end.idx <- i * 48

    # Create a temporary data frame with sliced data
    df.temp <- data.frame(
      mutation_subtype = count_df$SBS[start.idx:end.idx],
      value = count_df$value[start.idx:end.idx],
      overlap_type = count_df$overlap_type[start.idx:end.idx],
      title = rep(mutation_types[i], 48),
      stringsAsFactors = FALSE
    )

    # Add the temporary data frame to the list
    df.plot.data[[i]] <- df.temp
  }

  plots <- list()

  # Define axis text settings based on input parameters
  axis.text.x <- if (draw_x_axis_labels) {
    element_text(size = 4,
                 family = "Helvetica",
                 colour = "#888888",
                 angle = 90,
                 vjust = 0.5,
                 hjust = 0.5)
  } else {
    element_blank()
  }

  axis.text.y <- if (draw_y_axis_labels) {
    element_text(size = 5,
                 family = "Helvetica",
                 colour = "#888888")
  } else {
    element_blank()
  }

  # Extract unique 'overlap_type' values from all data frames in the list
  # where the 'value' column is not zero
  overlap_levels <- unique(unlist(lapply(df.plot.data, function(df) {
    # Filter out rows where the 'value' column is zero
    non_zero_df <- df[df$value != 0, ]
    # Return 'overlap_type' values from filtered data
    non_zero_df$overlap_type
  })))

  # Predefined order of levels
  preferred_order <- c("DO", "CO_MUT", "SO_MUT")

  # Filter the levels based on their presence in the extracted unique values
  overlap_levels <- preferred_order[preferred_order %in% overlap_levels]

  # Loop through each data frame in the list and create plots
  for (i in seq_along(df.plot.data)) {
    df.plot.data[[i]]$overlap_type <- factor(df.plot.data[[i]]$overlap_type,
                                             levels = overlap_levels)
  }

  # Assuming df.plot.data is a list of data frames
  df.plot.data <- lapply(df.plot.data, function(df) {
    # Filter out rows where 'overlap_type' is NA
    df <- df[!is.na(df$overlap_type), ]
    return(df)
  })

  # Initialise the plot
  plots <- list()

  # Compile the plot
  for (i in 1:6) {
    if (show_overlap_type == TRUE) {
      p <- ggplot(df.plot.data[[i]], aes_string(x = 'mutation_subtype',
                                                y = 'value',
                                                pattern = 'overlap_type')) +

        ggpattern::geom_bar_pattern(stat = "identity", width = 0.9,
                        fill = plot.colors[i], pattern_density = 0.005,
                        pattern_spacing = 0.05, pattern_colour = "#07d942",
                        pattern_fill = "#07d942") +

        scale_pattern_manual(
          values = c("CO_MUT" = "none",
                    "SO_MUT" = "stripe",
                    "DO" = "crosshatch"),
          labels = c("CO_MUT" = "Concordant",
                    "SO_MUT" = "Single-Read",
                    "DO" = "Discordant")) +

        ylab(NULL) + xlab(NULL) +
        scale_y_continuous(limits = ylim,
                          expand = c(0, 0)) +
        theme(axis.text.x = axis.text.x,
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.y = element_line(colour = "#DCDCDC",
                                                linetype = "solid"),
              panel.grid.major.x = element_blank(),
              panel.background = element_blank())

        # Conditional logic to hide legend except for the 4th plot
        if (i != 4) {
          p <- p + theme(legend.position = "none")
        } else {
          # Customize the legend appearance for the 4th plot
          p <- p +
            guides(pattern = guide_legend(
              title = "Overlap Type",
              title.position = "top", label.position = "right",
              override.aes = list(fill = "#f1f5f2"))) +
            theme(
              legend.title = element_text(size = 6, face = "bold"),
              legend.text = element_text(size = 5),
              legend.key.size = unit(0.63, "cm"),
              legend.spacing = unit(0.3, "cm"),
              legend.margin = margin(0.1, 0.1, 0.1, 0.1),
              legend.box.margin = margin(0.1, 0.1, 0.1, 0.1)
            )
        }
    } else if (show_overlap_type == FALSE) {

      # Compile the plot
      p <- ggplot(df.plot.data[[i]],
                  aes_string(x = 'mutation_subtype', y = 'value')) +
        geom_col(fill = plot.colors[i], width = 0.9)

      # Conditional logic to hide legend for all plots
      p <- p + theme(legend.position = "none")

      p <- p + ylab(NULL) + xlab(NULL) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        theme(axis.text.x = axis.text.x,
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major.y = element_line(colour = "#DCDCDC",
                                                linetype = "solid"),
              panel.grid.major.x = element_blank(),
              panel.background = element_blank())
    }

    # Draw top strip
    p <- p +
      theme(strip.background = element_rect(fill = plot.colors[i],
                                            colour = plot.colors[i]),
            strip.text.x = element_text(size = 7,
                                        family = "Helvetica",
                                        face = "bold",
                                        colour = strip.text.x.colors[i])) +
      facet_grid(~title)

    # Margins for the left-most mutation type sub-plot
    if (i == 1) {
      p <- p + theme(axis.text.y = axis.text.y,
                     plot.margin = unit(c(0.3,
                                          0,
                                          0.3,
                                          0.3), "cm"))

      # Margins for the right-most mutation type sub-plot
    } else if(i == 6) {
      p <- p + theme(axis.text.y = element_text(
                                                size = 5,
                                                family = "Helvetica",
                                                colour = "#FFFFFF00"),
      plot.margin = unit(c(0.3,
                           0.3,
                           0.3,
                           0), "cm"))

      # Margins for the middle mutation types
    } else {
      p <- p + theme(axis.text.y = element_text(size = 5,
                                                family = "Helvetica",
                                                colour = "#FFFFFF00"),
                     plot.margin = unit(c(0.3,
                                          0,
                                          0.3,
                                          0), "cm"))
    }
    plots[[i]] <- p
  }

  # Patchwork combine the plots
  plot_combined <- plots[[1]] + plots[[2]] +
    plots[[3]] + plots[[4]] +
    plots[[5]] + plots[[6]] +
    plot_layout(guides = 'collect', ncol = 6) &
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))

  # Use the tag label as a y-axis label
  p_output <- wrap_elements(plot_combined) +
    labs(tag = y_axis_title) +
    theme(
      plot.tag = element_text(size = rel(0.5), angle = 90),
      plot.tag.position = "left"
    )

  # Check if 'output_file' is not NULL and save the plot to the specified path
  if (!is.null(output_file)) {
    do.call("ggsave",
            c(list(plot = p_output, filename = output_file),
            ggsave_params))
    message("Plot saved to ", output_file)
  }

  return(p_output)

}

#' Helper function to process mutation-related information
#' @importFrom dplyr filter mutate group_by summarize ungroup as_tibble n
#' @importFrom plyr round_any
#' @importFrom stats setNames
#' @noRd
process_length_mut <- function(frag_obj, ref_type, downsample_ref) {
    # Ensure mutation data is checked
    check_mutation_in_metadata(frag_obj = frag_obj)
    check_mutation_status(as.data.frame(frag_obj))

    # Filter out discordant and non-MUT base (other) fragments
    size_table_inner <- dplyr::filter(as.data.frame(frag_obj), 
                                      !grepl("discordant", mutation_status), 
                                      grepl("MUT", mutation_status))

    # Determine the count for downsampling or normalization
    target_count <- nrow(size_table_inner)

    # Select the correct reference type fragments and normalize if necessary
    size_table_outer <- select_and_normalize(frag_obj = frag_obj,
                                             pattern = ifelse(
                                              ref_type == "locus_fragment",
                                              "REF", ref_type),
                                             downsample_ref = downsample_ref,
                                             target_count = target_count)

    # Aggregate and summarize data
    summary_size_table_outer <- dplyr::as_tibble(size_table_outer) %>%
                                dplyr::group_by(width) %>%
                                dplyr::summarize(count = n(), .groups = "drop")

    summary_size_table_inner <- dplyr::as_tibble(size_table_inner) %>%
                                dplyr::group_by(width) %>%
                                dplyr::summarize(count = n(), .groups = "drop")

    summary_size_table_inner$MUTANT <- "true"
    summary_size_table_outer$MUTANT <- "false"

    size_merged_df <- merge(summary_size_table_outer,
                            summary_size_table_inner,
                            all = TRUE)

    colnames(size_merged_df) <- c("SIZE", "COUNT", "MUTANT")

    sizeCharacterisationSummary <- size_merged_df %>%
      mutate(SIZE.ROUNDED = plyr::round_any(SIZE, accuracy = 5L)) %>%
      group_by(MUTANT, SIZE.ROUNDED) %>%
      summarise(COUNT = sum(COUNT), .groups = "drop_last") %>%
      ungroup() %>%
      mutate(TOTAL = sum(COUNT)) %>%
      mutate(PROPORTION = COUNT / TOTAL)

    result <- sizeCharacterisationSummary %>%
      mutate(MUTANT_LABEL = as.factor(ifelse(MUTANT,
                                             "Mutation Fragment",
                                             "Reference Fragment")))

    # Add an attribute named "mutational_info" with value TRUE
    attr(result, "mutational_info") <- TRUE

    return(result)
}

#' Helper function to select and optionally normalize REF fragment counts
#' @importFrom dplyr filter sample_n
#' @noRd
select_and_normalize <- function(frag_obj,
                                 pattern,
                                 downsample_ref = FALSE,
                                 target_count) {

  # Select the relevant reference fragments
  table <- as.data.frame(frag_obj) %>% filter(grepl(pattern,
                                                    mutation_status,
                                                    ignore.case = TRUE))
  # Downsample selected reference fragments to match the MUT fragment count
  if (downsample_ref) {
    set.seed(123)
    table <- dplyr::sample_n(table, target_count)
  }
  table
}

#' Function which generates fragment length plot with mutational data
#' @import ggplot2
#' @importFrom dplyr pull
#' @importFrom grid unit
#' @noRd
plot_length_mut <- function(x, ylim, output_file, ggsave_params) {
    # Plot Fragment lengths with integrated mutational information
    max_fraction <- max(dplyr::pull(x, PROPORTION))

    if (missing(ylim)) {
      ylim <- c(0, max_fraction * 1.1)
      message("ylim", " was set as: ", ylim[1], " - ", ylim[2])
    }

    p <- ggplot(data = x, aes(x = SIZE.ROUNDED,
                              y = PROPORTION,
                              fill = MUTANT_LABEL)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "cfDNA Fragment Length (bp)",
           y = "Proportion",
           fill = "cfDNA Fragment Type") +
      theme_classic() +
      theme(axis.text = element_text(size = 5),
            axis.title = element_text(size = 6, face = "bold"),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 7, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.position = c(0.85, 0.5)) +
      scale_x_continuous(limits = c(0.5, 500)) +
      geom_vline(xintercept = c(166, 166 * 2), linetype = "dashed",
                 alpha = 0.5) +
      annotate("text", x = 200, y = max(x$PROPORTION) * 1.1,
               label = "166bp", alpha = 0.5,
               size = 2) +
      annotate("text", x = 370, y = max(x$PROPORTION) * 1.1,
               label = "332bp", alpha = 0.5,
               size = 2) +
      scale_fill_manual(values = c("Mutation Fragment" = "#00BFC4",
                                   "Reference Fragment" = "#F8766D"))

  # Check if 'output_file' is not NULL and save the plot to the specified path
  if (!is.null(output_file)) {
    do.call("ggsave", c(list(plot = p, filename = output_file), ggsave_params))
    message("Plot saved to ", output_file)

  }
  return(p)
}

#' Helper function to process motif and mutation data
#' @importFrom dplyr mutate select full_join rename
#' @importFrom Biostrings DNAStringSet
#' @noRd
integrate_motif_mut <- function(frag_obj, downsample_ref, ref_type,
                                bsgenome_obj, motif_type, motif_length) {
  # Check if GRanges object has mutational data
  check_mutation_in_metadata(frag_obj = frag_obj)

  check_mutation_status(as.data.frame(frag_obj))

  # User defined parameter for downsampling fragments
  if (downsample_ref == TRUE) {
    # Downsample refeference base fragments to match mutation base fragments
    frag_obj <- adjust_ref_fragments(gr = frag_obj,
                                     ref_type = ref_type,
                                     downsample_ref = downsample_ref)
  }

  # Extract references and mutations based on mutation_status
  gr_ref <- frag_obj[grepl("REF|outer_fragment", frag_obj$mutation_status)]
  gr_mut <- frag_obj[grepl("MUT:read_pair_concordant|MUT:single_read",
                           frag_obj$mutation_status)]

  # Process motifs for reference and mutation fragments
  motif_calls_ref <- processMotif(frag_obj = gr_ref,
                                  bsgenome_obj = bsgenome_obj,
                                  motif_type = motif_type,
                                  motif_length = motif_length)

  motif_calls_mut <- processMotif(frag_obj = gr_mut,
                                  bsgenome_obj = bsgenome_obj,
                                  motif_type = motif_type,
                                  motif_length = motif_length)

  # Join reference and mutation motif calls
  result_frac <- full_join(motif_calls_ref, motif_calls_mut, by = "motif")

  # Rename with specific column names
  result_frac <- result_frac %>%
    dplyr::rename(
      motif = motif,
      n_ref = n.x,
      fraction_ref = fraction.x,
      n_mut = n.y,
      fraction_mut = fraction.y
    )

  # Recalculate fractional values
  result_frac <- result_frac %>%
    # Compute total sums for n_ref and n_mut
    mutate(
      total_n_ref = sum(n_ref),
      total_n_mut = sum(n_mut)
    ) %>%
  # Calculate fraction_ref and fraction_mut
    mutate(
      fraction_ref = n_ref / (total_n_ref + total_n_mut),
      fraction_mut = n_mut / (total_n_ref + total_n_mut)
    ) %>%
    # Remove the temporary total columns if they are no longer needed
    select(-total_n_ref, -total_n_mut)

  # Add an attribute named "mutational_info" with value TRUE
  attr(result_frac, "mutational_info") <- TRUE

  return(result_frac)
}

#' Helper function for downsampling reference base fragments in motifs
#' @importFrom S4Vectors mcols
#' @noRd
adjust_ref_fragments <- function(gr, ref_type, downsample_ref = FALSE) {
  # Return the original object if no downsampling is needed
  if (!downsample_ref) {
    return(gr)
  }

  # Set the seed for reproducibility
  set.seed(123)

  # Initialize ref_string based on ref_type
  if (ref_type == "locus_fragment") {
    ref_string <- "REF"
  } else if (ref_type == "outer_fragment") {
    ref_string <- "outer_fragment"
  }

  # Identify rows with REF (or outer_fragment) and MUT in the mutation_status
  ref_indices <- grep(ref_string, mcols(gr)$mutation_status)
  mut_indices <- grep("MUT", mcols(gr)$mutation_status)

  # Check if ref_indices is empty
  if (length(ref_indices) == 0) {
    stop(paste("Invalid use of the ref_type argument.\n",
               "Please verify the types of fragments",
               "present in your GRanges object."))
  }

  # Count the number of "REF" and "MUT" entries
  n_ref <- length(ref_indices)
  n_mut <- length(mut_indices)

  # Downsample "REF" entries if there are more "REF" than "MUT"
  if (n_ref > n_mut) {
    # Sample from ref_indices to match the number of "MUT" entries
    ref_indices_to_keep <- sample(ref_indices, n_mut)
    # Combine kept "REF" indices with "MUT" indices
    final_indices <- c(ref_indices_to_keep, mut_indices)
  } else {
    # If not downsampling, keep all indices
    final_indices <- c(ref_indices, mut_indices)
  }

  # Subset the GRanges object to include only the final indices
  return(gr[final_indices])
}

#' Helper function for plotting fragment motif and mutational data
#' Helper function for plotting fragment motif and mutational data
#' @import ggplot2
#' @importFrom dplyr mutate filter arrange
#' @importFrom stringr str_to_lower str_to_title str_extract
#' @importFrom tidyr pivot_longer
#' @importFrom grDevices adjustcolor
#' @noRd
plot_motif_mut <- function(
  x, ylim, x_title, plot_type, bar_color, motif_levels,
  output_file, ggsave_params) {

  # Print message
  message(paste(
    " The provided GRanges fragment object contains",
    "mutational information.\n",
    "The plot will be adapted to reflect the mutational fragment information."
  ))

  x <- tibble::as_tibble(x)
  plot_type <- stringr::str_to_lower(plot_type)

  # Create color mapping for REF and MUT using adjustcolor
  ref_colors <- setNames(adjustcolor(bar_color, alpha.f = 0.5),
                         paste0(names(bar_color), "_REF"))
  mut_colors <- setNames(adjustcolor(bar_color, alpha.f = 1.5),
                         paste0(names(bar_color), "_MUT"))
  all_colors <- c(ref_colors, mut_colors)

  # Reshape the data for stacking
  if (plot_type == "fraction") {
    x <- x %>%
      pivot_longer(cols = starts_with("fraction"),
                   names_to = "type", values_to = "value")
  } else {
    x <- x %>%
      pivot_longer(cols = starts_with("n_"),
                   names_to = "type", values_to = "value")
  }

  # Set ylim if not provided
  if (is.null(ylim)) {
    max_val <- max(x$value, na.rm = TRUE)
    ylim <- c(0, max_val * 1.2)
    message("ylim", " was set as: ", ylim[1], " - ", ylim[2])
  }

  # Extract the first letter of motif for color grouping and adjust type
  x <- x %>%
    mutate(group = stringr::str_extract(motif, "^[ATCG]"),
           type = ifelse(type == "fraction_ref" | type == "n_ref",
                         paste0(group, "_REF"), paste0(group, "_MUT"))) %>%
    filter(group %in% motif_levels) %>%
    arrange(match(group, motif_levels))  # Ensure correct ordering

  x$group <- factor(x$group, levels = motif_levels)
  x$motif <- factor(x$motif, levels = unique(x$motif))

  # Plotting with adjusted theme
  p <- ggplot(x, aes(x = motif, y = value, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = x_title, y = stringr::str_to_title(plot_type), fill = "Type") +
    scale_fill_manual(values = all_colors,
                      breaks = c("A_MUT", "A_REF"),
                      labels = c("MUT", "REF")) +
    scale_y_continuous(limits = ylim) +
    theme_motif_plot_mut() +
    guides(fill = guide_legend(
      title = "Motif Type",
      override.aes = list(fill = c("#757373", "lightgrey"))
  ))

  # Check if 'output_file' is not NULL and save the plot to the specified path
  if (!is.null(output_file)) {
    do.call("ggsave", c(list(plot = p, filename = output_file), ggsave_params))
    message("Plot saved to ", output_file)
  }

  return(p)
}

#' Helper function to process mutational information for CNV analysis
#' @importFrom plyranges join_overlap_left
#' @importFrom dplyr group_by_at summarize ungroup
#' @noRd
process_cnv_mut <- function(olap, frag_obj_mut) {
  # Join gene ranges with overlapping MUT/REF fragments
  olap_mut <- plyranges::join_overlap_left(olap, frag_obj_mut)

  # Convert to data frame for easier manipulation
  olap_mut_df <- as.data.frame(olap_mut)

  # Summarize fragment counts per gene
  olap_df <- suppressMessages({
    olap_mut_df %>%
      group_by_at(vars(1:(ncol(olap_mut_df) - 1))) %>%
      summarize(
        MUT_count = sum(grepl("MUT:read_pair_concordant|MUT:single_read",
                        mutation_status)),
        ALL_count = sum(grepl("[A-Za-z]", mutation_status))
      ) %>%
      ungroup() %>%
      as.data.frame()
  })

  return(olap_df)
}

#' Helper function for plotting CNV and mutational data
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 aes
#' @noRd
annotate_cnv_mut <- function(p, olap_df, segment_line_size,
                             output_file, ggsave_params, ...) {
  cat("CN plot with integrated mutational information:\n",
      "GRanges object holds data of cfDNA fragments.\n",
      "Fragments overlapping specified genes will be reported.\n",
      "Fragments with specific SNV bases are called Mutation Fragments.\n")

  if (any(olap_df$ALL_count > 0)) {

    cat("Overlapping fragments were detected for selected genes.\n")

    p <- p + geom_text_repel(data = olap_df,
      aes(x = .data$x_index,
          y = .data$logratio,
          label = paste(.data$SYMBOL, " Gene: ",
                        .data$MUT_count, "/",
                        .data$ALL_count,
                        " Fragments Carrying Candidate Mutations",
                        sep = "")),
      box.padding = 0.5,
      min.segment.length = 0,
      segment.linetype = 1,
      segment.size = segment_line_size / 2,
      segment.color = "black",
      max.overlaps = Inf,
      ...
    )
  } else {

    cat("No overlapping fragments detected for selected genes.\n",
        "Ensure gene identifiers are correct.\n",
        "Review or generate a new GRanges object.\n")

    p <- p + geom_text_repel(data = olap_df,
      aes(x = .data$x_index,
          y = .data$logratio,
          label = paste(.data$SYMBOL, " Gene: ",
                        "0/0",
                        " Fragments Carrying Candidate Mutations",
                        sep = "")),
      box.padding = 0.5,
      min.segment.length = 0,
      segment.linetype = 1,
      segment.size = segment_line_size / 2,
      segment.color = "black",
      max.overlaps = Inf,
      ...
    )
  }
  # Check if 'output_file' is not NULL and save the plot to the specified path
  if (!is.null(output_file)) {
    do.call("ggsave", c(list(plot = p, filename = output_file), ggsave_params))
    message("Plot saved to ", output_file)
  }

  return(p)
}