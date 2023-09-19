  mutationFile <- function(
                      bamfile,
                      use_names,
                      chromosome_to_keep,
                      strand_mode,
                      genome_label,
                      galp_flag,
                      galp_what,
                      galp_tag,
                      galp_mapqFilter,
                      mutation_file,
                      outer_downstream_lim,
                      outer_upstream_lim,
                      frag
                      ) {
      # Define mutation coordinates for readGAlignements()
      bed_full <- readMutationFile(mutation_file)

      bed_full <- subset(bed_full, chr %in% chromosome_to_keep)

      target_muts <- bed_full %>% select(chr, start, ref, alt)

      target_muts$target <- paste(target_muts$chr,
                                  target_muts$start,
                                  target_muts$ref,
                                  target_muts$alt,
                                  sep = ":")

      target_muts$locus <- paste(target_muts$chr, target_muts$start, sep = ":")

      target_muts <- target_muts %>% select(locus, target)

      which_processed <- make_granges(muts = bed_full)

      mismatch_df <- processMismatches(
                      which_processed = which_processed,
                      bed_full = bed_full,
                      bamfile = bamfile,
                      galp_flag = galp_flag,
                      galp_mapqFilter = galp_mapqFilter)


      # match the values and if they match by coordinate
      # add target_muts$target as target column 'target_mutation' in mismatch_df

      mismatch_df$locus <- sub(
                                "^([^:]+:[^:]+):.*",
                                "\\1",
                                mismatch_df$locus_info
                                )

      # Match values and add 'target_mutation' column to the first dataframe
      mismatch_df <- subset(mismatch_df, !grepl("chrX", locus_info))

      mismatch_df <- merge(mismatch_df, target_muts, by = "locus", all.x = TRUE)

      # GRanges to dataframe
      # This adds .1 to row names whenever there is a duplicate fragment ID
      metadata <- as.data.frame(names(frag))

      colnames(metadata) <- "fragment_event"

      metadata$fragment_event <- make.unique(
                                metadata$fragment_event, sep = ".")

      # Allow duplicate row names by adding .1 .2 .3
      fragment_gr_df <- suppressWarnings(
          as.data.frame(
            frag, use.outer.mcols = TRUE, row.names = metadata$fragment_event))

      fragment_gr_df$fragment_id <- rownames(fragment_gr_df)

      # Merge with fragmentomic annotation
      merged_df <- merge(
                      fragment_gr_df,
                      mismatch_df,
                      by = "fragment_id",
                      all.x = TRUE)

      merged_gr <- makeGRangesFromDataFrame(merged_df,
                                            keep.extra.columns = TRUE)
      gr_meta <- mcols(merged_gr)
      gr_meta$locus_info[is.na(gr_meta$locus_info)] <- "outer_fragment"
      gr_meta$locus_status[is.na(gr_meta$locus_status)] <- "outer_fragment"
      mcols(merged_gr) <- gr_meta

      return(merged_gr)

      # END, below is old code
      # need to still remove upstream/downstream coordinate params 
      # as these will be used in later functions for callX() and plotX()
      ##################################################################

  }


# Helper function - keep it
# mage granges for mismatch and  mutation positions
make_granges <- function(muts) {
  which_processed <- GRanges(seqnames = muts[, 1],
          ranges = IRanges(start = muts[, 2], end = muts[, 3]))

  return(which_processed)
}
