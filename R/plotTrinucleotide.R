#' Plot Trinucleotide SBS Data
#'
#' This function processes and plots trinucleotide data.
#' It first applies specified filters and transformations
#' to the data and then generates a visual representation
#' of the results.
#' 
#' The function handles data normalization, exclusion, and retention
#' based on provided column names, and it creates detailed plots with
#' options for customization of plot aesthetics.
#'
#' @import ggplot2
#' @import ggpattern
#' @import patchwork
#' @importFrom ggplot2 ggplot geom_bar aes_string scale_fill_manual theme
#' @importFrom ggpattern geom_bar_pattern
#' @importFrom dplyr select
#' @importFrom magrittr '%>%'
#'
#' @param trinuc_df DataFrame with trinucleotide data.
#' @param exclude_if_type_present Vector specifying mutation locus readpair
#' overlap types (CO_MUT, SO_MUT, CO_REF, SO_REF, DO, SO_OTHER, CO_OTHER)
#' whose non-zero presence triggers exclusion of loci.
#' For instance, using 'c("DO")' will exclude any loci
#' that contain even a single discordant read-pair overlap (DO).
#' @param retain_if_type_present Vector specifying mutation locus readpair
#' overlap types (CO_MUT, SO_MUT, CO_REF, SO_REF, DO, SO_OTHER, CO_OTHER)
#' whose non-zero presence is necessary to retain those mutation locus.
#' For instance, using 'c("CO")' will retain any loci
#' that contain even a single concordant read-pair overlap (CO_MUT).
#' @param remove_type Vector specifying mutation locus readpair overlap types 
#' (e.g., SO_MUT, CO_REF, DO) to set to 0 across all loci in the dataframe.
#' @param normalize_counts Logical; if TRUE, normalize SBS counts to sum to 1.
#' @param show_overlap_type Logical; if TRUE, show read-pair overlap types.
#' @param ylim Numeric; limits of the y-axis for the plot.
#' @param plot_title String; title for the plots.
#' @param y_axis_title String; title for the y-axis.
#' @param draw_x_axis_labels Logical; whether to draw x-axis labels.
#' @param draw_y_axis_labels Logical; whether to draw y-axis labels.
#' @param draw_y_axis_title Logical; whether to display a title for the y-axis.
#' @param output_file String; name and path of output pdf file.
#' @param ggsave_params A list of parameters to be passed to ggplot2::ggsave().
#'   This list can include any of the arguments that 'ggsave()' accepts.
#'   Default settings of 'ggsave()' will be used unless specified.
#'   Example: 'list(width = 10, height = 8, dpi = 300)'
#'
#' @return A trinucleotide SBS plot object and an optional pdf file.
#' @export
#' @examples
#' \dontrun{
#'  plotTrinucleotide(trinuc_df)
#' }
#'
plotTrinucleotide <- function(
  trinuc_df,
  exclude_if_type_present = NULL,
  retain_if_type_present = NULL,
  remove_type = NULL,
  normalize_counts = TRUE,
  show_overlap_type = TRUE,
  ylim = c(0, 0.5),
  plot_title = "Trinucleotide Profile",
  y_axis_title = "Fraction of Single Base Substitutions",
  draw_x_axis_labels = TRUE,
  draw_y_axis_labels = TRUE,
  draw_y_axis_title = TRUE,
  output_file = "./trinucleotide_profile.pdf",
  ggsave_params = list(width = 17,
                      height = 6,
                      units = "cm",
                      device = "pdf")) {

  # Process trinucleotide data
  count_df_processed <- processTrinucleotideData(
    trinuc_df = trinuc_df,
    exclude_if_type_present = exclude_if_type_present,
    retain_if_type_present = retain_if_type_present,
    remove_type = remove_type,
    normalize_counts = normalize_counts
  )

  # Plot the processed data
  plotTrinucData(
    count_df = count_df_processed,
    ylim = ylim,
    show_overlap_type = show_overlap_type,
    plot_title = plot_title,
    y_axis_title = y_axis_title,
    draw_x_axis_labels = draw_x_axis_labels,
    draw_y_axis_labels = draw_y_axis_labels,
    draw_y_axis_title = draw_y_axis_title,
    output_file = output_file,
    ggsave_params = ggsave_params
  )
}

