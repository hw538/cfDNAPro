#' getscore - helps to collect vplot score for a given fragment length range
#'
#' @param vmat - a vplot mat data frame from vplotR package with fragment length in rows and coverage signals in the columns
#' @param x - a vector with a fragment length range
#'
#' @return data frame with coverage signals cumulated for the given fragment length ranges
#' @export
#'
#' @examples
getscore <- function(vmat, x) {

  # collecting rows of vplot mat falling with in a given fragment length range
  
  cols <- which(as.numeric(colnames(vmat)) %in% c(x[1]:x[2]))

  # if given fragment length range has a signal then it is cumulated together
  if (length(cols) > 1) {
    res <- rowSums(vmat[, cols])
  } else {
    res <- vmat[, cols]
  }
  return(res)
}


#' @import VplotR
#' @import tidyr
#' @import zoo
#'
#' callMIDS - helps extract MIDS signals for a given cfDNA sequence Grange object in an given region of interest
#'
#' @param fragment_obj GRange object of the cfDNA sequence data
#' @param regions GRange object of the regions of interest
#' @param regnm string with a name for the regions of interest
#' @param xlimits gives the number of bp to look for signals from the midpoint of region of interest; deafults to -1000 to 1000 bp
#' @param fraglen gives the list of fragment length ranges to cumulate the signals; defaults to short 0 to 150 bp, medium 151 to 200 bp and long 201 to 1000 bp
#'
#' @return tibble with MIDS signals from region of interest, window scored and normalized
#' @export
#'
#' @examples
callMIDS <- function(fragment_obj, 
  regions,
  regnm,
  xlimits = c(-1000, 1000),
  fraglen = list(c(0, 150), c(151, 200), c(201, 1000))) {
  
  message("Started to extract midpoint distribution ...")
  
  # check fragment_obj and regions.
  
  if (class(fragment_obj)[1] != "GRanges" | class(regions)[1] != "GRanges") {
    
    message("Please provide a valid grange object for fragment_obj and regions ...")   
    return()
  } else {
    
    overlap <- findOverlaps(fragment_obj, regions)
    
    # check if there is overlap between fragment_obj and regions.
    
    if (length(overlap) > 0) {

      # collect cummulative vplot information for given fragment_obj in given regions.
      
      vplotfrag <- VplotR::plotVmat(fragment_obj, regions, xlim = xlimits, return_Vmat = TRUE)

      # collect cummulative vplot information for given fragment_obj aroound given regions for background correction.
      
      vplotfragbf <- VplotR::plotVmat(
        fragment_obj,
        flank(regions, 2000, start = TRUE),
        xlim = xlimits,
        return_Vmat = TRUE
      )
      vplotfragaf <- VplotR::plotVmat(
        fragment_obj,
        flank(regions, 2000, start = FALSE),
        xlim = xlimits,
        return_Vmat = TRUE
      )

      # local background correction.
      
      bgcorr <- vplotfrag - ((vplotfragbf + vplotfragaf) / 2)

      # cumulation of vplot scores for fragment length ranges
      
      final <- do.call(rbind.data.frame, lapply(fraglen, getscore, vmat = bgcorr))
      colnames(final) <- rownames(bgcorr)
      
    } else {
      
      # matrix with zeros if no overlap found
      
      final <- data.frame(matrix(0, length(fraglen), length(xlimits[1]:xlimits[2])-2))
    }

    # windowscoring is done to reduce the -1000 to 1000 bp signals to 197 points

    # this windowscoring widths have been choosen based on performance in this project and can be changed accordingly for future projects
        
    windowscored <- zoo::rollapply(t(final),width=10, FUN=function(x) {mean(x)},by=10)
    windowscored <- zoo::rollapply(windowscored,width=3, FUN=function(x) {mean(x)})
    windowscored <- data.frame(t(windowscored))
    xlims <- xlimits[[1]]:xlimits[[2]]
    xlimval <- xlims[c(-1, -length(xlims))]
    colnames(windowscored) <- round(zoo::rollapply(
      zoo::rollapply(
        xlimval,
        width = 10,
        FUN = function(x) {
          mean(x)
        },
        by = 10
      ),
      width = 3,
      FUN = function(x) {
        mean(x)
      }
    ))
    windowscored$fraglen <- as.character(paste0("from", lapply(fraglen, `[[`, 1),
      "to", lapply(fraglen, `[[`, 2)))
    windowscored$name <- rep(regnm,length(fraglen))
    
    # spread the data with the name and fraglen in column name so the data can be combined across samples
    # this is done for easy combining of data from all samples for different regions of interest
    
    windowscored_spread <- tidyr::gather(windowscored, key = "name", 
      value = "value", -fraglen) %>%
      tidyr::unite(col = "name.fraglen", name, fraglen, sep = ".") %>%
      tidyr::spread(key = name.fraglen, value = value)

    # normalization of the reads per million and regions per million
    
    RPM <- as.numeric(length(fragment_obj)))/1000000
    RegPM <- as.numeric(length(regions))/1000000
    windowscored_spread <- windowscored_spread / (RPM * RegPM)

    return(windowscored_spread)
  }
}
