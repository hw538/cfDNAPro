#' Title
#'
#' @param vmat 
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
getscore <- function(vmat, x) {
  cols = which(as.numeric(colnames(vmat)) %in% c(x[1]:x[2]))
  if (length(cols) > 1) {
    res = rowSums(vmat[, cols])
  } else {
    res = vmat[, cols]
  }
  rm(vmat)
  gc()
  return(res)
}


#' callMIDS - helps extract MIDS signals for a given cfDNA sequence Grange object in an given region of interest
#'
#' @param fragment_obj 
#' @param regions 
#' @param xlimits 
#' @param fraglen 
#'
#' @return
#' @export
#'
#' @examples
callMIDS <- function(fragment_obj, 
  regions,
  xlimits = c(-1000, 1000),
  fraglen = list(c(0, 150), c(151, 200), c(200, 1000))) {
  
  message("Started to extract midpoint distribution ...")
  
  # check fragment_obj and regions.
  
  if (class(fragment_obj)[1] != "GRanges" | class(regions)[1] != "GRanges") {
    
    message("Please provide a valid grange object for fragment_obj and regions ...")   
    return()
  } else {
    
    overlap = findOverlaps(fragment_obj, regions)
    if (length(overlap) > 0) {
      vplotfrag = VplotR::plotVmat(fragment_obj, regions, xlim = xlimits, return_Vmat = TRUE)
      vplotfragbf = VplotR::plotVmat(
        fragment_obj,
        flank(regions, 2000, start = TRUE),
        xlim = xlimits,
        return_Vmat = TRUE
      )
      vplotfragaf = VplotR::plotVmat(
        fragment_obj,
        flank(regions, 2000, start = FALSE),
        xlim = xlimits,
        return_Vmat = TRUE
      )
      bgcorr = vplotfrag - ((vplotfragbf + vplotfragaf) / 2)
      final = do.call(rbind.data.frame, lapply(fraglen, getscore, vmat = bgcorr))
      colnames(final) = rownames(bgcorr)
    } else {
      final = data.frame(matrix(0, length(fraglen), 1999))
    }
    windowscored=zoo::rollapply(t(final),width=10, FUN=function(x) {mean(x)},by=10)
    windowscored=zoo::rollapply(windowscored,width=3, FUN=function(x) {mean(x)})
    windowscored=data.frame(t(windowscored))
    xlims = xlimits[[1]]:xlimits[[2]]
    xlimval = xlims[c(-1,-length(xlims))]
    colnames(windowscored)[1:197]=round(zoo::rollapply(zoo::rollapply(xlimval,width=10, FUN=function(x) {mean(x)},by=10),width=3,FUN=function(x) {mean(x)}))
    windowscored$fraglen=as.character(paste0("from",lapply(fraglen, `[[`, 1),"to",lapply(fraglen, `[[`, 2)))
    windowscored_spread = tidyr::gather(windowscored, key = "name", value = "value", -fraglen) %>%
      tidyr::unite(col="name.fraglen", name, fraglen,sep = ".") %>% 
      tidyr::spread(key = name.fraglen, value=value)
    windowscored_spread = (windowscored_spread/(as.numeric(length(regions))*as.numeric(length(fragment_obj))))*1000000
    return(windowscored_spread)
  }
}
