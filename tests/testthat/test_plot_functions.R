context("Test plot functions")
library(cfDNAPro)

path = examplePath("groups_picard")
mode_df <- callMode(path = path)
size_df <- callSize(path = path)
metric_df <- callMetrics(path = path)
valley_df <- callValleyDistance(path = path)
peak_df <- callPeakDistance(path = path)

all_to_one_plot <- plotAllToOne(size_df)
metrics_plot <- plotMetrics(metric_df)
mode_plot <- plotMode(mode_df)
modesummary_plot <- plotModeSummary(mode_df)
peakdist_plot <- plotPeakDistance(peak_df)
valleydist_plot <- plotValleyDistance(valley_df)
singlegroup_plot <- plotSingleGroup(metric_df)


test_that("Test that plot functions are no problem",{
    expect_s3_class(all_to_one_plot$prop_plot, class = "ggplot")
    expect_s3_class(metrics_plot$median_prop_plot, class = "ggplot")
    expect_s3_class(mode_plot, class = "ggplot")
    expect_s3_class(modesummary_plot, class = "ggplot")
    expect_s3_class(peakdist_plot, class = "ggplot")
    expect_s3_class(valleydist_plot, class = "ggplot")
    expect_s3_class(singlegroup_plot$prop_plot, class = "ggplot")
    
})

#> Test passed 😀