context("Test call functions")
library(cfDNAPro)

path = examplePath("groups_picard")
mode_df <- callMode(path = path)
size_df <- callSize(path = path)
metric_df <- callMetrics(path = path)
valley_df <- callValleyDistance(path = path)
peak_df <- callPeakDistance(path = path)

expect_true(is.character(path))
expect_true(is.data.frame(mode_df))
expect_true(is.data.frame(size_df))
expect_true(is.data.frame(metric_df))
expect_true(is.data.frame(valley_df))
expect_true(is.data.frame(peak_df))





