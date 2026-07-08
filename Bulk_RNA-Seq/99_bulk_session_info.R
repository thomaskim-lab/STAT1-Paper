############################################################
### 99_bulk_session_info.R
############################################################

source("00_bulk_setup.R")

if (!requireNamespace("sessioninfo", quietly = TRUE)) {
  install.packages("sessioninfo")
}
library(sessioninfo)

sink(file.path(session_info_dir, "bulk_rnaseq_sessionInfo.txt"))

cat("Bulk RNA-seq session information\n")
cat("Generated on: ", as.character(Sys.time()), "\n\n")

cat("Base R sessionInfo()\n")
print(sessionInfo())

cat("\n\nsessioninfo::session_info()\n")
print(sessioninfo::session_info())

sink()
