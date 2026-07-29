# 02_Probe_Excel_Rows.R
library(openxlsx)

file_path <- "Drought fractional counts.xlsx"
cat("--- PEAKING AT FIRST 10 ROWS OF RAW EXCEL STRUCTURE ---\n")

for(i in 1:10) {
  row_data <- read.xlsx(file_path, sheet = 1, rows = i, colNames = FALSE)
  if(nrow(row_data) > 0) {
    cat("Row", i, "Excerpt:", paste(row_data[1, 1:min(5, ncol(row_data))], collapse=" | "), "\n")
  } else {
    cat("Row", i, "is completely EMPTY\n")
  }
}
