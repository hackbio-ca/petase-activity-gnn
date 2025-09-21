# Manual file inspection to debug the gitter file reading issue

gitter_path <- "Plate Images/unwashed/BHET25_2d_6_1.JPG.dat"

cat("=== Manual File Inspection ===\n")

if (!file.exists(gitter_path)) {
  stop("File not found: ", gitter_path)
}

# Read all lines
all_lines <- readLines(gitter_path)
cat("Total lines in file:", length(all_lines), "\n\n")

# Show first 10 lines with line numbers
cat("First 10 lines:\n")
for (i in seq_len(min(10, length(all_lines)))) {
  cat(sprintf("Line %2d: %s\n", i, all_lines[i]))
}

cat("\n=== Finding Header ===\n")

# Find comment lines
comment_lines <- which(startsWith(all_lines, "#"))
cat("Comment lines:", paste(comment_lines, collapse = ", "), "\n")

# Find non-comment lines
non_comment_lines <- which(!startsWith(all_lines, "#"))
cat("Non-comment lines start at:", non_comment_lines[seq_len(min(5, length(non_comment_lines)))], "\n")

# Look at the first few non-comment lines
cat("\nFirst few non-comment lines:\n")
for (i in seq_len(min(3, length(non_comment_lines)))) {
  line_num <- non_comment_lines[i]
  cat(sprintf("Line %2d: %s\n", line_num, all_lines[line_num]))
  
  # Parse this line
  cols <- strsplit(all_lines[line_num], "\t")[[1]]
  cat(sprintf("  -> %d columns: %s\n", length(cols), paste(cols, collapse = " | ")))
}

# Look for lines containing "xl", "xr", "yt", "yb"
cat("\n=== Looking for header pattern ===\n")
for (i in seq_along(all_lines)) {
  line <- all_lines[i]
  if (grepl("xl", line) && grepl("xr", line) && grepl("yt", line) && grepl("yb", line)) {
    cat(sprintf("Found header pattern at line %d: %s\n", i, line))
    
    # Parse columns
    cols <- strsplit(line, "\t")[[1]]
    cols <- trimws(cols)
    cat("Columns:", paste(cols, collapse = ", "), "\n")
    
    # Check for required columns
    required <- c("xl", "xr", "yt", "yb")
    found <- required %in% cols
    cat("Required columns found:", paste(required, "=", found, collapse = ", "), "\n")
    break
  }
}