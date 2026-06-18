extract_clone_data <- function(karyotype) {
  # Split into clones before removing brackets
  clones <- stringr::str_split(karyotype, "/")[[1]]
  lapply(clones, function(cl) {
    bracket <- stringr::str_extract(cl, "\\[[^\\]]+\\]")
    metaphases <- NA_real_
    if (!is.na(bracket)) {
      content <- stringr::str_replace_all(bracket, "\\[|\\]", "")
      content <- stringr::str_replace(content, "^cp", "")
      if (stringr::str_detect(content, "~")) {
        nums <- as.numeric(stringr::str_extract_all(content, "\\d+")[[1]])
        if (length(nums) >= 2) metaphases <- max(nums)
      } else {
        metaphases <- suppressWarnings(as.numeric(content))
      }
    }
    clean <- stringr::str_replace_all(cl, "\\[[^\\]]+\\]", "")
    raw_count <- stringr::str_extract(clean, "^\\d+(?:~\\d+)?")
    chrom_count <- NA_real_
    if (!is.na(raw_count) && raw_count != "") {
      if (stringr::str_detect(raw_count, "~")) {
        nums <- as.numeric(stringr::str_extract_all(raw_count, "\\d+")[[1]])
        if (length(nums) >= 2) chrom_count <- round(mean(nums))
      } else {
        chrom_count <- as.numeric(raw_count)
      }
    }
    list(chrom_count = chrom_count, metaphases = metaphases)
  })
}

# Chromosome count for a single karyotype string, taken from the most abnormal
# eligible clone (the one whose count is furthest from 46).
chromosome_count_from_karyotype <- function(karyotype) {
  clone_data <- extract_clone_data(karyotype)
  chrom_counts <- vapply(clone_data, `[[`, NA_real_, "chrom_count")
  metaphases <- vapply(clone_data, `[[`, NA_real_, "metaphases")
  # If any clone has brackets, only consider clones with >= 5 metaphases
  has_any_brackets <- any(!is.na(metaphases))
  if (has_any_brackets) {
    eligible <- !is.na(chrom_counts) & !is.na(metaphases) & metaphases >= 5
  } else {
    eligible <- !is.na(chrom_counts)
  }
  eligible_counts <- chrom_counts[eligible]
  if (length(eligible_counts) == 0) {
    return(NA_integer_)
  }
  # Most abnormal = furthest from 46
  distances <- abs(eligible_counts - 46)
  as.integer(eligible_counts[which.max(distances)])
}
