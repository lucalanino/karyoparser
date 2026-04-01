#' Classify Ploidy from Chromosome Count
#'
#' @param n Integer chromosome count (or NA).
#' @return Character string: one of "near_haploid", "low_hypodiploid",
#'   "high_hypodiploid", "diploid", "hyperdiploid", "other", or "unknown".
#' @keywords internal
ploidy_from_count <- function(n) {
  if (is.na(n)) {
    return("unknown")
  }
  if (n >= 23 && n <= 29) {
    "near_haploid"
  } else if (n >= 30 && n <= 33) {
    "low_hypodiploid"
  } else if (n >= 40 && n <= 45) {
    "high_hypodiploid"
  } else if (n == 46) {
    "diploid"
  } else if (n > 50) {
    "hyperdiploid"
  } else {
    "other"
  }
}

extract_clone_data <- function(karyotype) {
  # Split into clones before removing brackets
  clones <- stringr::str_split(karyotype, "/")[[1]]
  lapply(clones, function(cl) {
    # Extract metaphase count from bracket
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
    # Extract chromosome count
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

#' Classify Ploidy for a Karyotype String
#'
#' For composite karyotypes, filters to clones with >= 5 metaphases (when
#' brackets exist), then uses the most abnormal clone (furthest from 46).
#'
#' @param karyotype Single karyotype string.
#' @return A list with elements: ploidy (character), mixed (logical),
#'   chromosome_count (integer).
#' @keywords internal
ploidy_category <- function(karyotype) {
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
    return(list(
      ploidy = "unknown",
      mixed = FALSE,
      chromosome_count = NA_integer_
    ))
  }
  # Most abnormal = furthest from 46
  distances <- abs(eligible_counts - 46)
  most_abnormal <- eligible_counts[which.max(distances)]
  ploidy <- ploidy_from_count(most_abnormal)
  # Mixed if different ploidy categories across eligible clones
  all_ploidies <- vapply(eligible_counts, ploidy_from_count, character(1))
  mixed <- length(unique(all_ploidies)) > 1
  list(
    ploidy = ploidy,
    mixed = mixed,
    chromosome_count = as.integer(most_abnormal)
  )
}
