# Chromosome count per karyotype string, taken from the most abnormal
# eligible clone (the one whose count is furthest from 46) in each row. NA
# where no clone in a row is eligible. Vectorized: flattens all rows' clones
# into one (row, clone) pair per element so every regex runs once over the
# whole dataset, then picks the per-row argmax via order()+duplicated()
# instead of looping row by row (mirrors the flatten-and-group approach used
# in validate_karyotypes()).
chromosome_count_from_karyotype <- function(karyotypes) {
  n <- length(karyotypes)
  clones_list <- stringr::str_split(karyotypes, "/")
  flat_row <- rep(seq_len(n), lengths(clones_list))
  clones <- unlist(clones_list)

  bracket <- stringr::str_extract(clones, "\\[[^\\]]+\\]")
  has_bracket <- !is.na(bracket)
  content <- stringr::str_replace_all(bracket, "\\[|\\]", "")
  content <- stringr::str_replace(content, "^cp", "")
  is_range_meta <- has_bracket & stringr::str_detect(content, "~")
  range_max <- purrr::map_dbl(
    stringr::str_extract_all(content, "\\d+"),
    \(x) if (length(x) >= 2) max(as.numeric(x)) else NA_real_
  )
  metaphases <- dplyr::case_when(
    is_range_meta ~ range_max,
    has_bracket ~ suppressWarnings(as.numeric(content)),
    TRUE ~ NA_real_
  )

  clean <- stringr::str_replace_all(clones, "\\[[^\\]]+\\]", "")
  raw_count <- stringr::str_extract(clean, "^\\d+(?:~\\d+)?")
  is_range_count <- !is.na(raw_count) & stringr::str_detect(raw_count, "~")
  range_mean <- purrr::map_dbl(
    stringr::str_extract_all(raw_count, "\\d+"),
    \(x) if (length(x) >= 2) round(mean(as.numeric(x))) else NA_real_
  )
  chrom_count <- dplyr::case_when(
    is_range_count ~ range_mean,
    !is.na(raw_count) & raw_count != "" ~ suppressWarnings(as.numeric(
      raw_count
    )),
    TRUE ~ NA_real_
  )

  has_any_brackets <- .agg_any_by_group(has_bracket, flat_row, n)[flat_row]
  eligible <- ifelse(
    has_any_brackets,
    !is.na(chrom_count) & !is.na(metaphases) & metaphases >= 5,
    !is.na(chrom_count)
  )

  # Most abnormal = furthest from 46; ineligible clones are pinned to -Inf so
  # they never win the per-row argmax below unless every clone in a row is
  # ineligible, in which case the picked distance stays -Inf and the row's
  # result is NA.
  distance <- ifelse(eligible, abs(chrom_count - 46), -Inf)
  ord <- order(flat_row, -distance)
  picked <- !duplicated(flat_row[ord])
  picked_count <- chrom_count[ord][picked]
  picked_distance <- distance[ord][picked]

  result <- ifelse(is.infinite(picked_distance), NA_real_, picked_count)
  as.integer(result)
}
