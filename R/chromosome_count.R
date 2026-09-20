# Chromosome count per row, from the eligible clone furthest from 46. A clone
# is eligible at >= `min_metaphases`; rows where no clone carries a bracket
# skip the threshold entirely. NA when no clone is eligible.
chromosome_count_from_karyotype <- function(karyotypes, min_metaphases = 2) {
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
    !is.na(chrom_count) & !is.na(metaphases) & metaphases >= min_metaphases,
    !is.na(chrom_count)
  )

  # Ineligible clones are pinned to -Inf so they never win the argmax below;
  # a row where all clones are ineligible stays -Inf and resolves to NA.
  distance <- ifelse(eligible, abs(chrom_count - 46), -Inf)
  ord <- order(flat_row, -distance)
  picked <- !duplicated(flat_row[ord])
  picked_count <- chrom_count[ord][picked]
  picked_distance <- distance[ord][picked]

  result <- ifelse(is.infinite(picked_distance), NA_real_, picked_count)
  as.integer(result)
}
