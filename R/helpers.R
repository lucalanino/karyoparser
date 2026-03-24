truncate_str <- function(s, max_len = 40) {
  ifelse(nchar(s) > max_len, paste0(substr(s, 1, max_len - 3), "..."), s)
}

strip_bands <- function(band_str) {
  stringr::str_replace_all(band_str, "([pq]\\d+)\\.\\d+", "\\1")
}

normalize_token <- function(x) {
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "(^\\?|\\?$|~)", "")
  x <- strip_bands(x)
  out <- x
  # vectorized classification
  is_marker <- tidyr::replace_na(stringr::str_detect(x, "^[+-]?\\d*mar"), FALSE)
  out[is_marker] <- "marker_chromosomes"

  is_der <- tidyr::replace_na(stringr::str_detect(x, "^\\+?i?der\\("), FALSE)
  if (any(is_der)) {
    out[is_der] <- paste0(
      "der_",
      stringr::str_extract(x[is_der], "i?der\\([^)]+\\)")
    )
  }

  is_dic <- tidyr::replace_na(
    stringr::str_detect(x, "^(dic|idic|psu dic)\\("),
    FALSE
  )
  if (any(is_dic)) {
    head <- stringr::str_extract(x[is_dic], "^(dic|idic|psu dic)")
    paren <- stringr::str_extract(x[is_dic], "\\([^)]+\\)")
    out[is_dic] <- paste0(head, "_", paren)
  }
  out
}
