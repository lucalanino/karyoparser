strip_bands <- function(band_str) {
  stringr::str_replace_all(band_str, "([pq]\\d+)\\.\\d+", "\\1")
}

normalize_token <- function(x) {
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "(^\\?|\\?$|~)", "")
  x <- strip_bands(x)
  out <- x
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

build_sample_meta <- function(input_df, min_metaphases = 2) {
  counts_raw <- input_df |>
    dplyr::mutate(
      bracket = stringr::str_extract_all(original_karyotype, "\\[[^\\]]+\\]")
    ) |>
    tidyr::unnest_longer(bracket, keep_empty = TRUE) |>
    dplyr::mutate(
      content = ifelse(
        is.na(bracket),
        NA_character_,
        stringr::str_replace_all(bracket, "\\[|\\]", "")
      ),
      is_cp = !is.na(content) & stringr::str_detect(content, "^cp\\d+$"),
      is_range_count = !is.na(content) & stringr::str_detect(content, "~"),
      max_count = dplyr::case_when(
        is_cp ~ suppressWarnings(as.numeric(stringr::str_extract(
          content,
          "\\d+"
        ))),
        is_range_count ~ {
          nums <- stringr::str_extract_all(content, "\\d+")
          purrr::map_dbl(
            nums,
            ~ if (length(.x) >= 1) {
              suppressWarnings(max(as.numeric(.x)))
            } else {
              NA_real_
            }
          )
        },
        TRUE ~ suppressWarnings(as.numeric(content))
      )
    )

  counts_tbl <- counts_raw |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::summarise(
      total_metaphases = if (all(is.na(bracket))) {
        NA_integer_
      } else {
        as.integer(sum(tidyr::replace_na(max_count, 0)))
      },
      is_composite = any(is_cp, na.rm = TRUE),
      .groups = "drop"
    )

  input_df |>
    dplyr::mutate(
      head = stringr::str_extract(original_karyotype, "^\\d+(?:~\\d+)?"),
      is_range_karyotype = !is.na(head) & stringr::str_detect(head, "~"),
      cleaned = stringr::str_trim(stringr::str_replace_all(
        original_karyotype,
        "\\[[^\\]]+\\]",
        ""
      )),
      normal_karyotype = ifelse(
        !is_range_karyotype & cleaned %in% c("46,XX", "46,XY"),
        1L,
        0L
      ),
      chromosome_count = chromosome_count_from_karyotype(
        original_karyotype,
        min_metaphases = min_metaphases
      )
    ) |>
    dplyr::left_join(counts_tbl, by = ".pk_row_id")
}

build_clone_tokens <- function(sample_meta) {
  clones_tbl <- sample_meta |>
    dplyr::select(
      .pk_row_id,
      original_karyotype,
      is_composite,
      is_range_karyotype
    ) |>
    dplyr::mutate(
      clone_str = stringr::str_replace_all(
        original_karyotype,
        "\\[[^\\]]+\\]",
        ""
      )
    ) |>
    tidyr::separate_longer_delim(clone_str, delim = "/") |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::mutate(clone_id = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(clone_str = stringr::str_trim(clone_str)) |>
    tidyr::separate(
      clone_str,
      into = c("chromosome_count", "aberrations"),
      sep = ",",
      fill = "right",
      extra = "merge",
      remove = FALSE
    ) |>
    dplyr::mutate(
      chromosome_count = stringr::str_trim(chromosome_count),
      aberrations = stringr::str_trim(aberrations)
    )

  tokens_tbl <- clones_tbl |>
    dplyr::mutate(aberrations = ifelse(is.na(aberrations), "", aberrations)) |>
    tidyr::separate_longer_delim(aberrations, delim = ",") |>
    dplyr::mutate(aberr_raw = stringr::str_trim(aberrations)) |>
    dplyr::select(
      .pk_row_id,
      clone_id,
      chromosome_count,
      aberr_raw,
      is_range_karyotype,
      is_composite
    ) |>
    dplyr::group_by(.pk_row_id, clone_id) |>
    dplyr::mutate(token_index = dplyr::row_number()) |>
    dplyr::ungroup()

  range_tokens <- sample_meta |>
    dplyr::filter(is_range_karyotype) |>
    dplyr::transmute(
      .pk_row_id,
      clone_id = NA_integer_,
      chromosome_count = NA_character_,
      aberr_raw = "chromosome_count_range",
      token_index = NA_integer_,
      is_range_karyotype = TRUE,
      is_composite = is_composite
    )

  dplyr::bind_rows(tokens_tbl, range_tokens) |>
    dplyr::mutate(aberr_norm = normalize_token(aberr_raw))
}

# Rules fire independently; mutual exclusivity, where wanted, is encoded in
# the regex itself.
match_rules <- function(tokens_tbl, rules, rule_flag_names) {
  match_text <- strip_bands(tokens_tbl$aberr_raw)
  ids <- tokens_tbl$.pk_row_id
  out <- tibble::tibble(.pk_row_id = unique(ids))
  for (nm in .sanitize_flag_name(rule_flag_names)) {
    out[[nm]] <- 0L
  }
  for (j in seq_len(nrow(rules))) {
    hit <- stringr::str_detect(match_text, rules$regex[j])
    hit[is.na(hit)] <- FALSE
    nm <- .sanitize_flag_name(rules$flag_name[j])
    out[[nm]] <- pmax(
      out[[nm]],
      as.integer(out$.pk_row_id %in% unique(ids[hit]))
    )
  }
  out
}
