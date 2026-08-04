compute_aneuploidy <- function(tokens_tbl, chroms) {
  mono_cols <- paste0("mono", chroms)
  tris_cols <- paste0("tris", chroms)

  aneuploidy_tbl <- tokens_tbl |>
    dplyr::filter(aberr_raw != "") |>
    dplyr::transmute(
      .pk_row_id,
      sign = dplyr::case_when(
        stringr::str_detect(aberr_raw, "^-") ~ "-",
        stringr::str_detect(aberr_raw, "^\\+") ~ "+",
        TRUE ~ ""
      ),
      chrom = stringr::str_extract(aberr_raw, "(?<=^[-+])[0-9XY]+"),
      is_mono_tri = sign %in% c("-", "+") & !is.na(chrom),
      prefix = dplyr::case_when(
        sign == "-" ~ "mono",
        sign == "+" ~ "tris",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(is_mono_tri) |>
    dplyr::mutate(flag = paste0(prefix, chrom), value = 1L) |>
    dplyr::distinct(.pk_row_id, flag, .keep_all = TRUE) |>
    dplyr::select(.pk_row_id, flag, value) |>
    tidyr::pivot_wider(names_from = flag, values_from = value, values_fill = 0L)

  for (nm in c(mono_cols, tris_cols)) {
    if (!nm %in% names(aneuploidy_tbl)) aneuploidy_tbl[[nm]] <- 0L
  }
  aneuploidy_tbl
}

# Universal structural flags: fire independently of any rule match on the same
# token. Built around the indicator tokens in .aberr_indicators_paren /
# .aberr_indicators_bare (karyoparser-package.R), the single source of truth
# for these strings, rather than retyping them here; each flag adds only the
# anchoring/boundary it needs on top of its shared token.
.general_flag_patterns <- c(
  general_dicentric = paste0(.aberr_indicators_paren[["dic"]], "\\("),
  general_isodicentric = paste0(.aberr_indicators_paren[["idic"]], "\\("),
  general_isochromosome = paste0(.aberr_indicators_paren[["i"]], "\\("),
  general_pseudodicentric = paste0(.aberr_indicators_paren[["psu_dic"]], "\\("),
  general_ring = paste0("\\b", .aberr_indicators_paren[["r"]], "\\("),
  general_insertion = paste0(.aberr_indicators_paren[["ins"]], "\\("),
  general_duplication = paste0(.aberr_indicators_paren[["dup"]], "\\("),
  general_triplication = paste0(.aberr_indicators_paren[["trp"]], "\\("),
  general_translocation = paste0(
    "^[?~]?",
    .aberr_indicators_paren[["t"]],
    "\\([0-9XY]"
  ),
  general_addition = paste0(.aberr_indicators_paren[["add"]], "\\("),
  general_inversion = paste0(.aberr_indicators_paren[["inv"]], "\\("),
  general_deletion = paste0(.aberr_indicators_paren[["del"]], "\\("),
  general_marker = paste0("\\b", .aberr_indicators_bare[["mar"]], "\\b"),
  general_derivative = paste0(
    "^\\+?i?",
    .aberr_indicators_paren[["der"]],
    "\\("
  )
)

# Structural flags counted by the monosomal-karyotype rule (excludes a lone marker).
.general_flags_for_monosomal <- setdiff(
  names(.general_flag_patterns),
  "general_marker"
)

compute_general_flags <- function(tokens_tbl) {
  match_text <- strip_bands(tokens_tbl$aberr_raw)
  ids <- tokens_tbl$.pk_row_id
  out <- tibble::tibble(.pk_row_id = unique(ids))
  for (nm in names(.general_flag_patterns)) {
    hit <- stringr::str_detect(match_text, .general_flag_patterns[[nm]])
    hit[is.na(hit)] <- FALSE
    out[[nm]] <- as.integer(out$.pk_row_id %in% unique(ids[hit]))
  }
  out
}

compute_comma_counts <- function(tokens_tbl) {
  sex_first_regex <- paste0(
    "^(",
    paste(.sex_complements, collapse = "|"),
    ")$"
  )

  clone_counts <- tokens_tbl |>
    dplyr::filter(!is.na(clone_id)) |>
    dplyr::group_by(.pk_row_id, clone_id) |>
    dplyr::summarise(
      tokens_nonempty = sum(aberr_raw != ""),
      first_token = dplyr::first(aberr_raw[aberr_raw != ""]),
      first_is_sex = !is.na(first_token) &
        stringr::str_detect(first_token, sex_first_regex),
      has_idem = any(aberr_raw == "idem", na.rm = TRUE),
      raw_comma_count = as.integer(
        tokens_nonempty - ifelse(first_is_sex, 1L, 0L)
      ),
      .groups = "drop"
    )

  stemline_counts <- clone_counts |>
    dplyr::filter(clone_id == 1) |>
    dplyr::select(.pk_row_id, stemline_aberr_count = raw_comma_count)

  clone_counts <- clone_counts |>
    dplyr::left_join(stemline_counts, by = ".pk_row_id")

  clone_counts |>
    dplyr::mutate(
      comma_count_aberrations = dplyr::case_when(
        has_idem &
          clone_id > 1 &
          !is.na(stemline_aberr_count) ~ raw_comma_count +
          stemline_aberr_count -
          1L,
        TRUE ~ raw_comma_count
      )
    ) |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::summarise(
      comma_count_aberrations = max(tidyr::replace_na(
        comma_count_aberrations,
        0L
      )),
      .groups = "drop"
    )
}

compute_unique_counts <- function(tokens_tbl) {
  tokens_tbl |>
    dplyr::filter(
      !(aberr_norm %in%
        c("", "chromosome_count_range", "idem", "sl", .sex_complements))
    ) |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::summarise(
      n_unique_aberr = dplyr::n_distinct(aberr_norm),
      .groups = "drop"
    )
}
