# All chromosome arms (autosomes 1-22 plus X/Y, each p and q) for which
# unbalanced-derived partial losses get their own `unbal_partial_loss_<arm>`
# column.
.unbal_partial_loss_arms <- paste0(
  rep(c(as.character(1:22), "X", "Y"), each = 2),
  c("p", "q")
)

# Canonicalize a der's t(...) partner list: sort partners by chromosome (keeping
# breakpoints aligned) and build a signature so reciprocal ders of the SAME
# translocation group together. Two-partner ders also expose p1/p2/bb1/bb2 for
# loss derivation; multi-partner (three-way+) ders leave those NA, since loss is
# only derived for the simple two-partner case.
.canonical_der_t <- function(chroms_str, bands_str) {
  chroms <- strsplit(chroms_str, ";", fixed = TRUE)[[1]]
  bands <- if (is.na(bands_str)) {
    rep(NA_character_, length(chroms))
  } else {
    b <- strsplit(bands_str, ";", fixed = TRUE)[[1]]
    if (length(b) == length(chroms)) b else rep(NA_character_, length(chroms))
  }
  ord <- order(vapply(chroms, .chrom_sort_key, integer(1)))
  chroms <- chroms[ord]
  bands <- bands[ord]
  n <- length(chroms)
  sig <- paste0(
    "t(",
    paste(chroms, collapse = ";"),
    ")(",
    paste(tidyr::replace_na(bands, ""), collapse = ";"),
    ")"
  )
  two <- n == 2L
  list(
    n = n,
    sig = sig,
    chroms = chroms,
    p1 = if (two) chroms[[1]] else NA_character_,
    p2 = if (two) chroms[[2]] else NA_character_,
    bb1 = if (two) bands[[1]] else NA_character_,
    bb2 = if (two) bands[[2]] else NA_character_
  )
}

# TRUE when the donor multiset covers the partner multiset -- i.e. every partner
# chromosome appears as a centromere donor at least as many times as it occurs
# in the translocation. This is the complete reciprocal set (balanced); a lone
# der leaves some partner un-donored (unbalanced). Generalizes the two-partner
# rule (heterologous: both partners donored; homologous t(a;a): two der(a)) to
# any number of partners.
.donors_cover_partners <- function(donors, partners) {
  all(vapply(
    unique(partners),
    function(chrom) sum(donors == chrom) >= sum(partners == chrom),
    logical(1)
  ))
}

# Extract der(...)t(...) tokens with a canonical translocation signature and a
# per-(row, clone, signature) `balanced_pair` flag (TRUE when the reciprocal set
# is complete -- every partner appears as a centromere donor). Handles any number
# of translocation partners (two-way and three-way+). Shared by
# classify_translocation_balance() and derive_unbalanced_loss().
.der_translocations <- function(tokens_tbl) {
  der_t <- tokens_tbl |>
    dplyr::filter(!is.na(clone_id), aberr_raw != "") |>
    dplyr::mutate(
      raw = strip_bands(aberr_raw),
      is_der = stringr::str_detect(raw, "^\\+?i?der\\("),
      der_inner = stringr::str_match(raw, "^\\+?i?der\\(([^)]+)\\)")[, 2],
      t_chroms = stringr::str_match(
        raw,
        "t\\(([0-9XY]+(?:;[0-9XY]+)+)\\)"
      )[, 2],
      t_bands = stringr::str_match(
        raw,
        "t\\([0-9XY]+(?:;[0-9XY]+)+\\)\\(([^)]+)\\)"
      )[, 2]
    ) |>
    dplyr::filter(
      is_der,
      !is.na(t_chroms),
      !is.na(der_inner),
      !stringr::str_detect(der_inner, ";")
    )

  if (nrow(der_t) == 0) {
    return(
      der_t |>
        dplyr::mutate(
          n_partners = integer(0),
          p1 = character(0),
          p2 = character(0),
          bb1 = character(0),
          bb2 = character(0),
          donor = character(0),
          t_sig = character(0),
          balanced_pair = logical(0)
        )
    )
  }

  canon <- purrr::pmap(
    list(der_t$t_chroms, der_t$t_bands),
    .canonical_der_t
  )

  der_t |>
    dplyr::mutate(
      donor = der_inner,
      n_partners = purrr::map_int(canon, "n"),
      t_sig = purrr::map_chr(canon, "sig"),
      p1 = purrr::map_chr(canon, "p1"),
      p2 = purrr::map_chr(canon, "p2"),
      bb1 = purrr::map_chr(canon, "bb1"),
      bb2 = purrr::map_chr(canon, "bb2"),
      partners = purrr::map(canon, "chroms")
    ) |>
    dplyr::group_by(.pk_row_id, clone_id, t_sig) |>
    dplyr::mutate(
      balanced_pair = .donors_cover_partners(donor, partners[[1]])
    ) |>
    dplyr::ungroup()
}

# Classify each row as carrying balanced and/or unbalanced translocations.
#
# Balanced  := a bare t(...) token, OR a complete der()t() reciprocal set where
#              every partner chromosome appears as a centromere donor within a
#              clone (two-way pair or three-way+ set).
# Unbalanced := a der()t() whose reciprocal set is incomplete -- e.g. a lone
#              der(a)t(a;b) or a lone der(a)t(a;b;c) -- OR a whole-arm der(a;b)
#              (single derivative). A row may be both (mixed).
#
# These flags are independent of `derivative_chromosome`, which still fires.
classify_translocation_balance <- function(
  tokens_tbl,
  der_t = .der_translocations(tokens_tbl)
) {
  toks <- tokens_tbl |>
    dplyr::filter(!is.na(clone_id), aberr_raw != "") |>
    dplyr::mutate(
      raw = strip_bands(aberr_raw),
      is_der = stringr::str_detect(raw, "^\\+?i?der\\("),
      der_inner = stringr::str_match(raw, "^\\+?i?der\\(([^)]+)\\)")[, 2],
      is_bare_t = stringr::str_detect(raw, "^[?~]?t\\([0-9XY]")
    )

  empty <- tibble::tibble(
    .pk_row_id = integer(),
    balanced_translocation = integer(),
    unbalanced_translocation = integer()
  )

  bare_rows <- toks$.pk_row_id[toks$is_bare_t]

  # Whole-arm der(a;b): a single derivative, always unbalanced.
  whole_arm_rows <- toks$.pk_row_id[
    toks$is_der &
      !is.na(toks$der_inner) &
      stringr::str_detect(toks$der_inner, ";")
  ]

  der_balanced_rows <- der_t$.pk_row_id[der_t$balanced_pair]
  der_unbalanced_rows <- der_t$.pk_row_id[!der_t$balanced_pair]

  balanced_ids <- unique(c(bare_rows, der_balanced_rows))
  unbalanced_ids <- unique(c(der_unbalanced_rows, whole_arm_rows))
  all_ids <- sort(unique(c(balanced_ids, unbalanced_ids)))

  if (length(all_ids) == 0) {
    return(empty)
  }

  tibble::tibble(
    .pk_row_id = all_ids,
    balanced_translocation = as.integer(all_ids %in% balanced_ids),
    unbalanced_translocation = as.integer(all_ids %in% unbalanced_ids)
  )
}

# Derive implied partial losses from lone (unbalanced) der(a)t(a;b)(bp_a;bp_b):
#   - chromosome a loses material distal to bp_a -> arm(bp_a)
#   - partner b loses its centromere-side material  -> opposite arm of bp_b
# Emitted into dedicated `unbal_partial_loss_<arm>` columns (one per chromosome arm,
# autosomes plus X/Y) plus a generic `unbal_partial_loss` flag. These are kept
# SEPARATE from del()/mono* signals -- an unbalanced-derived 5q loss is not
# conflated with a true del(5q).
#
# Scope: only SIMPLE single-junction ders are resolved -- exactly one two-partner
# t() whose centromere donor is one of the partners. Multi-junction chains
# (der with >1 t(), three-way t(a;b;c), or a der named after a non-partner
# chromosome) cannot be resolved by this two-partner model, so no loss is
# derived for them (they are still flagged unbalanced by the classifier).
# Copy-number-aware gains are out of scope (they need whole-karyotype reasoning).
derive_unbalanced_loss <- function(der_t) {
  cols <- paste0("unbal_partial_loss_", .unbal_partial_loss_arms)
  empty <- tibble::tibble(.pk_row_id = integer())
  for (nm in c(cols, "unbal_partial_loss")) {
    empty[[nm]] <- integer()
  }

  dt <- der_t |>
    dplyr::filter(
      !balanced_pair,
      n_partners == 2L,
      stringr::str_count(raw, "t\\(") == 1L,
      donor == p1 | donor == p2
    )
  if (nrow(dt) == 0) {
    return(empty)
  }

  segs <- dt |>
    dplyr::mutate(
      donor_is_p1 = donor == p1,
      bp_a = dplyr::if_else(donor_is_p1, bb1, bb2),
      bp_b = dplyr::if_else(donor_is_p1, bb2, bb1),
      b_chrom = dplyr::if_else(donor_is_p1, p2, p1),
      a_seg = .arm_segment(donor, .arm_of(bp_a)),
      b_seg = .arm_segment(b_chrom, .opp_arm(.arm_of(bp_b)))
    )

  long <- dplyr::bind_rows(
    segs |> dplyr::transmute(.pk_row_id, seg = a_seg),
    segs |> dplyr::transmute(.pk_row_id, seg = b_seg)
  ) |>
    dplyr::filter(!is.na(seg))

  if (nrow(long) == 0) {
    return(empty)
  }

  any_loss <- long |>
    dplyr::distinct(.pk_row_id) |>
    dplyr::mutate(unbal_partial_loss = 1L)

  named <- long |>
    dplyr::filter(seg %in% .unbal_partial_loss_arms) |>
    dplyr::distinct(.pk_row_id, seg) |>
    dplyr::mutate(col = paste0("unbal_partial_loss_", seg), value = 1L) |>
    dplyr::select(.pk_row_id, col, value) |>
    tidyr::pivot_wider(
      names_from = col,
      values_from = value,
      values_fill = 0L
    )

  out <- any_loss |> dplyr::left_join(named, by = ".pk_row_id")
  for (nm in cols) {
    if (!nm %in% names(out)) out[[nm]] <- 0L
  }
  out |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(cols),
      ~ tidyr::replace_na(., 0L)
    ))
}

# Numeric sort key for a chromosome label (X -> 23, Y -> 24).
.chrom_sort_key <- function(chrom) {
  dplyr::case_when(
    chrom == "X" ~ 23L,
    chrom == "Y" ~ 24L,
    TRUE ~ suppressWarnings(as.integer(chrom))
  )
}

# First arm letter (p/q) of a band string, or NA.
.arm_of <- function(band) {
  stringr::str_extract(band, "^[pq]")
}

# Opposite chromosome arm.
.opp_arm <- function(arm) {
  dplyr::case_when(arm == "p" ~ "q", arm == "q" ~ "p", TRUE ~ NA_character_)
}

# Combine chromosome + arm into a segment key (e.g. "5q"), NA if arm unknown.
.arm_segment <- function(chrom, arm) {
  dplyr::if_else(is.na(arm) | is.na(chrom), NA_character_, paste0(chrom, arm))
}
