#' Validate and Create a Custom Rules Table
#'
#' Validates a data frame against the required schema for karyotype
#' parsing rules and returns it as a `karyo_rules` object accepted by
#' [parse_karyo()].
#'
#' @param rules A data frame with columns: `flag_name` (character), `regex`
#'   (character), `category` (character). Each rule fires independently when its
#'   `regex` matches a token; there is no priority/competition between rules, so
#'   mutual exclusivity (where wanted) must be encoded in the `regex` itself.
#'   `flag_name` must be unique -- each maps to exactly one output column, so
#'   alternative patterns for the same flag must be combined into a single
#'   `regex` (e.g. with `|`) rather than given as separate rows. `flag_name` can
#'   use ISCN-style punctuation (e.g. `"t(9;22)(q34;q11)"`) for readability --
#'   the output column name is a sanitized version, with runs of
#'   non-alphanumeric characters replaced by `_` (e.g. `t_9_22_q34_q11`).
#'   `flag_name` values that sanitize to the same column name are rejected.
#'
#' @return A `karyo_rules` object (tibble subclass) accepted by [parse_karyo()].
#' @export
validate_rules <- function(rules) {
  if (!is.data.frame(rules)) {
    stop("`rules` must be a data frame.", call. = FALSE)
  }
  required <- c("flag_name", "regex", "category")
  missing_cols <- setdiff(required, names(rules))
  if (length(missing_cols) > 0L) {
    stop(
      "`rules` missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(rules) == 0L) {
    stop("`rules` must have at least one row.", call. = FALSE)
  }
  dupes <- unique(rules$flag_name[duplicated(rules$flag_name)])
  if (length(dupes) > 0L) {
    stop(
      "`rules` has duplicate `flag_name` value(s): ",
      paste(dupes, collapse = ", "),
      ". Each flag_name must map to exactly one rule -- combine alternative ",
      "patterns into a single regex instead.",
      call. = FALSE
    )
  }
  sanitized <- .sanitize_flag_name(rules$flag_name)
  collides <- duplicated(sanitized) | duplicated(sanitized, fromLast = TRUE)
  if (any(collides)) {
    stop(
      "`rules` has `flag_name` value(s) that collide once sanitized into an ",
      "output column name: ",
      paste(unique(rules$flag_name[collides]), collapse = ", "),
      ". Rename one of the conflicting flag_name values so they remain ",
      "distinct after stripping non-alphanumeric characters.",
      call. = FALSE
    )
  }
  structure(rules, class = c("karyo_rules", class(rules)))
}

# Output column name for a rule's flag_name: non-alphanumeric runs -> "_",
# trimmed of leading/trailing "_" (e.g. "t(9;22)(q34;q11)" -> "t_9_22_q34_q11").
.sanitize_flag_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

#' Default Parsing Rules for Myeloid Neoplasms
#'
#' A `karyo_rules` tibble of regex-based rules for myeloid-specific recurrent
#' lesions: specific translocations/inversions, variable-partner translocations,
#' and chromosome-arm-specific deletions/additions. Passed to [parse_karyo()]
#' via the `rules` argument by default. Use [validate_rules()] to build a custom
#' rule set with the same schema.
#'
#' General, disease-agnostic structural aberrations (dicentrics, rings,
#' insertions, generic translocations/deletions, derivatives, etc.) are NOT
#' rules: they are detected universally by the parser and surfaced as the
#' `general_*` output columns, independently of the rule set in use and of any
#' specific rule firing on the same token.
#'
#' @format A `karyo_rules` tibble with columns: `flag_name`, `regex`,
#'   `category`.
#' @seealso [validate_rules()]
#' @export
myeloid_rules <- validate_rules(tibble::tribble(
  ~flag_name          , ~regex                                                                                        , ~category             ,
  "t(15;17)(q24;q21)" , "t\\(15;17\\)\\((q24|q22);q21\\)|t\\(17;15\\)\\(q21;(q24|q22)\\)"                             , "specific_tx"         ,
  "t(8;21)(q22;q22)"  , "t\\(8;21\\)\\((q22|q21(\\.3)?);q22\\)|t\\(21;8\\)\\(q22;(q22|q21(\\.3)?)\\)"                 , "specific_tx"         ,
  "inv(16)(p13q22)"   , "inv\\(16\\)\\(p13q22\\)"                                                                     , "specific_tx"         ,
  "t(16;16)(p13;q22)" , "t\\(16;16\\)\\(p13;q22\\)"                                                                   , "specific_tx"         ,
  "t(9;11)(p21;q23)"  , "t\\(9;11\\)\\(p21;q23\\)|t\\(11;9\\)\\(q23;p21\\)"                                           , "specific_tx"         ,
  "t(6;9)(p22;q34)"   , "t\\(6;9\\)\\(p22;q34\\)|t\\(9;6\\)\\(q34;p22\\)"                                             , "specific_tx"         ,
  "inv(3)(q21q26)"    , "inv\\(3\\)\\(q21q26\\)"                                                                      , "specific_tx"         ,
  "t(3;3)(q21;q26)"   , "t\\(3;3\\)\\(q21;q26\\)"                                                                     , "specific_tx"         ,
  "t(9;22)(q34;q11)"  , "t\\(9;22\\)\\(q34;q11\\)|t\\(22;9\\)\\(q11;q34\\)"                                           , "specific_tx"         ,
  "t(1;3)(p36;q21)"   , "t\\(1;3\\)\\(p36;q21\\)|t\\(3;1\\)\\(q21;p36\\)"                                             , "specific_tx"         ,
  "t(1;22)(p13;q13)"  , "t\\(1;22\\)\\(p13;q13\\)|t\\(22;1\\)\\(q13;p13\\)"                                           , "specific_tx"         ,
  "t(3;5)(q25;q35)"   , "t\\(3;5\\)\\(q25;q35\\)|t\\(5;3\\)\\(q35;q25\\)"                                             , "specific_tx"         ,
  "t(5;11)(q35;p15)"  , "t\\(5;11\\)\\(q35;p15\\)|t\\(11;5\\)\\(p15;q35\\)"                                           , "specific_tx"         ,
  "t(7;12)(q36;p13)"  , "t\\(7;12\\)\\(q36;p13\\)|t\\(12;7\\)\\(p13;q36\\)"                                           , "specific_tx"         ,
  "t(8;16)(p11;p13)"  , "t\\(8;16\\)\\(p11;p13\\)|t\\(16;8\\)\\(p13;p11\\)"                                           , "specific_tx"         ,
  "t(10;11)(p12;q14)" , "t\\(10;11\\)\\(p12;q14\\)|t\\(11;10\\)\\(q14;p12\\)"                                         , "specific_tx"         ,
  "t(11;12)(p15;p13)" , "t\\(11;12\\)\\(p15;p13\\)|t\\(12;11\\)\\(p13;p15\\)"                                         , "specific_tx"         ,
  "t(16;21)(p11;q22)" , "t\\(16;21\\)\\(p11;q22\\)|t\\(21;16\\)\\(q22;p11\\)"                                         , "specific_tx"         ,
  "t(16;21)(q24;q22)" , "t\\(16;21\\)\\(q24;q22\\)|t\\(21;16\\)\\(q22;q24\\)"                                         , "specific_tx"         ,
  "inv(16)(p13q24)"   , "inv\\(16\\)\\(p13q24\\)"                                                                     , "specific_tx"         ,
  "t(6;9)_other"      , "t\\(6;9\\)\\((?!p22;q34\\))p[^;)]*;q[^)]*\\)|t\\(9;6\\)\\((?!q34;p22\\))q[^;)]*;p[^)]*\\)"   , "specific_tx"         ,
  "t(9;11)_other"     , "t\\(9;11\\)\\((?!p21;q23\\))p[^;)]*;q[^)]*\\)|t\\(11;9\\)\\((?!q23;p21\\))q[^;)]*;p[^)]*\\)" , "specific_tx"         ,
  "t(9;22)_other"     , "t\\(9;22\\)\\((?!q34;q11\\))q[^;)]*;q[^)]*\\)|t\\(22;9\\)\\((?!q11;q34\\))q[^;)]*;q[^)]*\\)" , "specific_tx"         ,
  "inv(3)_other"      , "inv\\(3\\)\\((?!q21q26\\))q\\d+q\\d+\\)"                                                     , "specific_tx"         ,
  "t(3;3)_other"      , "t\\(3;3\\)\\((?!q21;q26\\))q[^;)]*;q[^)]*\\)"                                                , "specific_tx"         ,
  "t(v;11p15)"        , "t\\([0-9XY]+;11\\)\\([^)]+;p15\\)|t\\(11;[0-9XY]+\\)\\(p15;[^)]+\\)"                         , "variable_partner"    ,
  "t(v;11q23)"        , "t\\([0-9XY]+;11\\)\\([^)]+;q23\\)|t\\(11;[0-9XY]+\\)\\(q23;[^)]+\\)"                         , "variable_partner"    ,
  "t(3q26;v)"         , "t\\(3;[0-9XY]+\\)\\(q26;[^)]+\\)|t\\([0-9XY]+;3\\)\\([^)]+;q26\\)"                           , "variable_partner"    ,
  "del(5q)"           , "del\\(5q|del\\(5\\)\\(q"                                                                     , "chromosome_specific" ,
  "t(5q)"             , "t\\(5;[^)]+\\)\\(q[^;)]*;[^)]*\\)|t\\([^;]+;5\\)\\([^;]*;q[^)]*\\)"                          , "chromosome_specific" ,
  "add(5q)"           , "add\\(5q|add\\(5\\)\\(q"                                                                     , "chromosome_specific" ,
  "del(7q)"           , "del\\(7q|del\\(7\\)\\(q"                                                                     , "chromosome_specific" ,
  "del(12p)"          , "del\\(12p|del\\(12\\)\\(p"                                                                   , "chromosome_specific" ,
  "t(12p)"            , "t\\(12;[^)]+\\)\\(p[^;)]*;[^)]*\\)|t\\([^;]+;12\\)\\([^;]*;p[^)]*\\)"                        , "chromosome_specific" ,
  "add(12p)"          , "add\\(12p|add\\(12\\)\\(p"                                                                   , "chromosome_specific" ,
  "del(13q)"          , "del\\(13q|del\\(13\\)\\(q"                                                                   , "chromosome_specific" ,
  "i(17q)"            , "i\\(17q|i\\(17\\)\\(q"                                                                       , "chromosome_specific" ,
  "add(17p)"          , "add\\(17p|add\\(17\\)\\(p"                                                                   , "chromosome_specific" ,
  "del(17p)"          , "del\\(17p|del\\(17\\)\\(p"                                                                   , "chromosome_specific" ,
  "del(20q)"          , "del\\(20q|del\\(20\\)\\(q"                                                                   , "chromosome_specific" ,
  "del(11q)"          , "del\\(11q|del\\(11\\)\\(q"                                                                   , "chromosome_specific" ,
  "idic(X)(q13)"      , "idic\\(X\\)\\(q13"                                                                           , "chromosome_specific"
))
