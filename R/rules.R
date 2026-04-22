#' Validate and Create a Custom Rules Table
#'
#' Validates a data frame against the required schema for karyotype parsing rules
#' and returns it as a `karyo_rules` object accepted by [parse_karyo()].
#'
#' @param rules A data frame with columns: `flag_name` (character), `regex`
#'   (character), `category` (character), `priority` (positive numeric),
#'   `counts_for_monosomal` (logical), `competition_group` (character).
#'
#' @return A `karyo_rules` object (tibble subclass) accepted by [parse_karyo()].
#' @export
validate_rules <- function(rules) {
  if (!is.data.frame(rules)) {
    stop("`rules` must be a data frame.", call. = FALSE)
  }
  required <- c(
    "flag_name",
    "regex",
    "category",
    "priority",
    "counts_for_monosomal",
    "competition_group"
  )
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
  if (!is.numeric(rules$priority) && !is.integer(rules$priority)) {
    stop("`priority` must be numeric.", call. = FALSE)
  }
  if (any(rules$priority < 1)) {
    stop("All `priority` values must be >= 1.", call. = FALSE)
  }
  if (!is.logical(rules$counts_for_monosomal)) {
    stop("`counts_for_monosomal` must be logical.", call. = FALSE)
  }
  structure(rules, class = c("karyo_rules", class(rules)))
}

#' Default Parsing Rules for Myeloid Neoplasms
#'
#' A `karyo_rules` tibble of regex-based rules for identifying chromosomal
#' aberrations relevant to myeloid neoplasms. Passed to [parse_karyo()] via
#' the `rules` argument by default. Use [validate_rules()] to build a custom
#' rule set with the same schema.
#'
#' @format A `karyo_rules` tibble with columns: `flag_name`, `regex`,
#'   `category`, `priority`, `counts_for_monosomal`, `competition_group`.
#' @export
myeloid_rules <- validate_rules(tibble::tribble(
  ~flag_name              , ~regex                                                                        , ~category             , ~priority , ~counts_for_monosomal , ~competition_group ,
  # Specific translocations (both orientations)
  "t(15;17)(q24;q21)"     , "t\\(15;17\\)\\((q24|q22);q21\\)|t\\(17;15\\)\\(q21;(q24|q22)\\)"             , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(8;21)(q22;q22)"      , "t\\(8;21\\)\\((q22|q21(\\.3)?);q22\\)|t\\(21;8\\)\\(q22;(q22|q21(\\.3)?)\\)" , "specific_tx"         ,       100 , FALSE                 , "translocation"    ,
  "inv(16)(p13q22)"       , "inv\\(16\\)\\(p13q22\\)"                                                     , "specific_tx"         ,       100 , FALSE                 , "inversion"        ,
  "t(16;16)(p13;q22)"     , "t\\(16;16\\)\\(p13;q22\\)"                                                   , "specific_tx"         ,       100 , FALSE                 , "translocation"    ,
  "t(9;11)(p21;q23)"      , "t\\(9;11\\)\\(p21;q23\\)|t\\(11;9\\)\\(q23;p21\\)"                           , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(6;9)(p22;q34)"       , "t\\(6;9\\)\\(p22;q34\\)|t\\(9;6\\)\\(q34;p22\\)"                             , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "inv(3)(q21q26)"        , "inv\\(3\\)\\(q21q26\\)"                                                      , "specific_tx"         ,       100 , TRUE                  , "inversion"        ,
  "t(3;3)(q21;q26)"       , "t\\(3;3\\)\\(q21;q26\\)"                                                     , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(9;22)(q34;q11)"      , "t\\(9;22\\)\\(q34;q11\\)|t\\(22;9\\)\\(q11;q34\\)"                           , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(1;3)(p36;q21)"       , "t\\(1;3\\)\\(p36;q21\\)|t\\(3;1\\)\\(q21;p36\\)"                             , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(1;22)(p13;q13)"      , "t\\(1;22\\)\\(p13;q13\\)|t\\(22;1\\)\\(q13;p13\\)"                           , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(3;5)(q25;q35)"       , "t\\(3;5\\)\\(q25;q35\\)|t\\(5;3\\)\\(q35;q25\\)"                             , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(5;11)(q35;p15)"      , "t\\(5;11\\)\\(q35;p15\\)|t\\(11;5\\)\\(p15;q35\\)"                           , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(7;12)(q36;p13)"      , "t\\(7;12\\)\\(q36;p13\\)|t\\(12;7\\)\\(p13;q36\\)"                           , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(8;16)(p11;p13)"      , "t\\(8;16\\)\\(p11;p13\\)|t\\(16;8\\)\\(p13;p11\\)"                           , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(10;11)(p12;q14)"     , "t\\(10;11\\)\\(p12;q14\\)|t\\(11;10\\)\\(q14;p12\\)"                         , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(11;12)(p15;p13)"     , "t\\(11;12\\)\\(p15;p13\\)|t\\(12;11\\)\\(p13;p15\\)"                         , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(16;21)(p11;q22)"     , "t\\(16;21\\)\\(p11;q22\\)|t\\(21;16\\)\\(q22;p11\\)"                         , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "t(16;21)(q24;q22)"     , "t\\(16;21\\)\\(q24;q22\\)|t\\(21;16\\)\\(q22;q24\\)"                         , "specific_tx"         ,       100 , TRUE                  , "translocation"    ,
  "inv(16)(p13q24)"       , "inv\\(16\\)\\(p13q24\\)"                                                     , "specific_tx"         ,       100 , TRUE                  , "inversion"        ,
  # Specific translocations - variant bands (same chromosomes/arms, different sub-bands)
  "t(6;9)_other"          , "t\\(6;9\\)\\(p[^;)]*;q[^)]*\\)|t\\(9;6\\)\\(q[^;)]*;p[^)]*\\)"               , "specific_tx"         ,        95 , TRUE                  , "translocation"    ,
  "t(9;11)_other"         , "t\\(9;11\\)\\(p[^;)]*;q[^)]*\\)|t\\(11;9\\)\\(q[^;)]*;p[^)]*\\)"             , "specific_tx"         ,        95 , TRUE                  , "translocation"    ,
  "t(9;22)_other"         , "t\\(9;22\\)\\(q[^;)]*;q[^)]*\\)|t\\(22;9\\)\\(q[^;)]*;q[^)]*\\)"             , "specific_tx"         ,        95 , TRUE                  , "translocation"    ,
  "inv(3)_other"          , "inv\\(3\\)\\(q\\d+q\\d+\\)"                                                  , "specific_tx"         ,        95 , TRUE                  , "inversion"        ,
  "t(3;3)_other"          , "t\\(3;3\\)\\(q[^;)]*;q[^)]*\\)"                                              , "specific_tx"         ,        95 , TRUE                  , "translocation"    ,
  # Variable partner
  "t(v;11p15)"            , "t\\([0-9XY]+;11\\)\\([^)]+;p15\\)|t\\(11;[0-9XY]+\\)\\(p15;[^)]+\\)"         , "variable_partner"    ,        90 , TRUE                  , "translocation"    ,
  "t(v;11q23)"            , "t\\([0-9XY]+;11\\)\\([^)]+;q23\\)|t\\(11;[0-9XY]+\\)\\(q23;[^)]+\\)"         , "variable_partner"    ,        90 , TRUE                  , "translocation"    ,
  "t(3q26;v)"             , "t\\(3;[0-9XY]+\\)\\(q26;[^)]+\\)|t\\([0-9XY]+;3\\)\\([^)]+;q26\\)"           , "variable_partner"    ,        90 , TRUE                  , "translocation"    ,
  # Chromosome-specific
  "del(5q)"               , "del\\(5q|del\\(5\\)\\(q"                                                     , "chromosome_specific" ,        85 , TRUE                  , "deletion"         ,
  "t(5q)"                 , "t\\(5;[^)]+\\)\\(q[^;)]*;[^)]*\\)|t\\([^;]+;5\\)\\([^;]*;q[^)]*\\)"          , "chromosome_specific" ,        85 , TRUE                  , "translocation"    ,
  "add(5q)"               , "add\\(5q|add\\(5\\)\\(q"                                                     , "chromosome_specific" ,        85 , TRUE                  , "addition"         ,
  "del(7q)"               , "del\\(7q|del\\(7\\)\\(q"                                                     , "chromosome_specific" ,        85 , TRUE                  , "deletion"         ,
  "del(12p)"              , "del\\(12p|del\\(12\\)\\(p"                                                   , "chromosome_specific" ,        85 , TRUE                  , "deletion"         ,
  "t(12p)"                , "t\\(12;[^)]+\\)\\(p[^;)]*;[^)]*\\)|t\\([^;]+;12\\)\\([^;]*;p[^)]*\\)"        , "chromosome_specific" ,        85 , TRUE                  , "translocation"    ,
  "add(12p)"              , "add\\(12p|add\\(12\\)\\(p"                                                   , "chromosome_specific" ,        85 , TRUE                  , "addition"         ,
  "del(13q)"              , "del\\(13q|del\\(13\\)\\(q"                                                   , "chromosome_specific" ,        85 , TRUE                  , "deletion"         ,
  "i(17q)"                , "i\\(17q|i\\(17\\)\\(q"                                                       , "chromosome_specific" ,        85 , TRUE                  , "isochromosome"    ,
  "add(17p)"              , "add\\(17p|add\\(17\\)\\(p"                                                   , "chromosome_specific" ,        85 , TRUE                  , "addition"         ,
  "del(17p)"              , "del\\(17p|del\\(17\\)\\(p"                                                   , "chromosome_specific" ,        85 , TRUE                  , "deletion"         ,
  "del(20q)"              , "del\\(20q|del\\(20\\)\\(q"                                                   , "chromosome_specific" ,        85 , TRUE                  , "deletion"         ,
  "del(11q)"              , "del\\(11q|del\\(11\\)\\(q"                                                   , "chromosome_specific" ,        85 , TRUE                  , "deletion"         ,
  # Dicentric / idic / psu dic
  "idic(X)(q13)"          , "idic\\(X\\)\\(q13"                                                           , "chromosome_specific" ,        95 , TRUE                  , "isodicentric"     ,
  "dicentric"             , "dic\\("                                                                      , "general"             ,        80 , TRUE                  , "dicentric"        ,
  "isodicentric"          , "idic\\("                                                                     , "general"             ,        85 , TRUE                  , "isodicentric"     ,
  "pseudodicentric"       , "psu dic\\("                                                                  , "general"             ,        85 , TRUE                  , "pseudodicentric"  ,
  # General structural
  "ring_chromosome"       , "\\br\\("                                                                     , "general"             ,        70 , TRUE                  , "ring"             ,
  "insertion"             , "ins\\("                                                                      , "general"             ,        70 , TRUE                  , "insertion"        ,
  "duplication"           , "dup\\("                                                                      , "general"             ,        70 , TRUE                  , "duplication"      ,
  "triplication"          , "trp\\("                                                                      , "general"             ,        70 , TRUE                  , "triplication"     ,
  "general_translocation" , "^[?~]?t\\([0-9XY]"                                                           , "general"             ,        65 , TRUE                  , "translocation"    ,
  "general_addition"      , "add\\("                                                                      , "general"             ,        60 , TRUE                  , "addition"         ,
  "general_inversion"     , "inv\\("                                                                      , "general"             ,        60 , TRUE                  , "inversion"        ,
  "general_deletion"      , "del\\("                                                                      , "general"             ,        60 , TRUE                  , "deletion"         ,
  "marker_chromosome"     , "\\bmar\\b"                                                                   , "general"             ,        60 , FALSE                 , "marker"           ,
  "derivative_chromosome" , "^\\+?i?der\\("                                                               , "general"             ,        55 , TRUE                  , "derivative"
))
