pk <- function(...) parse_karyo(..., on_issues = "warn", verbose = FALSE)

has_issue <- function(karyotypes, issue_type) {
  a <- karyoparser:::.assess_karyotypes(karyotypes)
  as.integer(
    seq_along(karyotypes) %in%
      a$reported_issues$row_index[a$reported_issues$issue_type == issue_type]
  )
}
