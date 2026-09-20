#' Synthetic Example Karyotypes
#'
#' A synthetic dataset of 101 ISCN karyotype strings for trying out
#' [check_karyo()], [preprocess_karyo()], and [parse_karyo()]. Covers normal
#' karyotypes, one example per `myeloid_rules` flag, general structural
#' aberrations with no matching specific rule, aneuploidy-only cases,
#' balanced/unbalanced translocations, complex and monosomal karyotypes,
#' multi-clone/chimeric/ploidy variety, and a range of messy or out-of-scope
#' strings (leading dots, HTML entities, FISH suffixes, `mos`/`ncSCA`/
#' constitutional markers, ...) to demonstrate `check_karyo()` and
#' `preprocess_karyo()`.
#'
#' @format A tibble with 101 rows and columns:
#' \describe{
#'   \item{sample_id}{character. Synthetic sample identifier (`"EX001"`, ...).}
#'   \item{karyotype}{character. Raw ISCN karyotype string.}
#' }
#' @seealso [parse_karyo()], [check_karyo()], [preprocess_karyo()]
#' @export
example_karyotypes <- tibble::tribble(
  ~sample_id , ~karyotype                                                 ,

  # Normal karyotypes
  "EX001"    , "46,XX[20]"                                                ,
  "EX002"    , "46,XY[20]"                                                ,
  "EX003"    , "46,XY"                                                    ,

  # One example per myeloid_rules flag (specific translocations/inversions)
  "EX004"    , "46,XX,t(15;17)(q24;q21)[20]"                              ,
  "EX005"    , "46,XY,t(8;21)(q22;q22)[18]"                               ,
  "EX006"    , "46,XX,inv(16)(p13q22)[15]"                                ,
  "EX007"    , "46,XY,t(16;16)(p13;q22)[12]"                              ,
  "EX008"    , "46,XX,t(9;11)(p21;q23)[14]"                               ,
  "EX009"    , "46,XY,t(6;9)(p22;q34)[16]"                                ,
  "EX010"    , "46,XX,inv(3)(q21q26)[10]"                                 ,
  "EX011"    , "46,XY,t(3;3)(q21;q26)[11]"                                ,
  "EX012"    , "46,XX,t(9;22)(q34;q11)[20]"                               ,
  "EX013"    , "46,XY,t(1;3)(p36;q21)[13]"                                ,
  "EX014"    , "46,XX,t(1;22)(p13;q13)[9]"                                ,
  "EX015"    , "46,XY,t(3;5)(q25;q35)[10]"                                ,
  "EX016"    , "46,XX,t(5;11)(q35;p15)[12]"                               ,
  "EX017"    , "46,XY,t(7;12)(q36;p13)[15]"                               ,
  "EX018"    , "46,XX,t(8;16)(p11;p13)[10]"                               ,
  "EX019"    , "46,XY,t(10;11)(p12;q14)[8]"                               ,
  "EX020"    , "46,XX,t(11;12)(p15;p13)[9]"                               ,
  "EX021"    , "46,XY,t(16;21)(p11;q22)[10]"                              ,
  "EX022"    , "46,XX,t(16;21)(q24;q22)[11]"                              ,
  "EX023"    , "46,XY,inv(16)(p13q24)[12]"                                ,

  # "_other" non-canonical breakpoint variants
  "EX024"    , "46,XX,t(6;9)(p23;q33)[10]"                                ,
  "EX025"    , "46,XY,t(9;11)(p22;q24)[9]"                                ,
  "EX026"    , "46,XX,t(9;22)(q34;q12)[10]"                               ,
  "EX027"    , "46,XY,inv(3)(q22q27)[8]"                                  ,
  "EX028"    , "46,XX,t(3;3)(q22;q27)[9]"                                 ,

  # Variable-partner families
  "EX029"    , "46,XY,t(4;11)(q21;p15)[10]"                               ,
  "EX030"    , "46,XX,t(2;11)(p21;q23)[11]"                               ,
  "EX031"    , "46,XY,t(3;21)(q26;q11)[10]"                               ,

  # Chromosome-arm-specific deletions/additions/isochromosome
  "EX032"    , "46,XX,del(3)(q21q26)[11]"                                 ,
  "EX033"    , "46,XX,del(5)(q13q33)[15]"                                 ,
  "EX034"    , "46,XY,t(5;7)(q31;p15)[10]"                                ,
  "EX035"    , "46,XX,add(5)(q13)[9]"                                     ,
  "EX036"    , "46,XY,del(7)(q22q36)[12]"                                 ,
  "EX037"    , "46,XX,del(12)(p11p13)[10]"                                ,
  "EX038"    , "46,XY,t(12;9)(p13;q34)[9]"                                ,
  "EX039"    , "46,XX,add(12)(p11)[8]"                                    ,
  "EX040"    , "46,XY,del(13)(q14q22)[10]"                                ,
  "EX041"    , "46,XX,i(17)(q10)[15]"                                     ,
  "EX042"    , "46,XY,add(17)(p11)[10]"                                   ,
  "EX043"    , "46,XX,del(17)(p11p13)[9]"                                 ,
  "EX044"    , "46,XY,del(20)(q11q13)[10]"                                ,
  "EX045"    , "46,XX,del(11)(q14q23)[8]"                                 ,
  "EX046"    , "46,X,idic(X)(q13)[10]"                                    ,

  # General structural flags with no matching specific rule
  "EX047"    , "46,XY,dic(9;22)(p11;q11)[10]"                             ,
  "EX048"    , "46,XX,psu dic(15;17)(q10;q10)[10]"                        ,
  "EX049"    , "46,XY,r(7)[10]"                                           ,
  "EX050"    , "46,XX,ins(2;3)(p21;q13q26)[10]"                           ,
  "EX051"    , "46,XY,dup(1)(q21q32)[10]"                                 ,
  "EX052"    , "46,XX,trp(8)(q22q24)[10]"                                 ,
  "EX053"    , "45,XY,-7,+mar[10]"                                        ,
  "EX054"    , "46,XX,i(9)(q10)[10]"                                      ,
  "EX055"    , "46,XY,inv(9)(p11q13)[10]"                                 ,
  "EX056"    , "46,XX,idic(15)(q10)[10]"                                  ,

  # Aneuploidy only
  "EX057"    , "47,XY,+8[15]"                                             ,
  "EX058"    , "45,XY,-7[15]"                                             ,
  "EX059"    , "45,XX,-5[12]"                                             ,
  "EX060"    , "47,XX,+21[18]"                                            ,
  "EX061"    , "47,XY,+9[10]"                                             ,
  "EX062"    , "45,XX,-X[10]"                                             ,
  "EX063"    , "47,XX,+Y[8]"                                              ,
  "EX064"    , "45,XY,-18[9]"                                             ,

  # Balanced translocations (reciprocal der pairs)
  "EX065"    , "46,XX,der(5)t(5;17)(q11;q11),der(17)t(5;17)(q11;q11)[10]" ,
  "EX066"    , "46,XY,der(1)t(1;19)(q23;p13),der(19)t(1;19)(q23;p13)[10]" ,

  # Unbalanced translocations / derived partial loss
  "EX067"    , "45,XX,-7,der(8)t(8;21)(q22;q22)[10]"                      ,
  "EX068"    , "46,XX,der(5)t(5;8)(q11;q11)[10]"                          ,
  "EX069"    , "46,XY,der(9)t(9;22)(q34;q11)[10]"                         ,
  "EX070"    , "46,XX,der(1)t(2;3)(p11;q22)[10]"                          ,
  "EX071"    , "45,XX,der(9)t(9;22)(q34;q11)[10]"                         ,

  # Complex karyotypes (>= 3 distinct aberrations)
  "EX072"    , "46,XX,del(5)(q13),del(7)(q22),add(11)(p11)[10]"           ,
  "EX073"    , "46,XY,t(9;22)(q34;q11),del(20)(q11q13),+8[12]"            ,
  "EX074"    , "45,XX,-7,del(5)(q13),add(21)(p11)[9]"                     ,

  # Monosomal karyotypes
  "EX075"    , "45,XX,-5,-7[10]"                                          ,
  "EX076"    , "45,XY,-7,del(5)(q13)[10]"                                 ,
  "EX077"    , "43,XX,-5,-7,-17[8]"                                       ,

  # Messy/dirty strings -- demonstrate check_karyo() / preprocess_karyo()
  "EX078"    , ".46,XX,del(5)(q13)[10]"                                   ,
  "EX079"    , "46, XX,+8[10]"                                            ,
  "EX080"    , "46,XX,PSU DIC(15;22)[10]"                                 ,
  "EX081"    , ",46,,XX,del(7)(q22)[10],"                                 ,
  "EX082"    , "47,XY,\uFF0B8[10]"                                        ,
  "EX083"    , "46,XX[cp 20]"                                             ,
  "EX084"    , "46,XX,t(9;22)(q34;q11)[20] .Clinical note here"           ,
  "EX085"    , "46,XX,t(9;22)(q34;q11)[15] &lt;AML&gt;"                   ,
  "EX086"    , ".//46,XX,+8[10]"                                          ,
  "EX087"    , "46,XX,t(9;22)[15]/46,XX[3] .nuc ish(PDGFRA x3)[20/200]"   ,
  "EX088"    , "46,XX der(7)t(7;12)(q36;q24)[10]"                         ,
  "EX089"    , "46,XX,del(5)(q13),del(7)(q22)[10]/46,idem,+8[5]"          ,

  # Out-of-scope rows -- flagged by check_karyo(); mos/constitutional/triple-//
  # are unfixable, ncSCA is fixable (token stripped, remaining clone parsed)
  "EX090"    , "46,XX//47//48"                                            ,
  "EX091"    , "mos 47,XXY[10]/46,XY[5]"                                  ,
  "EX092"    , "ncSCA[4]/46,XY[11]"                                       ,
  "EX093"    , "47,XXYc[20]"                                              ,
  "EX094"    , "46,XX &amp; notes"                                        ,

  # Multi-clone / chimeric / ploidy variety
  "EX095"    , "46,XX[10]/46,idem[5]"                                     ,
  "EX096"    , "46,XX,t(9;22)(q34;q11)[15]//46,XX[5]"                     ,
  "EX097"    , "69,XXX[10]"                                               ,
  "EX098"    , "45,X[10]"                                                 ,
  "EX099"    , "47,XY,+21[10]/46,XY[10]"                                  ,
  "EX100"    , "46,XX,dup(1)(q21q32),t(9;22)(q34;q11)[10]"                ,
  "EX101"    , "46,XY,der(9)t(9;22;11)(q34;q11;q23)[10]"
)
