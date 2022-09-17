

#' Multiple gsub
#'
#' Run gsub iteratively
#' @param pattern Character string containing a regular expression (or
#'          character string for ‘fixed = TRUE’) to be matched in the
#'          given character vector.  Coerced by ‘as.character’ to a
#'          character string if possible.
#' @param replacement a replacement for matched pattern.
#'          Coerced to character if possible.  For ‘fixed = FALSE’ this
#'          can include backreferences ‘"\1"’ to ‘"\9"’ to parenthesized
#'          subexpressions of ‘pattern’.  For ‘perl = TRUE’ only, it can
#'          also contain ‘"\U"’ or ‘"\L"’ to convert the rest of the
#'          replacement to upper or lower case and ‘"\E"’ to end case
#'          conversion.
#' @param x a character vector where matches are sought, or an object
#'          which can be coerced by ‘as.character’ to a character vector.
#'          Long vectors are supported.
#' @param ignore.case if ‘FALSE’, the pattern matching is _case sensitive_
#'          and if ‘TRUE’, case is ignored during matching.
#' @param perl logical.  Should Perl-compatible regexps be used?
#' @param fixed logical.  If ‘TRUE’, ‘pattern’ is a string to be matched as
#'          is.  Overrides all conflicting arguments.
#' @param useBytes logical.  If ‘TRUE’ the matching is done byte-by-byte rather
#'          than character-by-character.  See ‘Details’.
#' @return Substituted text
#' @export
msub <- function(pattern, replacement, x, ignore.case=FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE) {
    lgt1i <- function(item, i) {
        if (length(item) > 1) {
            item[i]
        } else {
            item
        }
    }
    for (i in seq_along(pattern)) {
        x = gsub(pattern[i],
                 replacement[i],
                 x,
                 ignore.case=lgt1i(ignore.case, i),
                 perl=lgt1i(perl, i),
                 fixed=lgt1i(fixed, i),
                 useBytes=lgt1i(useBytes, i))
    }
    return(x)
}
