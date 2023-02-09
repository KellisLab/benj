


#' Load RDS file from fname or re-compute and save
#'
#' @param fname RDS filename
#' @param func Function to compute RDS
#' @param expiration Date object used for timestamping old data
#' @export
load_cached <- function(fname, func, expiration="2000-01-01", ...) {
    lk = filelock::lock(paste0(fname, ".lock"))
    obj = NULL
    if (file.exists(fname)) {
        if (file.info(fname)$mtime < as.POSIXct(expiration)) {
            print(paste0("File \"", fname, "\" is out of date, recomputing"))
        } else {
            obj = readRDS(fname)
        }
    }
    if (is.null(obj)) {
        obj = func(...)
        saveRDS(obj, fname)
    }
    filelock::unlock(lk)
    return(obj)
}

#' Generate a loader for a specified expiration date
#'
#' @param expiration Time for rds expiration
#' @return A load_or_compute closure with after filled in
#' @export
get_cache_loader <- function(expiration) {
    return(function(fname, func, ...) {
        load_cached(fname, func, expiration, ...)
    })
}

#' Reverse complement
#'
#' @export
reverseComplement <- function(x) {
    return(rev(chartr(old="ACGT", new="TGCA", x)))
}
