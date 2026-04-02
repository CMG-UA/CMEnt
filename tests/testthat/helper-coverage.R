is_covr_run <- function() {
    identical(tolower(Sys.getenv("R_COVR")), "true")
}

skip_if_covr_expensive <- function(reason = "Skipping expensive integration tests under covr.") {
    if (is_covr_run()) {
        testthat::skip(reason)
    }
}
