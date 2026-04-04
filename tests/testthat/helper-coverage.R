is_covr_run <- function() {
    identical(tolower(Sys.getenv("R_COVR")), "true")
}