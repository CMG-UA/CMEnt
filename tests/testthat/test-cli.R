options("CMEnt.verbose" = 0)

test_that("CLI command resolver supports subcommands and legacy flag-only invocations", {
    find_invocation <- CMEnt:::.resolveCMEntCLIInvocation(
        c("findDMRsFromSeeds", "--beta", "beta.tsv")
    )
    expect_equal(find_invocation$command, "findDMRsFromSeeds")
    expect_equal(find_invocation$command_args, c("--beta", "beta.tsv"))

    viewer_invocation <- CMEnt:::.resolveCMEntCLIInvocation(
        c("launchCMEntViewer", "--output_prefix", "results/test")
    )
    expect_equal(viewer_invocation$command, "launchCMEntViewer")
    expect_equal(viewer_invocation$command_args, c("--output_prefix", "results/test"))

    legacy_invocation <- CMEnt:::.resolveCMEntCLIInvocation(
        c("--beta", "beta.tsv", "--seeds_file", "seeds.tsv")
    )
    expect_equal(legacy_invocation$command, "findDMRsFromSeeds")
    expect_equal(
        legacy_invocation$command_args,
        c("--beta", "beta.tsv", "--seeds_file", "seeds.tsv")
    )

    help_invocation <- CMEnt:::.resolveCMEntCLIInvocation(character())
    expect_true(help_invocation$top_level_help)
})

test_that("top-level CLI help lists supported subcommands", {
    help_text <- CMEnt:::.topLevelCMEntCLIHelp("inst/bin/run_cment.R")

    expect_match(help_text, "findDMRsFromSeeds", fixed = TRUE)
    expect_match(help_text, "launchCMEntViewer", fixed = TRUE)
    expect_match(help_text, "defaults to findDMRsFromSeeds", fixed = TRUE)
})

test_that("dispatcher prints command-specific help for the viewer command", {
    skip_if_not_installed("optparse")

    help_output <- capture.output(
        CMEnt:::.runCMEntCLI(
            c("help", "launchCMEntViewer"),
            script_name = "run_cment.R"
        )
    )

    expect_true(any(grepl("launchCMEntViewer \\[options\\]", help_output)))
    expect_true(any(grepl("--output_prefix", help_output, fixed = TRUE)))
})

test_that("dispatcher routes launchCMEntViewer arguments through the CLI parser", {
    skip_if_not_installed("optparse")

    captured_args <- NULL
    local_mocked_bindings(
        launchCMEntViewerCLI = function(args) {
            captured_args <<- args
            NULL
        },
        .package = "CMEnt"
    )

    expect_no_error(
        CMEnt:::.runCMEntCLI(
            c(
                "launchCMEntViewer",
                "--output_prefix", "results/test",
                "--launch_browser", "FALSE",
                "--port", "3456",
                "--host", "0.0.0.0",
                "--diagnostic", "TRUE"
            ),
            script_name = "run_cment.R"
        )
    )

    expect_equal(captured_args$output_prefix, "results/test")
    expect_false(captured_args$launch_browser)
    expect_equal(captured_args$port, 3456L)
    expect_equal(captured_args$host, "0.0.0.0")
    expect_true(captured_args$diagnostic)
})
