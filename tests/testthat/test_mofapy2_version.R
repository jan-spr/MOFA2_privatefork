library(MOFA2)

# An explicit pin is passed throughout so these keep passing when .mofapy2_version is bumped.

test_that("matching and patch-level versions let training proceed", {
    expect_silent(.check_mofapy2_version("0.7.3", "0.7.3"))
    expect_message(.check_mofapy2_version("0.7.4", "0.7.3"), "patch release ahead")
    expect_warning(.check_mofapy2_version("0.7.2", "0.7.3"), "patch release behind")
})

test_that("minor and major mismatches stop training, in either direction", {
    expect_error(.check_mofapy2_version("0.8.0", "0.7.3"), "minor version ahead of")
    expect_error(.check_mofapy2_version("0.6.9", "0.7.3"), "minor version behind")
    expect_error(.check_mofapy2_version("1.0.0", "0.7.3"), "major version ahead of")
    expect_error(.check_mofapy2_version("0.7.3", "1.0.0"), "major version behind")
})

test_that("awkward version strings are handled", {
    # ordered numerically, not as strings: 10 is ahead of 9
    expect_message(.check_mofapy2_version("0.7.10", "0.7.9"), "patch release ahead")
    # package_version("0.7")$patch is NA
    expect_warning(.check_mofapy2_version("0.7", "0.7.3"), "patch release behind")
    # PEP 440 builds are not valid R version strings
    expect_error(.check_mofapy2_version("0.7.4.dev0", "0.7.3"), "Only numeric versions")
})
