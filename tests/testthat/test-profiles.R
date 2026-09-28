test_that("the profile counts insertions per site and base", {
    tfbs <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1001, 2001), width = 201)
    )
    ins <- GenomicRanges::GRanges("chr1",
        IRanges::IRanges(c(1101, 2101), width = 1),
        score = c(2, 4)
    )
    prof <- accProfile(ins, tfbs)
    expect_true(all(c("x", "avg_ins") %in% names(prof)))
    expect_equal(prof[x == 0, avg_ins], 3)
    expect_equal(sum(prof$avg_ins), 3)
})

test_that("an uncovered motif returns NULL", {
    tfbs <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1001, width = 201))
    ins <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1101, width = 1),
        score = 1
    )
    expect_null(accProfile(ins, tfbs))
})
