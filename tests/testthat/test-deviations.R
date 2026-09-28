test_that("a flat profile has a deviation of one", {
    prof <- data.table::data.table(x = -250:250, avg_ins = 1)
    expect_equal(accDeviationScore(prof), 1)
})

test_that("a central dip lowers the deviation", {
    prof <- data.table::data.table(x = -250:250, avg_ins = 1)
    prof[abs(x) <= 25, avg_ins := 0.2]
    expect_lt(accDeviationScore(prof), 1)
})

test_that("the deviation ignores the scale of the profile", {
    prof <- data.table::data.table(x = -250:250, avg_ins = 1)
    prof[abs(x) <= 25, avg_ins := 0.2]
    scaled <- data.table::copy(prof)[, avg_ins := avg_ins * 17]
    expect_equal(accDeviationScore(prof), accDeviationScore(scaled))
})

test_that("row z-scores are centred", {
    z <- computeRowZScore(matrix(c(1, 2, 3, 10, 20, 30), nrow = 2,
        byrow = TRUE
    ))
    expect_equal(rowSums(z), c(0, 0))
})
