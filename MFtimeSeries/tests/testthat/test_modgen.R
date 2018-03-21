context("Test model generation")


m <- modgen(n=2, p=1)

test_that("modgen produces a stable 2-dim system of order 1", {
    expect_equal(all(abs(eigen(m$a)$values)<1), TRUE)
})

