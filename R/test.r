lmtest <- function()
{
    set.seed(1234)
    x <- matrix(rnorm(30), 10, 3)
    y <- rnorm(10)

    lm_fit(x, y)
}
