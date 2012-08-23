KnitPost <- function(input, base.url = "/") {
    require(knitr)

    opts_knit$set(base.url = base.url)
    fig.path <- paste0("images/", sub(".Rmd$", "", basename(input)), "/")
    output   <- sub(".Rmd", ".markdown", input)

    opts_chunk$set(fig.path = fig.path)
    opts_chunk$set(fig.cap = "center")

    render_jekyll()
    knit(input, output = output, envir = parent.frame())
}