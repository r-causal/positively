# The layer of a figure that draws a given aesthetic, found by what it carries
# rather than by the position it was added in.
built_layer <- function(plot, aesthetic) {
  built <- ggplot2::ggplot_build(plot)
  Filter(function(layer) aesthetic %in% names(layer), built$data)[[1]]
}

# A histogram is the layer carrying a count.
histogram_bars <- function(plot) {
  built_layer(plot, "count")
}

# The fill and the colour a stock geom_histogram draws with, read off a
# reference figure. ggplot2 4.0 resolves a geom's defaults from the theme, so
# `default_aes` holds the rule rather than the value it settles on, and
# `get_geom_defaults()` postdates the declared ggplot2 floor of 3.5.0. Building
# a reference resolves the defaults through the pipeline the views under test go
# through, so what it returns is the stock appearance whatever version draws it.
stock_histogram_aes <- function() {
  reference <- ggplot2::ggplot(mapping = ggplot2::aes(x = c(0, 1, 2, 3))) +
    ggplot2::geom_histogram(bins = 2)
  bars <- histogram_bars(reference)
  list(fill = unique(bars$fill), colour = unique(bars$colour))
}
