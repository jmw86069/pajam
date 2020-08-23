# pajam version 0.0.1.900

## new functions

* `launch_pajam()` launches an R-shiny app, to help navigate the
data. The app uses really nice interactive features from the
ComplexHeatmap package.

# pajam version 0.0.0.900

## initial release

* `proteinatlas_heatmap()` the main function which creates a heatmap
using `ComplexHeatmap::Heatmap()`.
* `proteinatlas_expr_fdb11` is the central expression data matrix used
for the heatmap.
* `proteinatlas_genesets_fdb11` is a `list` of useful genesets that
can be used as annotations beside heatmap rows.
