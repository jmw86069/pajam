# pajam version 0.0.3.900

## changes to existing functions

* `proteinatlas_heatmap()` new argument `gene_im_colors`
for custom color ranges. Also `gene_im` data is expected
to have values `c(-1, 0, 1)` instead of just `c(0, 1)`
previously. Also the `gene_im` color legend is displayed
by default.

# pajam version 0.0.2.900

## changes to existing functions

* `pajam_shiny_ui()` fixed a simple typo in `selected_genes` that
defines the initial starting set of genes displayed in the shiny
app.
* `'launch_pajam()` via `pajam_shiny_ui()` will honor the
value of `selected_genes` if defined in the global environment
`.GlobalEnv`, otherwise it will use the default set of genes.

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
