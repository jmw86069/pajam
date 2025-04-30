# TODO pajam

## 28feb2024

* Download various newly available files: https://www.proteinatlas.org/about/download/

   * Subcellular location data. Specific eye toward nuclear location, areas within
     the nucleus, nucleolus, etc.
   * normal_tissue.tsv - categorical annotations: Not detected, Low, Medium, High
     for each tissue, and cell types within tissue.
     Not sure how best to use this data, but it could augment the expression matrix data
     for example put "*" in heatmap expression cells where expression is detected
     "Medium" or "High", to show which expression levels are considered "detected".
     Currently heatmaps show relative expression, apparently many low values may be
     annotated "Not detected".
