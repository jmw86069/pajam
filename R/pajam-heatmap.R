

#' proteinatlas.org expression heatmap
#' 
#' proteinatlas.org expression heatmap
#' 
#' This function takes proteinatlas expression data `expr`,
#' and creates a heatmap of expression using a subset of
#' genes provided as `genes`.
#' 
#' By default, columns in `expr` are split by type, where
#' `colnames(expr)` are expected to have suffix `" - Type"`
#' at the end of each column name. If the columns cannot be
#' split accordingly, then all columns are assigned one
#' split name `"Expression"`.
#' 
#' To customize data for individual samples, the `expr` data
#' should be filtered before calling this function.
#' 
#' By default, the row-centering method centers each row by
#' the row minimum expression using only `"Tissue"` samples,
#' so that expression will be displayed relative to the
#' lowest tissue expression. An optional set of control samples
#' can be provided with argument `controlSamples`.
#' 
#' @family pajam heatmap
#' 
#' @return `Heatmap` produced by `ComplexHeatmap::Heatmap()` by default,
#'    when `return_type="heatmap"`; when `return_type="list"` it
#'    returns a `list` with components used in the heatmap, perhaps
#'    most important is the actual expression data matrix after
#'    expression centering, expression filtering, and sample
#'    subsetting operations, as relevant. The `list` also includes
#'    the Heatmap under element `"hm"`, so it can be plotted
#'    using `ComplexHeatmap::draw()`.
#' 
#' @param expr `numeric` matrix containing gene rows, and
#'    biological sample columns.
#' @param genes `character` vector of genes to display. These values
#'    are matched case-insensitive to the start of values in
#'    `rownames(expr)`, in order to facilitate pattern-matching.
#'    Specifically, the input is passed to
#'    `jamba::provigrep(c(genes), rownames(expr))` which iteratively
#'    matches each input `genes` using case-insensitive `grep()`.
#' @param type `character` vector of one or more expression types
#'    to include, where `"all"` includes all columns in `expr`.
#' @param row_cex,column_cex `numeric` value used to adjust row and
#'    column text size, where `1` is the default size. Text size is
#'    auto-scaled based upon the number of rows and columns being
#'    displayed; these values adjust the auto-scaled font size.
#' @param ramp `character` name of a color ramp, or a `character`
#'    vector of R colors, used to create a color ramp.
#' @param lens `numeric` value used to adjust the intensity of
#'    the color ramp, where values > 0 make colors change earlier,
#'    and values < 0 make colors change later.
#' @param color_ceiling `numeric` value indicating the maximum used
#'    for the color gradient, using units as displayed on the heatmap
#'    (exponentiated.) If `NULL` then the maximum value is used.
#' @param cluster_columns,cluster_rows `logical` indicating whether to
#'    cluster columns and rows, respectively.
#' @param centered `logical` indicating whether data should be centered
#'    by the median expression across all samples. Note that median
#'    expression in many cases is zero, for genes not widely expressed
#'    across all samples being displayed.
#' @param column_filter,row_filter numeric value which hides columns or
#'    rows when the column or row maximum expression is not at or above
#'    this numeric threshold. The value filtered is the expression value
#'    indicated on the heatmap, the normal expression value, not
#'    log2-transformed.
#' @param controlSamples `character` vector, or `NULL`, indicating
#'    specific samples to use when centering data by row. When
#'    `controlSamples=NULL`, the default is to use all Tissue samples.
#' @param gene_names `logical` indicating whether to display full gene
#'    names, provided using `genejam::freshenGenes()`.
#' @param gene_im `matrix` intended to be an incidence matrix, whose
#'    values are `0` and `1`, with gene symbols as rownames, and various
#'    columns indicating different annotations. When a rowname in the
#'    heatmap is not found in the incidence matrix,
#'    it subtitutes `0` to fill the empty space. See examples to see how
#'    to use the `proteinatlas_genesets_fdb11` data.
#' @param gene_im_colors `function` or `vector` that contains data
#'    to define colors for `gene_im`, using methods sufficient for
#'    `ComplexHeatmap::HeatmapAnnotation()`. When `gene_im_colors` is
#'    `NULL` a default function is used which is intended only for
#'    `c(0, 1)`.
#' @param left_annotation optional heatmap annotation as produced
#'    by `ComplexHeatmap::HeatmapAnnotation()`, as an alternative to
#'    supplying `gene_im`. Note than when `gene_im` is supplied,
#'    it will replace `left_annotation`. This option is available
#'    only to allow supplying a customized side annotation.
#' @param fill_missing `logical` indicating whether the input `genes` should
#'    bo used as-is with no pattern matching, and by substituting `0`
#'    for any missing entries. This argument is useful when trying to
#'    align this heatmap with another existing `Heatmap` object, where
#'    the rows must be exactly aligned.
#' @param border `logical` indicating whether to draw a heatmap border.
#' @param useCenterGroups `logical` used when `centered=TRUE`, and when
#'    `controlSamples` are provided, will center each sample type
#'    independently, instead of centering all samples versus `"Tissue"`
#'    type. We decided by default to center all samples versus `"Tissue"`
#'    by default, because some immune-specific genes are highly expressed
#'    in all `"Blood"` samples but not `"Tissue"`, and it was visually
#'    confusing.
#' @param trim_columns `logical` indicating whether to trim the colnames
#'    to remove the sample type suffix.
#' @param rowStatsFunc `function` used to define the per-row value used
#'    when `centered=TRUE`, passed to `jamma::centerGeneData()`.
#'    The default is to use the `matrixStats::rowMins()` which returns
#'    the minimum expression value per row.
#' @param return_type `character` string indicating the type of return
#'    object: `"heatmap"` returns the `Heatmap` object sufficient to
#'    be plotted using `ComplexHeatmap::draw()`; `"list"` returns `list`
#'    of relevant data components used to produce the heatmap, sufficient
#'    for reviewing more details. The `list` output also includes data
#'    actually used in the heatmap, after the expression centering,
#'    expression filtering, and sample subsetting, as relevant.
#'    if relevant.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `ComplexHeatmap::Heatmap()`
#'    for other customizations.
#' 
#' @examples 
#' test_genes <- c("DKK1","DKK4","CXCL12","IL6R","MET",
#'    "HK2","FTL","FTH1","STAT1","STAT3","CDKN1B");
#' proteinatlas_heatmap(genes=test_genes,
#'    type="Blood",
#'    centered=FALSE,
#'    cluster_rows=TRUE,
#'    cluster_columns=TRUE,
#'    row_filter=2)
#' 
#' # use proteinatlas_genesets_fdb11
#' use_im <- c("secreted_proteins",
#'    "membrane_proteins",
#'    "NOT_membrane_secreted",
#'    "TFs");
#' proteinatlas_im <- list2im_opt(proteinatlas_genesets_fdb11[use_im]);
#' test_genes <- c("DKK1","DKK4","CXCL12","IL6R","MET",
#'    "HK2","FTL","FTH1","STAT1","STAT3","CDKN1B");
#' proteinatlas_heatmap(genes=test_genes,
#'    type="Blood",
#'    centered=TRUE,
#'    gene_im=proteinatlas_im);
#' 
#' @import jamba
#' 
#' @export
proteinatlas_heatmap <- function
(expr=proteinatlas_expr_fdb11,
   genes=NULL,
   samples=NULL,
   type=c("all", "Tissue", "Cell", "Blood", "Brain"),
   row_cex=1,
   column_cex=1,
   ramp="Reds",
   lens=0,
   color_ceiling=NULL,
   cluster_columns=FALSE,
   cluster_rows=FALSE,
   column_split=NULL,
   row_split=NULL,
   centered=FALSE,
   row_filter=0,
   column_filter=0,
   controlSamples=NULL,
   gene_names=FALSE,
   gene_im=NULL,
   gene_im_colors=NULL,
   left_annotation=NULL,
   fill_missing=TRUE,
   border=TRUE,
   useCenterGroups=TRUE,
   trim_columns=FALSE,
   rowStatsFunc=matrixStats::rowMins,
   return_type=c("heatmap", "list"),
   verbose=FALSE,
   ...)
{
   type <- match.arg(type,
      several.ok=TRUE);
   return_type <- match.arg(return_type);
   
   if (length(genes) == 0) {
      use_genes <- jamba::mixedSort(rownames(expr));
   } else {
      if (all(genes %in% rownames(expr))) {
         use_genes <- genes;
      } else {
         if (fill_missing) {
            missing_genes <- setdiff(genes, rownames(expr));
            expr_new <- jamba::rbindList(
               lapply(jamba::nameVector(missing_genes), function(g){
                  g1 <- expr[1,,drop=FALSE] * 0;
                  rownames(g1) <- g;
                  g1;
               }));
            expr <- rbind(expr, expr_new);
            use_genes <- genes;
         } else {
            gene_grep <- gsub("^([A-Za-z])", "^\\1", genes);
            use_genes <- jamba::provigrep(gene_grep,
               sortFunc=jamba::mixedSort,
               jamba::mixedSort(rownames(expr)));
         }
      }
   }
   if (length(use_genes) == 0) {
      stop("No usable genes found.");
   }
   expr <- expr[use_genes,,drop=FALSE];
   
   if (gene_names) {
      gene_labels <- genejam::freshenGenes(use_genes,
         final=c("GENENAME"),
         empty_rule="original")$GENENAME;
      gene_labels <- ifelse(gene_labels == use_genes,
         use_genes,
         paste0(use_genes, " - ", gene_labels));
      #jamba::printDebug("mean gene nchar:", mean(nchar(use_genes)));
      #jamba::printDebug("mean label nchar:", mean(nchar(gene_labels)));
      nchar_ratio <- (mean(nchar(use_genes)))^(1/4) / (mean(nchar(gene_labels)))^(1/4);
      row_cex <- row_cex * nchar_ratio;
      #jamba::printDebug("nchar_ratio:", nchar_ratio);
   } else {
      gene_labels <- use_genes;
   }
   
   if (length(samples) > 0) {
      if (all(samples %in% colnames(expr))) {
         use_samples <- samples;
      } else {
         use_samples <- jamba::provigrep(samples, colnames(expr));
      }
      if (length(use_samples) == 0) {
         stop("No usable samples found.");
      }
      expr <- expr[,use_samples,drop=FALSE];
   }
   
   expr_m <- log2(1 + jamba::rmNA(naValue=0, expr));
   if (centered) {
      if (length(controlSamples) == 0) {
         controlSamples <- jamba::vigrep("- Tissue", colnames(expr));
      }
      if (useCenterGroups) {
         centerGroups <- gsub("^.+ - ", "", colnames(expr));
      } else {
         centerGroups <- NULL;
      }
      expr_m <- jamma::centerGeneData(expr_m,
         rowStatsFunc=rowStatsFunc,
         centerGroups=centerGroups,
         mean=FALSE,
         controlSamples=controlSamples);
   }
   
   if ("all" %in% type) {
   } else {
      use_colnames <- jamba::provigrep(paste0(" - .*", type),
         colnames(expr));
      expr <- expr[,use_colnames,drop=FALSE];
   }
   if (length(column_split) == 0) {
      column_split <- gsub("^.+ - ", "", colnames(expr));
      if (length(unique(column_split)) == ncol(expr)) {
         column_split <- rep("Expression", ncol(expr));
      }
      column_split <- factor(column_split,
         levels=unique(column_split));
   }
   if (length(unique(column_split)) > 1) {
      #u_colnames <- gsub(" - .+$", "", colnames(expr));
      #dupe_u_colnames <- names(tcount(u_colnames, minCount=2));
      ordered_colnames <- jamba::provigrep(c("tissue", "."), colnames(expr));
      u_colnames <- gsub("cells", "cell",
         gsub(" - .+$", "", 
            ordered_colnames));
      dupe_u_colnames <- names(jamba::tcount(u_colnames, minCount=2));
      if (length(dupe_u_colnames) > 0) {
         keep_colnames <- (!colnames(expr) %in% 
               ordered_colnames[match(dupe_u_colnames, u_colnames)]);
         if (verbose) {
            jamba::printDebug("proteinatlas_heatmap(): ",
               "Dropped duplicate colnames: ",
               sep=", ",
               setdiff(colnames(expr),
                  colnames(expr)[keep_colnames]));
         }
         
         expr <- expr[,keep_colnames,drop=FALSE];
         column_split <- factor(column_split[keep_colnames]);
      }
   }
   
   expr_m <- expr_m[rownames(expr),colnames(expr),drop=FALSE];
   
   if (length(row_filter) > 0 && row_filter != 0) {
      row_keep <- (matrixStats::rowMaxs(abs(expr_m)) >= log2(row_filter));
      gene_labels <- gene_labels[row_keep];
      use_genes <- use_genes[row_keep];
      expr_m <- expr_m[row_keep,,drop=FALSE];
      # version 0.0.4.900 also subset row_split if supplied
      if (length(row_split) > 0) {
         row_split <- factor(row_split[row_keep]);
      }
   }
   
   if (length(column_filter) > 0 && column_filter != 0) {
      column_keep <- (matrixStats::colMaxs(abs(expr_m)) >= log2(column_filter));
      if (length(column_split) == ncol(expr_m)) {
         column_split <- factor(column_split[column_keep]);
      }
      expr_m <- expr_m[,column_keep,drop=FALSE];
   }
   
   ## Automatic fontsize
   row_fontsize <- jamba::noiseFloor(
      row_cex * 60/(nrow(expr_m))^(1/2),
      ceiling=20,
      minimum=2);
   column_fontsize <- jamba::noiseFloor(
      column_cex * 60/(ncol(expr_m))^(1/2),
      ceiling=20,
      minimum=2);
   
   ## optional row annotations
   if (length(gene_im) > 0) {
      if (!all(use_genes %in% rownames(gene_im))) {
         missing_genes <- setdiff(use_genes, rownames(gene_im));
         if (verbose) {
            jamba::printDebug("proteinatlas_heatmap(): ",
               "Adding missing_genes to gene_im:",
               missing_genes);
         }
         gene_im_new <- jamba::rbindList(
            lapply(jamba::nameVector(missing_genes), function(g){
               g1 <- gene_im[1,,drop=FALSE] * 0;
               rownames(g1) <- g;
               g1;
            }));
         gene_im <- rbind(gene_im, gene_im_new);
      }
      if (length(gene_im_colors) == 0) {
         gene_im_at <- jamba::mixedSort(keepNegative=TRUE,
            unique(as.vector(unlist(gene_im))));
         if (length(gene_im_at) > 5) {
            gene_im_at <- NULL;#pretty(range(gene_im_at, na.rm=TRUE));
         }
         gene_im_colors <- circlize::colorRamp2(
            breaks=c(-1, 0, 1),
            colors=c("dodgerblue3", "white", "red3"))
      } else {
         gene_im_at <- NULL;
      }
      left_annotation <- ComplexHeatmap::rowAnnotation(
         name="gene_im",
         gene_im=gene_im[use_genes,,drop=FALSE],
         border=TRUE,
         annotation_name_gp=grid::gpar(fontsize=column_fontsize),
         col=list(gene_im=gene_im_colors),
         show_legend=TRUE,
         annotation_legend_param=list(
            border=TRUE,
            color_bar="discrete",
            at=gene_im_at)
      );
   }
   
   ## determine useful color range
   if (centered) {
      if (length(color_ceiling) == 0) {
         fc_max <- max(abs(expr_m), na.rm=TRUE);
      } else {
         fc_max <- log2(abs(color_ceiling));
      }
      if ("Reds" %in% ramp) {
         ramp <- "RdBu_r";
      }
      colbr <- circlize::colorRamp2(
         breaks=seq(from=-fc_max, to=fc_max, length.out=11),
         colors=jamba::getColorRamp(ramp,
            divergent=TRUE,
            n=11,
            lens=lens));
      label_at <- seq(from=-ceiling(fc_max), to=ceiling(fc_max));
      label_at <- seq(from=ceiling(min(expr_m, na.rm=TRUE)),
         to=floor(min(c(fc_max, max(expr_m, na.rm=TRUE)))));
      label_text <- 2^(abs(label_at)) * sign(label_at+1e-10);
      name <- "proteinatlas.org\ncentered\nexpression";
   } else {
      if (length(color_ceiling) == 0) {
         expr_max <- 2^max(expr_m, na.rm=TRUE)-1;
      } else {
         expr_max <- color_ceiling;
      }
      colbr <- circlize::colorRamp2(
         breaks=seq(from=0, to=ceiling(log2(1+expr_max)), length.out=11),
         colors=c("white",
            jamba::getColorRamp(ramp,
               n=10,
               lens=lens,
               trimRamp=c(0, 1))));
      label_x <- seq(from=0, to=ceiling(log2(1+expr_max)));
      label_at <- c(0, log2(1+2^label_x))
      label_text <- c(0, 2^label_x);
      name <- "proteinatlas.org\nexpression";
   }
   
   if (trim_columns) {
      column_labels <- gsub(" - .+", "", colnames(expr_m));
   } else {
      column_labels <- colnames(expr_m);
   }
   
   hm <- ComplexHeatmap::Heatmap(
      border=border,
      expr_m,
      name=name,
      row_names_gp=grid::gpar(fontsize=row_fontsize),
      cluster_rows=cluster_rows,
      cluster_row_slices=FALSE,
      column_labels=column_labels,
      column_names_gp=grid::gpar(fontsize=column_fontsize),
      cluster_columns=cluster_columns,
      cluster_column_slices=FALSE,
      column_split=column_split,
      row_split=row_split,
      row_labels=gene_labels,
      left_annotation=left_annotation,
      heatmap_legend_param=list(
         border=TRUE,
         at=label_at,
         labels=label_text,
         color_bar="discrete"),
      col=colbr,
      ...);
   if ("heatmap" %in% return_type) {
      return(hm);
   }
   ret_list <- list(
      expr_m=expr_m,
      hm=hm,
      col=colbr,
      row_split=row_split,
      column_split=column_split,
      column_labels=column_labels
   );
   if (length(gene_im) > 0) {
      ret_list$gene_im <- gene_im[use_genes,,drop=FALSE];
      ret_list$left_annotation <- left_annotation;
   }
   
}

#' List to incidence matrix
#' 
#' List to incidence matrix
#' 
#' This function takes a list of `character` vectors and returns
#' a `numeric` `matrix` whose rownames are entries from the
#' `character` vectors, colnames are `names` from the `setlist`,
#' and whose values are `0` or `1`.
#' 
#' @family pajam list
#' 
#' @examples
#' setlist <- list(setA=LETTERS[1:10], setB=LETTERS[7:14]);
#' list2im_opt(setlist);
#' 
#' @export
list2im_opt <- function
(setlist,
 ...)
{
   setnamesunion <- Reduce("union", setlist);
   setlistim <- do.call(cbind, lapply(setlist, function(i){
      i_match <- match(i, setnamesunion);
      j <- rep(0, length(setnamesunion));
      j[i_match] <- 1;
      j;
   }))
   rownames(setlistim) <- setnamesunion;
   return(setlistim);
}

