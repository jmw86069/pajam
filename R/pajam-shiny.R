# pajam-shiny.R
# functions to enable an R-shiny app

#' R-shiny UI for Protein Atlas visualization
#' 
#' @import shiny
#' @import shinydashboard
#' @import shinydashboardPlus
#' @import shinyWidgets
#'
#' @param ... additional arguments are ignored.
#' 
#' @export
pajam_shiny_ui <- function
(...)
{
   ##
   # header
   header <- dashboardHeaderPlus(
      title=tagList("pajam v0.1 Protein Atlas Shiny")
   );
   
   all_genes <- rownames(proteinatlas_expr_fdb11);
   selected_genes <- c("DKK1","DKK4","CXCL12","IL6R","MET",
      "HK2","FTL","FTH1","STAT1","STAT3","CDKN1B");
   all_sampletypes <- c("all", "Tissue", "Cell", "Blood", "Brain");
   selected_sampletypes <- c("Tissue", "Cell");
   all_annotations <- names(proteinatlas_genesets_fdb11);
   selected_annotations <- c("secreted_proteins",
      "membrane_proteins",
      "NOT_membrane_secreted",
      "TFs");

   sidebar <- dashboardSidebar(
      sidebarMenu(
         id="tabs",
         menuItem(
            text="Protein Atlas Heatmap",
            tabName="heatmapplot",
            icon=icon("bezier-curve")),
         selectizeInput(inputId="filter_genes",
            label="Genes to display",
            choices=all_genes,
            selected=selected_genes,
            multiple=TRUE,
            options = list(
               splitOn=I("(function() { return /[,;\\n\\r ]+/; })()"),
               create=I("function(input, callback){ return { value: input, text: input }; }")
            )
         ),
         selectizeInput(inputId="filter_sampletypes",
            label="Sample types",
            choices=all_sampletypes,
            selected=selected_sampletypes,
            multiple=TRUE),
         selectizeInput(inputId="use_im",
            label="Include annotations:",
            choices=all_annotations,
            selected=selected_annotations,
            multiple=TRUE),
         selectInput(inputId="cluster_rows",
            label="Cluster rows:",
            choices=c(TRUE, FALSE),
            selected=TRUE),
         selectInput(inputId="cluster_columns",
            label="Cluster columns:",
            choices=c(TRUE, FALSE),
            selected=TRUE),
         actionButton(
            inputId="apply_filter",
            label="Update"),
         #menuItem(
         #   text="Interactome Filtering",
         #   tabName="filtering",
         #   icon=icon("th-list", lib="glyphicon")),
         #menuItem(
         #   text="Exploratory Plots",
         #   tabName="exploratoryplots",
         #   icon=icon("chart-area")),
         menuItem(
            text="Help",
            tabName="guides",
            icon=icon("info"))
         #menuItem(
         #   text="Samples and Data",
         #   tabName="samplesdata",
         #   icon=icon("table"))
      )
   );

   
   heatmapTab <- fluidPage(
      fluidRow(
         column(
            width=6,
            style="padding:0px",
            shinydashboard::box(
               width="100%",
               height="100%",
               fluidRow(
                  column(
                     width=12,
                     shinyjqui::jqui_resizable(
                        plotOutput("main_heatmap",
                           height="500px",
                           width="100%",
                           brush="ht_brush",
                           click="ht_click")
                     )
                  )
               )
            )
         ),
         column(width=6,
            style="padding:0px",
            shinydashboard::box(
               width="100%",
               height="100%",
               fluidRow(
                  column(
                     width=12,
                     shinyjqui::jqui_resizable(
                        plotOutput("sub_heatmap",
                           height="500px",
                           width="100%")
                     )
                  )
               )
            )
         )
      ),
      verbatimTextOutput("ht_click_content")
   );
   
   # define guides tab
   guidesTab <- fluidPage(
      tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
      sidebarLayout(
         mainPanel(
            width=7,
            tabBox(
               width=12,
               tabPanel(
                  title="Protein Atlas Interactive Heatmap",
                  uiOutput("pajam_guide")
               ),
               tabPanel(
                  title="Protein Atlas Data",
                  uiOutput("pajam_data")
               )
            )
         ),
         sidebarPanel(
            width=5,
            "ProteinAtlas.org data is provided in an interactive
            heatmap to allow querying and visualization of specific
            genes, in a variety of human tissue sources and cell lines.",
            tags$ul(
               tags$li(
                  strong(style="color:firebrick",
                     "The BioRxiv record is available:"),
                  br(),
                  a("D.E.Gordon, et al",
                     "BioRxiv (Preprint). ",
                     em(" A SARS-CoV-2-Human Protein-Protein Interaction Map Reveals Drug Targets and Potential Drug-Repurposing."),
                     href="https://www.biorxiv.org/content/10.1101/2020.03.22.002386v1.full.pdf")
               )
            ),
            tags$p("Relevant R version info:"),
            tags$ul(
               tags$li(
                  strong(style="color:black", R.version.string)
               ),
               tags$li(
                  strong(style="color:black", "jamba:"),
                  as.character(packageVersion("jamba"))
               ),
               tags$li(
                  strong(style="color:black", "pajam:"),
                  as.character(packageVersion("pajam"))
               ),
               tags$li(
                  strong(style="color:black", "colorjam:"),
                  as.character(packageVersion("colorjam"))
               ),
               tags$li(
                  strong(style="color:black", "shiny:"),
                  as.character(packageVersion("shiny"))
               ),
               tags$li(
                  strong(style="color:black", "ComplexHeatmap:"),
                  as.character(packageVersion("ComplexHeatmap"))
               )
            )
         )
      )
   );

   # dashboard body
   body <- dashboardBody(
      #shinyjs::useShinyjs(),
      #setShadow(class="box"),
      #setShadow(class="boxPlus"),
      tabItems(
         tabItem(tabName="heatmapplot", heatmapTab),
         tabItem(tabName="guides", guidesTab)
         #tabItem(tabName="filtering", filteringTab)
      )
   );

   # assemble the page
   dp <- dashboardPage(
      header,
      sidebar,
      body,
      skin="blue");
   dp;
   
}


get_pajam_guides <- function
(...)
{
   # interactome_guide (in the central panel of the help page)
   pajam_data <- fluidPage(
      h1("About the pajam Protein Atlas content",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$h3("Available content"),
         tags$p(
            "ProteinAtlas.org provides downloadable content for human genes
            which have protein and transcript expression data across a
            range of human tissues and cell lines."),
         tags$p(
            "ProteinAtlas.org also provides several gene annotations we
            found useful to add context to the heatmaps, a selection of
            which are described below:"),
         tags$ul(
            tags$li(
               strong("Secreted", style="color:navy"),
               " - indicates the protein product is secreted from human cells."
            ),
            tags$li(
               strong("Membrane", style="color:navy"),
               " - indicates the protein product is transmembrane-associated."
            ),
            tags$li(
               strong("NOT_secteted_membrane", style="color:navy"),
               " - indicates a protein product is neither secreted, nor
               membrane-associated, in any known context. Note this annotation
               does not mean the protein cannot be cytoplasmic, it means the
               protein is not *always* cytoplasmic."
            ),
            tags$li(
               strong("FDA", style="color:navy"),
               " - indicates the protein product is the target of
               an FDA-approved therapeutic agent."
            ),
            tags$li(
               strong("TF", style="color:navy"),
               " - indicates the protein product is a transcription factor."
            )
         )
      )
   );
   
   # pajam_guide
   pajam_guide <- fluidPage(
      h1("About the Protein Atlas heatmap",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$h3("How to navigate the pajam heatmap"),
         tags$p(
            "The heatmap is visualized using the Bioconductor
            package ComplexHeatmap."),
         tags$ul(
            tags$li(
               strong("Zoom", style="color:navy"),
               " - Click an drag to create a rectangle on the left
               heatmap, to zoom into the same region for display
               on the right side."
            ),
            tags$li(
               strong("Pan", style="color:navy"),
               " - Pan by click-and-drag with a mouse, or using the arrow buttons."
            ),
            tags$li(
               strong("Hover", style="color:navy"),
               " - Hovering a node or edge will reveal tooltip text, with more
               information depending upon the data available. Gene nodes will
               include the full gene name, aliases, and optionally the official
               gene symbol. Pathways will show the full pathway name. Colorized
               edges will show the Compound, Drug Target, Drug Status, as
               provided by the authors."
            )
         )
      )
   )
   return(list(
      pajam_guide=pajam_guide,
      pajam_data=pajam_data
   ));
}


#' R-shiny UI for Protein Atlas visualization
#' 
#' @export
pajam_shiny_server <- function
(input,
 output,
 session,
 shiny_env=new.env(),
 ...)
{
   #
   options("warn"=-1);
   #output$sashimiplot_guide <- renderUI(sashimiplot_guide);
   #output$sashimiplotviz_guide <- renderUI(sashimiplotviz_guide);
   pajam_guides <- get_pajam_guides();
   
   output$pajam_guide <- renderUI(pajam_guides$pajam_guide);
   output$pajam_data <- renderUI(pajam_guides$pajam_data);

   # update the "Update" button when something has changed
   observe({
      ## only enable the button when a change is made
      filter_genes <- intersect(input$filter_genes,
         rownames(proteinatlas_expr_fdb11));
      filter_sampletypes <- input$filter_sampletypes;
      
      shinyjs::enable("apply_filter");
   });
   
   ## isolate() makes input refresh only when reactive() is triggered
   ## in this case input$apply_filter is the only trigger
   get_filters <- reactive({
      input$apply_filter;
      filter_genes <- isolate(input$filter_genes);
      filter_sampletypes <- isolate(input$filter_sampletypes);
      cluster_rows <- isolate(input$cluster_rows %in% "TRUE");
      cluster_columns <- isolate(input$cluster_columns %in% "TRUE");
      return(list(
         genes=filter_genes,
         sampletypes=filter_sampletypes,
         cluster_rows=cluster_rows,
         cluster_columns=cluster_columns));
   });

   get_im_data <- reactive({
      use_im <- input$use_im;
      if (length(use_im) > 0) {
         proteinatlas_im <- list2im_opt(proteinatlas_genesets_fdb11[use_im]);
      } else {
         proteinatlas_im <- NULL;
      }
      return(proteinatlas_im);
   });
   
   # function to return expression matrix to display
   # proteinatlas_expr_fdb11
   # proteinatlas_genesets_fdb11
   get_heatmap_data <- reactive({
      use_filters <- get_filters();
      proteinatlas_im <- get_im_data();
      ht <- proteinatlas_heatmap(genes=use_filters$genes,
         type=use_filters$sampletypes,
         centered=TRUE,
         gene_im=proteinatlas_im,
         cluster_rows=use_filters$cluster_rows,
         cluster_columns=use_filters$cluster_columns,
         row_filter=2);
      return(ht);
   });
   
   output$main_heatmap <- renderPlot({
      ht <- get_heatmap_data();
      shiny_env$ht <- draw(ht);
      shiny_env$ht_pos = ht_pos_on_device(shiny_env$ht);
   });
   
   output$sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
         grid.newpage()
         grid.text("No region is selected.", 0.5, 0.5)
      } else {
         lt <- ComplexHeatmap:::get_pos_from_brush(input$ht_brush);
         pos1 <- lt[[1]];
         pos2 <- lt[[2]];
         
         ht <- shiny_env$ht;
         pos <- selectArea(ht,
            mark=FALSE,
            pos1=pos1,
            pos2=pos2,
            verbose=FALSE,
            ht_pos=shiny_env$ht_pos);

         row_index <- unlist(pos[1, "row_index"]);
         column_index <- unlist(pos[, "column_index"]);
         column_slice_index <- rep(unlist(pos[, "column_slice"]),
            lengths(pos[, "column_index"]));
         printDebug("ht@column_title:", ht@column_title);
         printDebug("column_slice_index:");print(column_slice_index);
         column_slice <- names(column_order(ht))[column_slice_index];
         if (all(is.na(column_slice))) {
            column_slice <- NULL;
         }
         #printDebug("column_slice:");print(column_slice);
         m <- ht@ht_list[[1]]@matrix;
         ht_select <- Heatmap(
            border=TRUE,
            m[row_index, column_index, drop=FALSE],
            col=ht@ht_list[[1]]@matrix_color_mapping@col_fun,
            show_heatmap_legend=FALSE,
            column_split=column_slice,
            cluster_rows=FALSE,
            cluster_columns=FALSE);
         draw(ht_select);
      }
   });
   
   output$ht_click_content <- renderText({
      if (is.null(input$ht_click)) {
         "Not selected."
      } else {
         pos1 <- ComplexHeatmap:::get_pos_from_click(input$ht_click);
         jamba::printDebug("pos1:");
         print(sdim(pos1));
         
         ht <- shiny_env$ht;
         pos <- selectPosition(ht,
            mark=FALSE,
            pos=pos1,
            verbose=FALSE,
            ht_pos=shiny_env$ht_pos);
         
         row_index <- pos[1, "row_index"];
         column_index <- pos[1, "column_index"];
         m <- ht@ht_list[[1]]@matrix;
         v <- m[row_index, column_index];
         
         glue::glue("row index: {row_index}",
            "column index: {column_index}",
            "value: {v}", .sep = "\n");
      }
   });
   
}


#' Launch pajam shiny app
#' 
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import glue
#' @import ComplexHeatmap
#' 
#' @export
launch_pajam <- function
(...,
 width=2200,
 port=8888,
 host="0.0.0.0",
 options=list(width=width,
    host=host,
    port=port))
{
   suppressPackageStartupMessages(require(shiny));
   suppressPackageStartupMessages(require(shinydashboard));
   #suppressPackageStartupMessages(require(shinydashboardPlus));
   suppressPackageStartupMessages(require(ComplexHeatmap));
   options("shiny.host"=host);
   jamba::printDebug("launch_pajam(): ",
      "host:", host);
   jamba::printDebug("launch_pajam(): ",
      "port:", port);
   options("shiny.port"=as.numeric(port));
   
   shiny::shinyApp(ui=pajam_shiny_ui,
      server=pajam_shiny_server,
      #onStart=sashimiAppConstants,
      options=options
   );
}
