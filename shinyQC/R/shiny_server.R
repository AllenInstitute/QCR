#' Build the shiny server
#'
#' @param wkdir to Seurat RDS files to view
#'
#' @return shiny server
#'
#' @keywords internal
build_server = function(wkdir){
    ##
    server <- shinyServer(function(input, output, session) {

        ## Define theme of page
        .theme<- theme(
            axis.line = element_line(colour = 'gray', size = .75), 
            panel.background = element_blank(),  
            plot.background = element_blank()
        )

        ## Create a reactive variable that persists through observe() style functions
        data = reactiveValues()

        ## Observe the selection drop down, if a change to reactive variable input$file.name happens then run this block
        observeEvent(input$file.name, {
            ## Load the data only once
            if(grepl("RDS", input$file.name)){ 
                print(paste0("Loading: ", input$file.name)); 
                data$rna.data  = readRDS(file.path(wkdir, input$file.name)) 
                data$rna.table = data_summary(data$rna.data@meta.data)

                ## Dynamically update the df.fields box based on the newly created data$rna.table
                updateSelectizeInput(session = session,
                                    inputId  = "df.fields",
                                    choices  = sort(unique(colnames(data$rna.table))),
                                    selected = c("seurat_clusters"), 
                                    server = TRUE)
            }else{ data$rna.data = NULL }

            if(!is.null(data$rna.data)){
                ## Draw scatterplot
                output$scatterplot <- renderPlot({
                    scatter.plot(data$rna.table)
                })

                ## Draw umap
                output$umap <- renderPlot({
                    umap(data$rna.data)
                })

                ## Draw the data table
                output$filetable <- renderTable({
                    as.data.frame(table(data$rna.data@meta.data[,c("predicted.id")]))
                })
            }
        })

        ## React if df.fields is changed, this will cause side-effects to the Selected table
        df_filtered = reactive({
            data$rna.table[,union(c("seurat_clusters", "nFeature_RNA", "nCount_RNA"), input$df.fields)]
        })

        ## Check for changes to df.fields and then update Selected table
        observeEvent(input$df.fields, {
            selected.data.frame <<- brushedPoints(df_filtered(), input$scatter_brush)

            output$brush_info <- renderPrint({
                selected.data.frame
            })
        })

        ## Watch for changes to scatter.plot brush, and redraw umap with selections
        observeEvent(input$scatter_brush, {
            selected.data.frame <<- brushedPoints(df_filtered(), input$scatter_brush)

            output$brush_info <- renderPrint({
                selected.data.frame
            })

            output$umap <- renderPlot({
                umap(data$rna.data, selected.data.frame$seurat_clusters)
            })
        })

        # When a double-click happens, check if there's a brush on the plot.
        # If so, zoom to the brush bounds; if not, reset the zoom.
        observeEvent(input$umap_dblclick, {
            brush <- input$umap_brush
            range.x <- NULL; range.y <- NULL;
            if(!is.null(brush)) {
                range.x <- c(brush$xmin, brush$xmax)
                range.y <- c(brush$ymin, brush$ymax)
            } 
            output$umap <- renderPlot({
                umap(data$rna.data, selected.data.frame$seurat_clusters, range.x=range.x, range.y=range.y)
            })
        })

        ## Downloadable csv of selected dataset ----
        output$downloadData <- downloadHandler(
            filename = function() {
                file.path(wkdir, paste(gsub(".RDS", "", input$file.name), "filter.csv", sep = "_"))
            },
            content = function(file) {
                data = data$rna.data@meta.data %>% filter(seurat_clusters %in% brushedPoints(df_filtered(), input$scatter_brush)$seurat_clusters)
                write.csv(data, file, row.names = FALSE)
            }
        )
    })
    return(server)
}