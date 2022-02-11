#' Build the shiny server ui
#'
#' @param wkdir Path to Seurat RDS files to view
#'
#' @return shiny ui
#'
#' @keywords internal
build_ui = function(wkdir){
    ui <- shinyUI(pageWithSidebar(
        ## ---- Title ----
        headerPanel("RNASeq QC"),
        ## ---- Sidebar panel for user inputs ----
        sidebarPanel(
            selectInput(inputId  = 'file.name',
                label    = 'Seurat .RDS file',
                choices  = c("Choose a file" = "", list.files(path = wkdir))),
            selectizeInput(inputId  = 'df.fields',
                label    = 'Data.frame fields',
                choices  = NULL,
                multiple = TRUE),
            downloadButton("downloadData", "Save selected")),
        ## ---- Output ----			
        mainPanel(
            h3(textOutput('caption')),
            tabsetPanel(
                ## ---- Scatterplot, UMAP ----
                tabPanel("Visualize", 
                        splitLayout(cellWidths = c("50%", "50%"),
                            plotOutput("scatterplot", 
                                height = 280*2, width = 250*2, 
                                click = "scatter_click", 
                                brush = brushOpts(id = "scatter_brush", resetOnNew=TRUE)),
                            plotOutput("umap", 
                                height = 280*2, width = 250*2, 
                                click = "umap_click", 
                                dblclick = "umap_dblclick", 
                                brush = brushOpts(id = "umap_brush", resetOnNew=TRUE)))),
                ## ---- Table ----
                tabPanel("Summary", 
                        tableOutput("filetable")),
                ## ---- Selected elements ----
                fluidRow(column(width = 10, h4("Selected:"), verbatimTextOutput("brush_info")))
            )
        )
    ))
    return(ui)
}