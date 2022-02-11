#' Launch the shinyQC server
#'
#' This function initializes and launches a shiny server focused on `wkdir`. To connect locally
#' to a server running on a cluster environemnt use tunnelling: 
#'
#'      ssh -t AIBS.id@cluster.name -L 5000:localhost:5000
#'
#' @param wkdir Path to Seurat RDS files to view
#' @param port Server port that will be exposed for connections
#'
#' @import shiny
#' @import shinyjs
#' @import data.table
#' @import tidyr
#' @import dplyr
#' @import Seurat
#' @import ggplot2
#' @import gridExtra
#'
#' @export
shinyQC = function(wkdir=".", port=5000){
    
    ##
    options(shiny.maxRequestSize=200*1024^2)

    ##
    ui = build_ui(wkdir)
    server = build_server(wkdir)

    ##
    runApp(shinyApp(ui = ui, server = server), launch.browser = FALSE, port=port)
}
