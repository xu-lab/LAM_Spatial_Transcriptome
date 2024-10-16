#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# In this app...
# Load in RDS file
# Display the 2-D scatterplot (UMAP)
# - checkboxes for each cluster

# Violin plots of gene expression levels in each cell cluster
# - user chooses a gene

# Lists of top biomarkers for every cluster

# Names of clusters?

library(shiny)
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("LAM Spatial Transcriptome"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            # p("Please be mindful that RDS data may take more than 30 seconds to load and cluster."),
            # Which file is to be selected?
            selectInput("file", "Choose file:", choices = list.files("./spatialdata")),
            # Select dimensionality reduction type
            selectInput("red_type", "Select dimensionality reduction method for scatterplot:", choices = c("umap", "tsne")),
            
            tags$hr(style="margin: 30px -20px 20px; border: 0; border-top: 3px solid #C0C0C0;"),
            # Make the verbatimTextOutput scrollable: https://forum.posit.co/t/verbatimtextoutput-sizing-and-scrollable/1193
            tags$head(tags$style("#gene_list{color:black; font-size:12px;  overflow-y:scroll; max-height: 200px; background: ghostwhite;}")),
            
            # The biomarkers will appear here. The user can view the levels of one or more of them across each cluster
            uiOutput("biomarkers"),
        ),

        # Show a plot of the generated distribution
        # Want to add a dividing line:
        # https://stackoverflow.com/questions/69859938/vertical-divider-between-columns-in-shiny
        mainPanel(
          # Want to put the tissue image and the dimensionality reduction plot side by side
          # https://stackoverflow.com/questions/40438390/how-to-put-outputs-side-by-side-in-shiny
          fluidRow(
            column(width=5, h3("Cell tissue image")),
            column(width=7, h3("Scatter plot of cell clusters"),
                   style = 'border-left: 1px solid')
          ),
          fluidRow(
            column(width=5, imageOutput("tissue")),
            column(width=7, plotOutput("scatter"),
                   style = 'border-left: 1px solid'),
            style = 'border-bottom: 1px solid'
          ),
          
          # Visualize 1 biomarker at a time? No scrolling.
          fluidRow(
            column(width=6, h3("Biomarker expression levels across clusters")),
            column(width=6, h3("Biomarker expression levels in individual cells"),
                   style = 'border-left: 1px solid')
          ),
          fluidRow(
            column(width=6, plotOutput("violinPlot")),
            column(width=6, plotOutput("featurePlot"),
                   style = 'border-left: 1px solid')
          )
          
          # Visualize 6 biomarkers at a time? Extra scrolling will be needed
          # h3("Violin plots of biomarker expression levels across each cluster"),
          # plotOutput("violinPlot"),
          # h3("Visualize biomarker expression levels in individual cells"),
          # plotOutput("featurePlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  current_file = reactive(input$file) # which data to display
  print("Checkpoint 1")
  
  in_file = reactive( paste0("spatialdata/", current_file()) )
  data = reactive({readRDS(in_file())})
  print("Checkpoint 2")
  
  # The data may not be in the right format.
  # If the assay is of type Assay5, it needs to be converted to Seurat.
  data_seurat = reactive({
    # Create the Seurat object
    seurat_obj = CreateSeuratObject(counts = data()@assays[[1]])
    
    # Transfer metadata and other components
    seurat_obj@meta.data = data()@meta.data
    seurat_obj@active.ident = data()@active.ident
    seurat_obj@graphs = data()@graphs
    seurat_obj@reductions = data()@reductions
    seurat_obj@images = data()@images
    seurat_obj@commands = data()@commands
    print("Metadata imported!")
    
    # This next step is vitally important. Can't draw FeaturePlots without a "data" layer!
    seurat_obj@assays$RNA$data = data()@assays$SCT$data
    
    # https://github.com/satijalab/seurat/discussions/8328
    print(Assays(seurat_obj))
    print(Layers(seurat_obj))
    
    return(seurat_obj)
  })
  print("Checkpoint 3")
  
  # Rownames of the file are the biomarkers
  # Actually, the data should be found in the first layer of "assays"
  # Layer is referred to by number in case the layer names diverge across input files
  markers <- reactive(rownames(data_seurat()))
  
  print("Checkpoint 4")
  
  # Display the tissue image
  # https://shiny.posit.co/r/articles/build/images/
  output$tissue = renderImage({
    print("Checkpoint 5")
    list(src="tissue_hires_image.png", width="100%")
  }, deleteFile = FALSE)
  
  # Select biomarker(s) in the left sidebar
  output$biomarkers = renderUI({
    print("Checkpoint 6")
    # print(unique(markers$gene))
    # print(sort(unique(markers()$cluster)))
    
    selectizeInput("biomarkers",
                   "Select up to 4 biomarkers to visualize expression levels:",
                   choices = markers(),
                   multiple = TRUE,
                   options = list(maxItems = 4)
                   )
    
  })
  
  rtype = reactive(input$red_type)
  # To be used for scatterplot and feature plots
  print("Checkpoint 7")
  
  # Show the individual clusters in the top right
  output$scatter = renderPlot({
    print("Checkpoint 8")
    if(rtype() == "umap"){
      print("UMAP")
      # display <- RunUMAP(data(), dims = 1:10) # holding up time!
      DimPlot(data_seurat(), reduction = "umap", pt.size = 4, label = TRUE, label.size = 7)
    } else if(rtype() == "tsne"){
      display <- RunTSNE(data_seurat(), dims = 1:10)
      DimPlot(display, reduction = "tsne", pt.size = 4, label = TRUE, label.size = 7)
    } else {
      print("Uh oh...")
    }
  })
  
  selected_markers = reactive({input$biomarkers}) # which biomarkers will be plotted?
  # We can use this list of selected markers to plot both 
  print("Checkpoint 9")

  # Show the violin plots in the middle
  output$violinPlot <- renderPlot({
      print("Checkpoint 10")
      print(selected_markers())
      #print( unique(filtered_markers()[filtered_markers()$cluster == 0,]$gene) )
      if(is.null(selected_markers())){
        print("Problem! (10)")
      } else {
        print("No problem. (10)")
        # 2 columns
        VlnPlot(data_seurat(),
                features = selected_markers(),
                ncol = 2
                )
      }
  })
  
  # Show the feature plots at the bottom
  # https://satijalab.org/seurat/reference/featureplot
  output$featurePlot = renderPlot({
    print("Checkpoint 11")
    
    selected <- selected_markers()
    if (is.null(selected) || length(selected) == 0) {
      print("Problem! (11)")
      return(NULL)  # Stop rendering if no markers are selected
    } else {
      print("No problem. (11)")
    }
    
    seurat_data <- data_seurat()
    if (is.null(seurat_data)) {
      print("data_seurat is NULL")
      return(NULL)  # Prevent further errors
    }
    
    assays_available <- Assays(seurat_data)
    print(assays_available)
    
    if (!"RNA" %in% assays_available) {
      print("RNA assay does not exist")
      return(NULL)
    }
    
    # print(rownames(seurat_data@assays$RNA))  # Safely access the RNA assay
    
    FeaturePlot(seurat_data,
                features = selected,
                pt.size = 0.5,
                ncol = 2,
                reduction = rtype())
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
