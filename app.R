#最终版
library(shiny)
library(shinydashboard)

header <- dashboardHeader(title = "scRNA Sequence Analysis",#标题
                          titleWidth = "300",
                          disable = FALSE)

sidebar <- dashboardSidebar(#侧边栏
  sidebarMenu(
      menuItem("Upload",tabName = "Upload"),
      menuItem("Adjust",tabName = "Adjust",icon = icon("th")),
      menuItem("UMAP",tabName = "UMAP",icon = icon("th")))
)

Upload <- fluidRow(
  p("The file name must be the ",strong("SAME")," to the following.",style = "color:red;font-family:Impact,Charcoal,sanserif;"),
  hr(),
  fileInput("file1", "Please upload the features.tsv.gz",
            placeholder = "No file selected"),
  fileInput("file2", "Please upload the barcodes.tsv.gz",
            placeholder = "No file selected"),
  fileInput("file3","Please upload the matrix.mtx.gz",
            placeholder = "No file selected"),
  verbatimTextOutput('result'),
  verbatimTextOutput('result1'),
  verbatimTextOutput('result2'),
  actionButton('go',"Go",icon = icon('play-circle')),#实现按下这个Go按钮，开始计算UMAP图
)

Adjust <- fluidRow(
  box(hr(),
      sliderInput(inputId = "Range1",
                  label = "Range of nFeature_RNA",
                  min = 200,
                  max = 10000,
                  value = c(200,2500),
                  width = "200px"),
      helpText("Please select the range to filter nFeature_RNA"),
      status = "primary",br(),
      plotOutput(outputId = "before_filter"),
      plotOutput(outputId = "after_filter"),),
  
  box(hr(),
      numericInput(inputId = "Range2",
                   label = "Range of percent.mt",
                   value = 5,
                   width = "200px"),
      helpText("Please input the upper bound to filter high percent.mt cell"),
      status = "primary",br(),
      plotOutput("Loadings"),
      plotOutput("Heatmap")),
  box(hr(),
      numericInput(inputId = "PC_select",
                   label = "Select diminsions you want (How many PCs you want to see)",
                   value = "10",
                   width = "200px"),
      status = "primary",br(),
      plotOutput(outputId = "Elbow")),
) 

UMAP <- fluidRow(
  box(hr(),
      numericInput(inputId = "Resolutions",
                    label = "Resolutions of FindCluster: ",
                    value = "0.5",
                    step = "0.1",
                    width = "200px"),
       helpText("Format example: 1.0"),
       status = "warning",br(),
       hr(),
       downloadButton("save", "UMAP"),
       hr(),
       plotOutput(outputId = 'UMAP'),
       plotOutput(outputId = "tSNE"),
       plotOutput(outputId = "PCA")
      ),
  
  box(hr(),
      textInput(inputId = "genes",
                label = "genes of interested",
                value = ("MS4A1"),
                width = "200px"),
      helpText("name genes you're interested,separated by tab "),
      status = "warning",br(),
      actionButton('Search',"Search",icon = icon('play-circle')),
      #plotOutput("heat_map"),
      plotOutput(outputId = 'gene'))
)

body <- dashboardBody(#主界面
  tabItems(
    tabItem(tabName = "Upload",Upload),
    tabItem(tabName = "Adjust",Adjust),
    tabItem(tabName = "UMAP",UMAP))
)


shinyApp(
  ui = dashboardPage(skin = "purple",
    header, sidebar, body),
  
  server <- function(input, output) {
    library(dplyr)
    library(Seurat)
    setwd("C:/Users/luoxi/Desktop")
    options(shiny.maxRequestSize=100*1024^2)#100MB
    destDir <- 'C:/Users/luoxi/Desktop/ne'#这是一个本地新建的文件夹，将用户上传的表达矩阵先放到这个文件夹里
    Range1 <- reactive({
      cbind(input$Range1[1],input$Range1[2])
    })
    Range2 <- reactive({
      input$Range2
    })
    PC_select <- reactive({
      input$PC_select
    })
    Resolutions <- reactive({
      input$Resolutions
    })
    Genes <- reactive({
      c(input$genes)
    })
    
    output$result <- renderPrint({#实现把上传的feature文件放到ne那个文件夹里
      fea <- input$file1
      
      if (is.null(fea)) {
        cat("NOT FILE\n")
        return(FALSE)
      }
      cat("Reading file:", fea$name, "\n")
      cat("size:", fea$size, " Bytes, type:", fea$type, "\n")
      if (dir.exists(destDir)){
        cat("Copying file to:", destDir,"\n")
        result <- file.copy( fea$datapath,
                             file.path(destDir, fea$name) )
      } else {
        result <- FALSE
      }
      result
    })
    
    output$result1 <- renderPrint({#实现把上传的barcodes文件放到ne那个文件夹里
      bar <- input$file2
      
      if (is.null(bar)) {
        cat("NOT FILE\n")
        return(FALSE)
      }
      cat("Reading file:", bar$name, "\n")
      cat("size:", bar$size, " Bytes, type:", bar$type, "\n")
      if (dir.exists(destDir)){
        cat("Copying file to:", destDir,"\n")
        result1 <- file.copy( bar$datapath,
                              file.path(destDir, bar$name) )
      } else {
        result1 <- FALSE
      }
      result1
    })
    
    output$result2 <- renderPrint({#实现把上传的matrix文件放到ne那个文件夹里
      mat <- input$file3
      
      if (is.null(mat)) {
        cat("NOT FILE\n")
        return(FALSE)
      }
      cat("Reading file:", mat$name, "\n")
      cat("size:", mat$size, " Bytes, type:", mat$type, "\n")
      if (dir.exists(destDir)){
        cat("Copying file to:", destDir,"\n")
        result2 <- file.copy( mat$datapath,
                              file.path(destDir, mat$name) )
      } else {
        result2 <- FALSE
      }
      result2
    })
      
    observeEvent(input$go,{
      cat("you have success!\n")
      rna_data <- Read10X(data.dir = 'C:/Users/luoxi/Desktop/ne')
      rna <- CreateSeuratObject(counts = rna_data,min.cells = 3, min.features = 200)
      rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-") 
      VlnPlot <- VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      output$VlnPlot <- renderPlot(VlnPlot)
      plot1 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      before_filter<-plot1 +plot2
      output$before_filter <- renderPlot(before_filter)
          rna <- subset(rna, subset = nFeature_RNA > Range1()[1] & nFeature_RNA < Range1()[2] & percent.mt <  Range2())
      plot3 <- FeatureScatter(rna,feature1 = "nCount_RNA",feature2 = "percent.mt")
      plot4 <- FeatureScatter(rna,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
      after_filter <-plot3+plot4
      output$after_filter <- renderPlot(after_filter)
      rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
      rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
      all.genes <- rownames(rna)
      rna <- ScaleData(rna, features = all.genes)
      rna <- RunPCA(rna, features = VariableFeatures(object = rna))
      #Investigate the intrinsic dimensionality of the data using an elbow plot:
      Elbow <- ElbowPlot(rna)
      output$Elbow <- renderPlot(Elbow)
      Loadings <- VizDimLoadings(rna,dims = 1:2,reduction = "pca")
      output$Loadings <- renderPlot(Loadings)
      Heatmap <- DimHeatmap(rna,dims = 1:2,cells = 100,balanced = TRUE)
      output$Heatmap <- renderPlot(Heatmap)
      rna <- FindNeighbors(rna, dims = 1:PC_select())
      rna <- FindClusters(rna, resolution = Resolutions())
      PCA <- DimPlot(rna,reduction = "pca")
      output$PCA <- renderPlot(PCA)
      
      rna <- RunUMAP(rna, dims = 1:PC_select())
      UMAP <- DimPlot(rna, reduction = "umap",label = TRUE)
      output$UMAP <- renderPlot(UMAP)
      rna <- RunTSNE(rna)
      tSNE<- DimPlot(rna,reduction = "tsne")
      output$tSNE <- renderPlot(tSNE)
      
      
      observeEvent(input$Search,{
        gene_on_UMAP <- FeaturePlot(rna, features = c(Genes()))
        output$gene <- renderPlot(gene_on_UMAP)
      })
      
      output$save <- downloadHandler(
        file = "UMAP" ,
        content = function(file) {
          png(file = UMAP)
        })
    })
})