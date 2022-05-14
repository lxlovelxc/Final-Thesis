#核心功能版

library(shiny)


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Please upload the features.tsv.gz", accept = ".gz",
              placeholder = "No file selected"),
      fileInput("file2", "Please upload the barcodes.tsv.gz", accept = ".gz",
              placeholder = "No file selected"),
      fileInput("file3","Please upload the matrix.mtx.gz",accept = ".gz",
              placeholder = "No file selected"),
      actionButton('go',"Go",icon = icon('play-circle')),
      hr(),
      sliderInput(inputId = "Range1",
                label = "Range of nFeature_RNA",
                min = 200,
                max = 10000,
                value = c(200,2500),
                width = "200px"),
      helpText("Please select the range to filter nFeature_RNA"),
      hr(),
      numericInput(inputId = "Range2",
                label = "Range of percent.mt",
                value = 5,
                width = "200px"),
      helpText("Please input the upper bound to filter high percent.mt cell"),
      hr(),
      numericInput(inputId = "PC_select",
                label = "Select diminsions you want (How many PCs you want to see)",
                value = "10",
                width = "200px"
      ),
      hr(),
      numericInput(inputId = "Resolutions",
                label = "Resolutions of FindCluster: ",
                value = "0.5",
                step = "0.1",
                width = "200px"),
      helpText("Format example: 1.0"),
      hr(),
      textInput(inputId = "genes",
                label = "genes of interested",
                width = "200px"
      ),
      helpText("name genes you're interested,separated by tab "),
    ),
    mainPanel(
      verbatimTextOutput('result'),
      verbatimTextOutput('result1'),
      verbatimTextOutput('result2'),
      verbatimTextOutput('Range1'),
      verbatimTextOutput('Range2'),
      verbatimTextOutput('PC_select'),
      verbatimTextOutput('Resolutions'),
      #textOutput("Genes"),
      plotOutput('UMAP'),
    ),
  ),
)

server <- function(input, output) {
  options(shiny.maxRequestSize=100*1024^2)
  destDir <- 'C:/Users/luoxi/Desktop/ne'
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
    input$genes
  })
  output$Genes <-renderPrint(input$genes)
    output$result <- renderPrint({
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
  
    output$result1 <- renderPrint({
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
    
    output$result2 <- renderPrint({
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
    
    output$UMAP <- renderPlot({
      library(dplyr)
      library(Seurat)
      setwd("C:/Users/luoxi/Desktop")
      
      observeEvent(input$go,{
        cat("you have success!\n")
        rna_data <- Read10X(data.dir = 'C:/Users/luoxi/Desktop/ne')
        rna <- CreateSeuratObject(counts = rna_data,min.cells = 3, min.features = 200)
        #使用从MT-作为线粒体基因集
        rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-") 
        # Visualize QC metrics as a violin plot
        p1<-VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        cat('p1 need to be print now')
        print(p1)
        plot1 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        cat("p2 now")
        print(plot1)
        print(plot2)
        rna <- subset(rna, subset = nFeature_RNA > Range1()[1] & nFeature_RNA < Range1()[2] & percent.mt < Range2())#######这个需要确定
        rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
        rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
        all.genes <- rownames(rna)
        rna <- ScaleData(rna, features = all.genes)
        rna <- RunPCA(rna, features = VariableFeatures(object = rna))# dimensionality reduction using PCA
        p3 <- VizDimLoadings(rna, dims = 1:2, reduction = "pca")
        cat("p3 now")
        print(p3)
        p4<-DimPlot(rna, reduction = "pca")
        print(p4)
        p5<-DimHeatmap(rna, dims = 1, cells = 500, balanced = TRUE)
        print(p5)
        #rna <- JackStraw(rna, num.replicate = 100)
        #rna <- ScoreJackStraw(rna, dims = 1:20)
        #p6<-JackStrawPlot(rna, dims = 1:15)
        #print(p6)
        rna <- FindNeighbors(rna, dims = 1:PC_select())#Clustering using shared nearest neighbor(SNN)
        rna <- FindClusters(rna, resolution = Resolutions())
        rna <- RunUMAP(rna, dims = 1:PC_select()) 
        #DimPlot(rna, reduction = "umap",label = TRUE)        
        UMAP<-DimPlot(rna, reduction = "umap")
        print(UMAP)
        cat("you have done successful!")
        # find all markers of cluster 2
        cluster2.markers <- FindMarkers(rna, ident.1 = 2, min.pct = 0.25)
        head(cluster2.markers, n = 5)
        # find all markers distinguishing cluster 5 from clusters 0 and 3
        cluster5.markers <- FindMarkers(rna, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
        head(cluster5.markers, n = 5)
        # find markers for every cluster compared to all remaining cells, report only the positive
        # ones
        rna.markers <- FindAllMarkers(rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        rna.markers %>%
          group_by(cluster) %>%
          slice_max(n = 2, order_by = avg_log2FC)
        cluster0.markers <- FindMarkers(rna, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
        rna.markers %>%
          group_by(cluster) %>%
          top_n(n = 10, wt = avg_log2FC) -> top10
        heat_map <- DoHeatmap(rna, features = top10$gene) + NoLegend()
        print(heat_map)
        #p4 <- VlnPlot(rna, features = c("MS4A1", "CD79A"))
        #print(p4)
        #p5<- FeaturePlot(rna, features = c("MS4A1", "CD79A"))
        #print(p5)
        
      })
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
