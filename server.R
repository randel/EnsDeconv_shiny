# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize=50000*1024^2)

#library(BiocManager)
#options(repos = BiocManager::repositories())
########Library Packages
# library(knitr)
# library(markdown)
# library(shiny)
# library(shinythemes)
# library(shinycustomloader)
# library(shinyhelper)
# library(foreach)
# library(tidyverse)
# library(EnsDeconv)
# library(sparseMatrixStats)
# #library(scran)
# library(ggpubr)
# library(DeconRNASeq)
# library(dplyr)
# library(formattable)
# library(shinysky)
# library(shinycssloaders)
# library(shinyWidgets)
# library(R.utils)


# load(file = "data/human_lengths.rda")
# load(file = "data/db.RData")

#source file
# source("methods/Deconv_meth_helper.R")
# source("methods/my_bisque.R")
# source("methods/my_music_construct.R")
# source("methods/my_music_utils.R")
# source("methods/cosg.R")
# source("pkg/DeconRNASeq.R")
# source("pkg/hspe.R")
shinyServer(function(input, output,session){
  observeEvent(input$chooseref, {
    updateTabsetPanel(inputId = "params", selected = input$chooseref)
  })
  
  observeEvent(input$choosepara, {
    updateTabsetPanel(inputId = "para_params", selected = input$choosepara)
  })
  
  
  
  ########### Update reference ###########
  refdata <- reactive({
    req(input$metaref)
    if(str_detect(input$metaref$datapath,"csv")){
      metaref <- read.csv(input$metaref$datapath,
                          header = TRUE)
    }else if(str_detect(input$metaref$datapath,"rds")){
      metaref <- readRDS(input$metaref$datapath)
    }else if(str_detect(input$metaref$datapath,"txt")){
      metaref <-  read.delim(input$metaref$datapath)
    }
  })
  
  
  observeEvent(refdata(), {
    updateSelectInput(session, "columnsref", choices=colnames(refdata()))
    updateSelectInput(session, "columnssample", choices=colnames(refdata()))
  })
  
  
  
  
  
  
  ###############CRM output################
  #########################################
  observe_helpers()
  dcInput <- eventReactive(input$dcupdate, {
    
    if (input$data_source == "custom") {
      req(input$bulk)
      # Load bulk data from uploaded file
      if (str_detect(input$bulk$datapath, "csv")) {
        to_deconv <- read.csv(input$bulk$datapath, header = TRUE)
        if (sum(sapply(to_deconv, is.numeric)) != ncol(to_deconv)) {
          to_deconv <- to_deconv[!duplicated(to_deconv[, 1]), ]
          bulkgene <- to_deconv[, 1]
          to_deconv <- as.matrix(to_deconv[, -1])
          rownames(to_deconv) <- bulkgene
          rm(bulkgene)
        }
      } else if (str_detect(input$bulk$datapath, "rds")) {
        to_deconv <- readRDS(input$bulk$datapath)
      } else if (str_detect(input$bulk$datapath, "txt")) {
        to_deconv <- read.delim(input$bulk$datapath)
        if (sum(sapply(to_deconv, is.numeric)) != ncol(to_deconv)) {
          to_deconv <- to_deconv[!duplicated(to_deconv[, 1]), ]
          bulkgene <- to_deconv[, 1]
          to_deconv <- as.matrix(to_deconv[, -1])
          rownames(to_deconv) <- bulkgene
          rm(bulkgene)
        }
      }


    } else if (input$data_source == "demo") {
      to_deconv <- readRDS("./data/to_deconv.rds")
    }
    
    
    if(input$chooseref == "custom"){
      metaref <- refdata() 
      # load reference data
      if(str_detect(input$ref$datapath,"csv")){
        ref <- read.csv(input$ref$datapath,
                        header = TRUE)
        if(sum(sapply(ref, is.numeric)) != ncol(ref)){
          ref <- ref[!duplicated(ref[,1]), ]
          refgene <-  ref[,1]
          ref <- as.matrix(ref[,-1])
          rownames(ref) <- refgene
        }
      }else if(str_detect(input$ref$datapath,"rds")){
        ref <- readRDS(input$ref$datapath)
      }else if(str_detect(input$ref$datapath,"txt")){
        ref <-  read.delim(input$ref$datapath)
        if(sum(sapply(ref, is.numeric)) != ncol(ref)){
          ref <- ref[!duplicated(ref[,1]), ]
          refgene <-  ref[,1]
          ref <- as.matrix(ref[,-1])
          rownames(ref) <- refgene
        }
      }
      
      refname <- input$columnsref
      colnames(metaref)[which(colnames(metaref) == refname)] = "deconv_clust" 
      refname <- input$columnssample
      colnames(metaref)[which(colnames(metaref) == refname)] = "SubjectName" 
    }else if(input$chooseref == "tissue"){
      
      ref_list <- readRDS(paste0("./data/",input$localtissue,"_sig.rds"))

      if(input$localtissue == "Brain"){
        sub_ref <- sub(" .*", "", input$brainReferences)
        ref_list <- ref_list[sub_ref]

        # Detect whether bulk data contains gene name or Ensembl gene IDs
        if(sum(unique(sapply(rownames(to_deconv), str_length)) != 15) !=0){
          library(biomaRt)
          mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
          ref_list <- lapply(ref_list, function(cur){
            gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                              values = rownames(cur$ref_matrix), mart= mart)
            gene_IDs = gene_IDs %>% filter(hgnc_symbol!="")
            
            cur$ref_matrix = cur$ref_matrix[gene_IDs$ensembl_gene_id,]
            rownames(cur$ref_matrix) = gene_IDs$hgnc_symbol
          })
          
        }
      }else if(input$localtissue == "Blood"){
        if(sum(unique(sapply(rownames(to_deconv), str_length)) != 15) !=0){
          library(biomaRt)
          mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
          ref_list <- lapply(ref_list, function(cur){
            gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                              values = rownames(cur$ref_matrix), mart= mart)
            gene_IDs = gene_IDs %>% filter(hgnc_symbol!="")
            
            cur$ref_matrix = cur$ref_matrix[gene_IDs$ensembl_gene_id,]
            rownames(cur$ref_matrix) = gene_IDs$hgnc_symbol
          })
          
        }
      }
      
      
    }
    # else{
    #   metaref <- readRDS(paste0("./data/meta_",input$localblood,".rds"))
    #   ref <- readRDS(paste0("./data/ref_",input$localblood,".rds"))
    #   # Detect whether bulk data contains gene name or Ensembl gene IDs
    #   if(sum(unique(sapply(rownames(to_deconv), str_length)) != 15) !=0){
    #     library(biomaRt)
    #     mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    #     gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
    #                       values = rownames(ref), mart= mart)
    #     gene_IDs = gene_IDs %>% filter(hgnc_symbol!="")
    #     
    #     ref = ref[gene_IDs$ensembl_gene_id,]
    #     rownames(ref) = gene_IDs$hgnc_symbol
    #   }
    # }
    time_limit <- 60*input$time_tolerance
    
    
    if(input$chooseref != "tissue"){
      metaref$SamplesNames = colnames(ref)
      
      ref_list = list()
      ref_list$"ref" = list()
      ref_list$"ref"$ref_matrix = as.matrix(ref)
      ref_list$"ref"$meta_ref = metaref
      
      params = get_params(data_type = input$datatype, data_name = names(testdata$ref_list), n_markers = input$nmrk,Marker.Method = input$mrk,TNormalization = input$norm,CNormalization =input$norm ,dmeths = input$Deconv,Scale = input$scale)
    }else{
      params = get_params(data_type = input$datatype, data_name = names(ref_list), n_markers = input$nmrk,Marker.Method = "none",TNormalization = "none",CNormalization ="none" ,dmeths = input$Deconv,Scale = input$scale)
    }
    

    if(nrow(params)>1){
      #params
      if(input$choosepara == "FALSE"){
        tryCatch({
          withTimeout({
            res <- EnsDeconv(count_bulk = as.matrix(to_deconv), ref_list = ref_list, ncore = 4, parallel_comp = F, params = params,inrshiny = T)
          },timeout = time_limit)
        }, TimeoutException = function(ex) {
          stop("Function execution timed out!")
        })
        
      }else{
        tryCatch({
          withTimeout({
            res <- EnsDeconv(count_bulk = as.matrix(to_deconv), ref_list = ref_list, ncore = input$ncore, parallel_comp = T, params = params)
          },timeout = time_limit)
        }, TimeoutException = function(ex) {
          stop("Function execution timed out!")
        })
        
      }
      
      
      res_p = lapply(res[[2]], function(x) x[["a"]][["p_hat"]][[1]][[1]])
      res_p[["Ensemble"]] = res[["EnsDeconv"]][["ensemble_p"]]
      
      #res_new <- list(res_p = res_p, case_ctrl = input$case_ctrl)
      #outt = lapply(res, function(x) x[["estimate"]])
      res_p[sapply(res_p, is.null)] <- NULL
      res_p
    }else{
      if(input$choosepara == "FALSE"){
        tryCatch({
          withTimeout({
            res <- gen_all_res_list_rshiny(count_bulk = as.matrix(to_deconv), meta_bulk = NULL, ref_list = ref_list, true_frac =NULL, outpath =NULL, ncore =4, parallel_comp = F, params = params)
          }, timeout = time_limit)
        }, TimeoutException = function(ex) {
          stop("Function execution timed out!")
        })
        
      }else{
        tryCatch({
          withTimeout({
            res <- gen_all_res_list(count_bulk = as.matrix(to_deconv), meta_bulk = NULL, ref_list = ref_list, true_frac =NULL, outpath =NULL, ncore =input$ncore, parallel_comp = T, params = params)          
          },timeout = time_limit)
        }, TimeoutException = function(ex) {
          stop("Function execution timed out!")
        })
      }
      ind = sapply(res, function(x){
        length(x[["a"]][["p_hat"]][[1]])
      })
      res = res[which(ind == 1)]
      res_p = lapply(res, function(x) x[["a"]][["p_hat"]][[1]][[1]])
      res_p[sapply(res_p, is.null)] <- NULL
      res_p
    }
    
    
  })
  
  
  output$plots <- renderPlot({
    
    res<- dcInput()
    # n <- length(res)
    # print(res)
    # plot_output_list <- lapply(1:n, function(i) {
    #   
    #   plotname <- names(res)[[i]]
    #   plotOutput(plotname, height = 580, width = 550)
    #   
    #   
    # })
    # 
    # 
    # do.call(tagList, plot_output_list)
    # df <- list()
    # for (i in 1:length(res)) {
    u <- intersect(rownames(res[[1]]),rownames(res[["Ensemble"]]))
    d <- data.frame(Subject = u,Method = "EnsDeconv",res[["Ensemble"]][u,]/rowSums(res[["Ensemble"]][u,]))
    df <-reshape2::melt(d,id.vars = c('Subject','Method'))
    # }
    # df <- bind_rows(df)
    names(df) <- c("Subject","Method","CellType","p_hat")
    
    ggboxplot(df,x = "CellType", y ="p_hat",color = "CellType")+ylab("EnsDeconv estimated cell type proportion")+xlab("Cell type")+color_palette("jco")
  
    
  })
  
  output$summaryText <- renderPrint({
    res <- dcInput()
    print(res[["Ensemble"]])
  })
  
  output$heatmap <- renderPlot({
    res<- dcInput()
    pheatmap::pheatmap(res[["Ensemble"]],cluster_rows = F, cluster_cols = F)
  })
  
  # Downloadable csv of selected dataset ----
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     "outcome.RData"
  #     
  #   },
  #   content = function(file) {
  #     resout = dcInput()
  #     resout = lapply(resout, as.data.frame)
  #     resout <- lapply(resout, function(x) rownames_to_column(x, "sample"))
  #     save(resout, file = file)
  #   }
  # )
  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$dfileType == "RData") {
        "results.RData"
      } else {
        "results.xlsx"
      }
    },
    content = function(file) {
      resout = dcInput()
      resout = lapply(resout, as.data.frame)
      resout <- lapply(resout, function(x) rownames_to_column(x, "sample"))
      names(resout)[names(resout) == "Ensemble"] <- "EnsDeconv"
      
      if ("EnsDeconv" %in% names(resout)) {
        resout <- c(EnsDeconv = list(resout[["EnsDeconv"]]), resout[names(resout) != "EnsDeconv"])
      }
      
      resout <- bind_rows(lapply(names(resout), function(name) {
        df <- resout[[name]]
        df$scenario <- name
        df
      }), .id = "Index")
      if (input$dfileType == "RData") {
        save(resout, file = file)
      } else if (input$dfileType == "Excel") {
        lapply(seq_along(resout), function(i) {
          write.xlsx(resout, file = file)
        })
      }
    }
  )
  
  
  data <- reactive({
    req(input$metabulk)
    if(str_detect(input$metabulk$datapath,"csv")){
      metabulk <- read.csv(input$metabulk$datapath,
                           header = TRUE)
    }else if(str_detect(input$metabulk$datapath,"rds")){
      metabulk <- readRDS(input$metabulk$datapath)
    }else if(str_detect(input$metabulk$datapath,"RData")){
      metabulk <- load(input$metabulk$datapath)
    }
  })
  
  # filtereddata <- eventReactive({
  #    input$update
  #    data()
  # },  {
  #    req(data())
  #    if(is.null(input$columns) || input$columns == "")
  #       data() else 
  #          data()[, colnames(data()) %in% input$columns]
  # })
  
  observeEvent(data(), {
    updateSelectInput(session, "columns", choices=colnames(data()))
  })
  
})