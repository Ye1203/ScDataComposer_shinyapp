library(shiny)
library(Seurat)
library(shinyjs)
# =========================
# UI
# =========================
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Seurat Composer"),
  
  actionButton("add_path", "➕ Add Path"),
  
  tags$hr(),
  
  div(
    id = "path_container",
    uiOutput("path_ui_1")
  )
)

# =========================
# SERVER
# =========================
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    n = 1,
    seurat = list()
  )
  
  # ---- helper: load seurat (eventReactive style) ----
  load_seurat <- function(i) {
    eventReactive(input[[paste0("load_", i)]], {
      
      path <- input[[paste0("path_", i)]]
      req(path)
      
      showModal(
        modalDialog(
          title = "Loading Seurat Object",
          paste("Loading:", path),
          footer = NULL,
          easyClose = FALSE
        )
      )
      obj <- readRDS(path)
      removeModal()
      updateTextInput(
        session,
        paste0("path_", i),
        label = NULL,
        value = paste0("Loaded: ", path)
      )
      shinyjs::disable(paste0("path_", i))
      obj
    }, ignoreInit = TRUE)
  }
  
  
  seurat_loaded <- list()
  
  # ---- UI for a single path block ----
  make_path_ui <- function(i) {
    output[[paste0("path_ui_", i)]] <- renderUI({
      tagList(
        wellPanel(
          h4(paste("Path", i)),
          
          fluidRow(
            column(
              8,
              textInput(
                paste0("path_", i),
                "Seurat RDS Path",
                placeholder = "/path/to/seurat.rds"
              )
            ),
            column(
              2,
              br(),
              actionButton(paste0("load_", i), "Load")
            ),
            column(
              2,
              if (i > 1) {
                br()
                actionButton(
                  paste0("remove_", i),
                  "Delete",
                  class = "btn-danger"
                )
              }
            )
          ),
          
          uiOutput(paste0("meta_select_ui_", i)),
          uiOutput(paste0("subset_rename_ui_", i))
        )
      )
    })
  }
  
  # initialize first block
  make_path_ui(1)
  seurat_loaded[[1]] <- load_seurat(1)
  
  # ---- add new path ----
  observeEvent(input$add_path, {
    rv$n <- rv$n + 1
    i <- rv$n
    
    insertUI(
      selector = "#path_container",
      where = "beforeEnd",
      ui = uiOutput(paste0("path_ui_", i))
    )
    
    make_path_ui(i)
    seurat_loaded[[i]] <<- load_seurat(i)
  })
  
  # ---- remove path ----
  observe({
    lapply(2:rv$n, function(i) {
      observeEvent(input[[paste0("remove_", i)]], {
        removeUI(selector = paste0("#path_ui_", i))
        rv$seurat[[i]] <- NULL
      }, ignoreInit = TRUE)
    })
  })
  
  # ---- per-path logic ----
  observe({
    lapply(1:rv$n, function(i) {
      
      # ---- meta.data column selector (only after load) ----
      output[[paste0("meta_select_ui_", i)]] <- renderUI({
        obj <- seurat_loaded[[i]]()
        req(obj)
        
        selectInput(
          paste0("meta_col_", i),
          "Select meta.data column",
          choices = mdcols <- setdiff(names(sapply(obj@meta.data, function(x) !is.numeric(x))[sapply(obj@meta.data, function(x) !is.numeric(x))]), "CB")
        )
      })
      
      # ---- unique values + checkbox + rename ----
      output[[paste0("subset_rename_ui_", i)]] <- renderUI({
        obj <- seurat_loaded[[i]]()
        col <- input[[paste0("meta_col_", i)]]
        
        req(
          obj,
          col,
          col %in% colnames(obj@meta.data),
          nrow(obj@meta.data) > 0
        )
        
        vals <- sort(unique(obj@meta.data[[col]]))
        
        tagList(
          lapply(vals, function(val) {
            safe_val <- make.names(val)
            
            fluidRow(
              column(
                4,
                checkboxInput(
                  paste0("keep_", i, "_", safe_val),
                  label = val,
                  value = TRUE
                )
              ),
              column(
                8,
                textInput(
                  paste0("rename_", i, "_", safe_val),
                  label = NULL,
                  value = val
                )
              )
            )
          })
        )
      })
    })
  })
}

# =========================
# RUN APP
# =========================
shinyApp(ui = ui, server = server)
