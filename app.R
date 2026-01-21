library(shiny)
library(shinyjs)
library(Seurat)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Seurat Composer (ID-based)"),
  actionButton("add_path", "➕ Add Path"),
  tags$hr(),
  div(id = "path_container"),
  actionButton("generate_new", "Generate New RDS File", class = "btn-success")
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    paths = character(),
    seurat = list(),
    paths_values = list(),
    counter = 0
  )
  
  register_path <- function(id) {
    
    # Render the UI block for one path
    output[[paste0("path_ui_", id)]] <- renderUI({
      obj_loaded <- !is.null(rv$seurat[[id]])
      
      val <- if (obj_loaded) {
        paste0("Loaded: ", isolate(rv$paths_values[[id]] %||% ""))
      } else {
        isolate(rv$paths_values[[id]] %||% "")
      }
      
      wellPanel(
        h4(paste("Path", id)),
        fluidRow(
          column(
            6,
            textInput(
              paste0("path_", id),
              "Seurat RDS Path",
              value = val,
              placeholder = "/path/to/seurat.rds"
            )
          ),
          column(
            6,
            br(),
            div(style = "display: flex; gap: 5px; align-items: center;",
                actionButton(paste0("load_", id), "Load", disabled = obj_loaded),
                if (obj_loaded) {
                  actionButton(paste0("reset_", id), "Reset", class = "btn-warning")
                } else {
                  NULL
                },
                if (length(rv$paths) > 1) {
                  actionButton(paste0("remove_", id), "Delete", class = "btn-danger")
                } else {
                  NULL
                }
            )
          )
        ),
        if (obj_loaded) uiOutput(paste0("meta_ui_", id)),
        if (obj_loaded) uiOutput(paste0("subset_ui_", id))
      )
    })
    
    # After UI rendered, enable or disable the textInput dynamically
    observe({
      obj_loaded <- !is.null(rv$seurat[[id]])
      if (obj_loaded) {
        shinyjs::disable(paste0("path_", id))
      } else {
        shinyjs::enable(paste0("path_", id))
      }
    })
    
    # Load button with tryCatch
    observeEvent(input[[paste0("load_", id)]], {
      path <- input[[paste0("path_", id)]]
      req(path)
      
      showModal(modalDialog(title = "Loading", paste("Reading:", path), footer = NULL, easyClose = FALSE))
      
      obj <- NULL
      success <- FALSE
      tryCatch({
        obj <- readRDS(path)
        success <- TRUE
      }, error = function(e) {
        showModal(modalDialog(title = "Error",
                              paste("Failed to load Seurat object:\n", e$message),
                              easyClose = TRUE,
                              footer = modalButton("Close")))
      })
      
      if (success) {
        removeModal()
        rv$seurat[[id]] <- obj
        isolate({
          rv$paths_values[[id]] <- path
        })
      } else {
        removeModal()
      }
    }, ignoreInit = TRUE)
    
    # Reset button clears loaded data and keeps path editable
    observeEvent(input[[paste0("reset_", id)]], {
      rv$seurat[[id]] <- NULL
      rv$paths_values[[id]] <- NULL
    }, ignoreInit = TRUE)
    
    # Remove button deletes UI block and clears data
    observeEvent(input[[paste0("remove_", id)]], {
      isolate({
        rv$paths <- setdiff(rv$paths, id)
        rv$seurat[[id]] <- NULL
        rv$paths_values[[id]] <- NULL
      })
      removeUI(selector = paste0("#path_ui_", id))
    }, ignoreInit = TRUE)
    
    # Meta column selection UI (only if data loaded)
    output[[paste0("meta_ui_", id)]] <- renderUI({
      req(!is.null(rv$seurat[[id]]))
      obj <- rv$seurat[[id]]
      non_num <- names(obj@meta.data)[!sapply(obj@meta.data, is.numeric)]
      choices <- setdiff(non_num, "CB")
      selectInput(paste0("meta_col_", id), "Select meta.data column", choices = choices)
    })
    
    # Subset and rename UI based on selected meta column
    output[[paste0("subset_ui_", id)]] <- renderUI({
      req(!is.null(rv$seurat[[id]]))
      obj <- rv$seurat[[id]]
      col <- input[[paste0("meta_col_", id)]]
      req(col, col %in% colnames(obj@meta.data))
      
      vals <- sort(unique(obj@meta.data[[col]]))
      
      tagList(
        lapply(vals, function(val) {
          safe <- make.names(val)
          fluidRow(
            column(4,
                   checkboxInput(paste0("keep_", id, "_", safe), label = val, value = TRUE)
            ),
            column(8,
                   textInput(paste0("rename_", id, "_", safe), label = NULL, value = val)
            )
          )
        })
      )
    })
  }
  
  # Add new path block
  observeEvent(input$add_path, {
    new_id <- paste0("p", isolate(rv$counter) + 1)
    rv$counter <- isolate(rv$counter) + 1
    
    isolate({
      rv$paths <- c(rv$paths, new_id)
    })
    
    insertUI(selector = "#path_container",
             where = "beforeEnd",
             ui = uiOutput(paste0("path_ui_", new_id)))
    
    register_path(new_id)
  })
  
  # Initial path block at app start
  observeEvent(rv$counter, {
    if (rv$counter == 0) {
      new_id <- "p1"
      rv$counter <- 1
      rv$paths <- c(new_id)
      
      insertUI(selector = "#path_container",
               where = "beforeEnd",
               ui = uiOutput(paste0("path_ui_", new_id)))
      
      register_path(new_id)
    }
  }, once = TRUE)
}

shinyApp(ui, server)
