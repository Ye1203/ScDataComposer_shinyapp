library(shiny)
library(shinyjs)
library(Seurat)

# ================================
# UI for a single Seurat path module
# ================================
pathModuleUI <- function(id) {
  ns <- NS(id)
  
  # Wrap the whole module in a div with unique ID for removal
  div(
    id = paste0(id, "_wrapper"),
    wellPanel(
      h4(paste("RDS File", gsub(".*_", "", id))),
      
      fluidRow(
        column(
          6,
          textInput(
            ns("path"),
            "RDS File Path",
            value = "",
            placeholder = "/path/to/rds_file.rds"
          )
        ),
        column(
          6,
          br(),
          div(
            style = "display: flex; gap: 5px; align-items: center;",
            actionButton(ns("load"), "Load"),
            actionButton(ns("reset"), "Reset", class = "btn-warning"),
            actionButton(ns("remove"), "Delete", class = "btn-danger")
          )
        )
      ),
      
      uiOutput(ns("meta_ui")),
      uiOutput(ns("subset_ui"))
    )
  )
}

# ================================
# Server logic for a single module
# ================================
pathModuleServer <- function(id, remove_callback) {
  moduleServer(id, function(input, output, session) {
    
    # Local reactive state for this module only
    local_rv <- reactiveValues(
      obj = NULL,
      path_value = "",
      obj_loaded = FALSE
    )
    
    # Enable/disable inputs based on load state
    observe({
      if (local_rv$obj_loaded) {
        shinyjs::disable("path")
        shinyjs::disable("load")
        shinyjs::enable("reset")
      } else {
        shinyjs::enable("path")
        shinyjs::enable("load")
        shinyjs::disable("reset")
      }
    })
    
    # Load Seurat object from RDS
    observeEvent(input$load, {
      req(input$path)
      
      showModal(modalDialog(
        title = "Loading",
        paste("Reading:", input$path),
        footer = NULL,
        easyClose = FALSE
      ))
      
      tryCatch({
        obj <- readRDS(input$path)
        local_rv$obj <- obj
        local_rv$path_value <- input$path
        local_rv$obj_loaded <- TRUE
        removeModal()
      }, error = function(e) {
        removeModal()
        showModal(modalDialog(
          title = "Error",
          paste("Failed to load RDS file:\n", e$message),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      })
    })
    
    # Reset module state
    observeEvent(input$reset, {
      local_rv$obj <- NULL
      local_rv$obj_loaded <- FALSE
      updateTextInput(session, "path", value = "")
    })
    
    # Remove this module from UI and server
    observeEvent(input$remove, {
      remove_callback(id)
    })
    
    # UI for selecting a metadata column
    output$meta_ui <- renderUI({
      req(local_rv$obj_loaded)
      obj <- local_rv$obj
      
      # Keep only non-numeric metadata columns
      non_num <- names(obj@meta.data)[!sapply(obj@meta.data, is.numeric)]
      choices <- setdiff(non_num, "CB")
      
      selectInput(session$ns("meta_col"), "Select meta.data column", choices = choices)
    })
    
    # UI for subsetting and renaming cluster labels
    output$subset_ui <- renderUI({
      req(local_rv$obj_loaded, input$meta_col)
      obj <- local_rv$obj
      col <- input$meta_col
      
      if (is.null(col) || !col %in% colnames(obj@meta.data)) {
        return(NULL)
      }
      
      vals <- sort(unique(obj@meta.data[[col]]))
      
      tagList(
        lapply(vals, function(val) {
          safe <- make.names(val)
          fluidRow(
            column(
              4,
              checkboxInput(session$ns(paste0("keep_", safe)), label = val, value = TRUE)
            ),
            column(
              8,
              textInput(session$ns(paste0("rename_", safe)), label = NULL, value = val)
            )
          )
        })
      )
    })
    
    # Return reactive accessors to parent server
    list(
      obj = reactive({ local_rv$obj }),
      path = reactive({ local_rv$path_value }),
      obj_loaded = reactive({ local_rv$obj_loaded }),
      meta_col = reactive({ input$meta_col }),
      get_choices = reactive({
        req(local_rv$obj_loaded, input$meta_col)
        obj <- local_rv$obj
        col <- input$meta_col
        vals <- sort(unique(obj@meta.data[[col]]))
        
        choices <- list()
        for (val in vals) {
          safe <- make.names(val)
          keep_id <- paste0("keep_", safe)
          rename_id <- paste0("rename_", safe)
          
          if (!is.null(input[[keep_id]]) && input[[keep_id]]) {
            choices[[val]] <- input[[rename_id]]
          }
        }
        choices
      })
    )
  })
}

# ================================
# Main UI
# ================================
ui <- fluidPage(
  useShinyjs(),
  titlePanel("scRDS_DataCompiler"),
  actionButton("add_path", "➕ Add RDS file"),
  tags$hr(),
  div(id = "path_container"),
  actionButton("generate_new", "Generate New RDS File", class = "btn-success"),
  br(),br(),br()
)

# ================================
# Main Server
# ================================
server <- function(input, output, session) {
  
  # Reactive values to track modules
  rv <- reactiveValues(
    path_ids = character(),
    path_modules = list(),
    counter = 0,
    initialized = FALSE
  )
  
  # Add a new path module dynamically
  observeEvent(input$add_path, {
    new_id <- paste0("sample_", rv$counter + 1)
    rv$counter <- rv$counter + 1
    
    # Store the module ID
    isolate({
      rv$path_ids <- c(rv$path_ids, new_id)
    })
    
    # Insert the module UI
    insertUI(
      selector = "#path_container",
      where = "beforeEnd",
      ui = pathModuleUI(new_id)
    )
    
    # Define the removal callback function
    # This will be called when the module's delete button is clicked
    remove_callback <- function(module_id) {
      # Remove from reactive values
      isolate({
        rv$path_ids <- setdiff(rv$path_ids, module_id)
        rv$path_modules[[module_id]] <- NULL
      })
      
      # Remove the UI element
      # Use the wrapper ID we defined in pathModuleUI
      removeUI(
        selector = paste0("#", module_id, "_wrapper")
      )
    }
    
    # Initialize the server module
    module <- pathModuleServer(new_id, remove_callback)
    
    # Store the module reference
    isolate({
      rv$path_modules[[new_id]] <- module
    })
  })
  
  # Create the first module on app start (only once)
  # Using a flag to ensure it runs only once
  observe({
    if (!rv$initialized) {
      rv$initialized <- TRUE
      shinyjs::click("add_path")
    }
  })
  
  # Summary before generating new object
  observeEvent(input$generate_new, {
    req(length(rv$path_ids) > 0)
    
    compiler_log <- list()
    
    summary_blocks <- lapply(rv$path_ids, function(id) {
      module <- rv$path_modules[[id]]
      if (is.null(module) || !module$obj_loaded()) return(NULL)
      
      path_val <- module$path()
      meta_col <- module$meta_col()
      choices  <- module$get_choices()
      
      compiler_log[[id]] <<- list(
        path       = path_val,
        meta_col   = meta_col,
        rename_map = if (length(choices) > 0) choices else NULL
      )
      
      if (is.null(meta_col) || length(choices) == 0) {
        cluster_info <- tags$p(em("No clusters selected"))
      } else {
        cluster_info <- lapply(names(choices), function(val) {
          tags$tr(
            tags$td(val),
            tags$td("→"),
            tags$td(choices[[val]])
          )
        })
        
        cluster_info <- tags$table(
          class = "table table-condensed",
          tags$thead(
            tags$tr(
              tags$th("Original"),
              tags$th(""),
              tags$th("Renamed To")
            )
          ),
          tags$tbody(cluster_info)
        )
      }
      
      tagList(
        tags$hr(),
        tags$h4(tags$b("Input RDS File: "), path_val),
        tags$p(tags$b("Meta column: "), meta_col),
        cluster_info
      )
    })
    
    compiler_log <- Filter(Negate(is.null), compiler_log)
    rv$compiler_log <- compiler_log
    
    summary_blocks <- Filter(Negate(is.null), summary_blocks)
    
    if (length(summary_blocks) == 0) {
      showModal(modalDialog(
        title = "No Data Loaded",
        "Please load at least one RDS file before generating.",
        easyClose = TRUE
      ))
      return()
    }
    
    showModal(modalDialog(
      title = "Summary Before Generating New Object",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_generate", "Confirm", class = "btn btn-success")
      ),
      do.call(tagList, summary_blocks)
    ))
  })
  
  # Confirm generation of new RDS file
  observeEvent(input$confirm_generate, {
    removeModal()
    showModal(modalDialog(
      title = "Saving RDS file",
      size = "l",
      
      textInput("new_meta_name", "the meta.data column name in new RDS file after processing",
                value = "SampleID", placeholder = "SampleID_xxxx"),
      radioButtons(
        "merge_mode",
        "How should Input RDS files be interpreted?",
        choices = c(
          "Mode 1: Input RDS files contain multiple samples" = "rename_based",
          "Mode 2: Each Input RDS file corresponds to a single sample" = "path_is_sample"
        )
      ),
      uiOutput("merge_description"),
      hr(),
      textInput("save_path", "Enter path to save new RDS file",
                value ="",
                placeholder = "path/to/new_file.rds",
                width = "100%"),
      uiOutput("save_path_error"),
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_saving", "Save", class = "btn btn-success")
      )
    ))
  })
  
  output$merge_description <- renderUI({
    req(input$merge_mode)
    
    if (input$merge_mode == "rename_based") {
      HTML(
        "
      - Each RDS file may contain multiple samples.<br>
      - The new metadata column will be created based on your renamed cluster labels after subsetting.<br><br>
      <b>Example:</b><br>
      &nbsp;&nbsp;&nbsp;&bull; Input RDS files 1: sample1, sample2<br>
      &nbsp;&nbsp;&nbsp;&bull; Input RDS files 2: sample3, sample4<br><br>
      <b>After merging, new metadata column will contain:</b><br>
      sample1, sample2, sample3, sample4"
      )
      
    } else {
      HTML(
        "
      - Each RDS file represents one independent sample.<br>
      - The new metadata column will be created using Path names.<br>
      - Cluster renaming does NOT affect sample identity.<br><br>
      <b>Example:</b><br>
      &nbsp;&nbsp;&nbsp;&bull; Input RDS files 1 (sample1): biological cluster1<br>
      &nbsp;&nbsp;&nbsp;&bull; Input RDS files 2 (sample2): biological cluster1<br><br>
      <b>After merging, new metadata column will contain:</b><br>
      sample1, sample2"
      )
    }
  })
  
  output$save_path_error <- renderUI({
    if (is.null(input$save_path) || input$save_path == "") return(NULL)
    if (!grepl("\\.rds$", input$save_path, ignore.case = TRUE)) {
      tags$div(style = "color: red;", "Error: path must end with .rds.")
    } else if (!dir.exists(dirname(input$save_path))) {
      tags$div(style = "color: orange;", "Warning: directory does not exist, folders will be generated automatically.")
    } else {
      tags$div(style = "color: green;", "Correct Path, Click \"Save\" to start processing.")
    }
  })
  
  # Confirm saving the new RDS file
  observeEvent(input$confirm_saving, {
    req(input$save_path)
    removeModal()
    
    if (!grepl("\\.rds$", input$save_path, ignore.case = TRUE)) {
      showModal(modalDialog(
        title = "Saving Path Input Error",
        "Saving Path must end with .rds",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return()
    }
    
    showModal(modalDialog(
      title = "Generating New RDS file",
      paste("Subseting, Renaming and Merging Datasets..."),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      dir_path <- dirname(input$save_path)
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
      }
      new_meta_name <- input$new_meta_name
      merge_mode <- input$merge_mode
      
      processed_list <- list()
      add_ids <- c()
      
      for (id in rv$path_ids) {
        
        module <- rv$path_modules[[id]]
        if (is.null(module) || !module$obj_loaded()) next
        
        obj <- module$obj()
        meta_col <- module$meta_col()
        choices <- module$get_choices()
        
        # Skip if nothing selected
        if (is.null(meta_col) || length(choices) == 0) next
        
        # Keep only selected groups
        keep_vals <- names(choices)
        cells_keep <- rownames(obj@meta.data)[obj@meta.data[[meta_col]] %in% keep_vals]
        obj_sub <- subset(obj, cells = cells_keep)
        
        # Rename cluster labels
        new_labels <- choices[obj_sub@meta.data[[meta_col]]]
        obj_sub@meta.data[, meta_col] <- as.character(unname(new_labels))
        
        if (merge_mode == "rename_based") {
          obj_sub[[new_meta_name]] <- obj_sub@meta.data[, meta_col]
        } else if (merge_mode == "path_is_sample") {
          obj_sub[[new_meta_name]] <- id
        } else {
          obj_sub[[new_meta_name]] <- id
        }
        
        # Store processed object
        processed_list[[id]] <- obj_sub
        add_ids <- c(add_ids, id)
      }
      
      
      # Ensure at least one dataset remains
      if (length(processed_list) == 0) {
        stop("No valid datasets available after subsetting.")
      }
      # Merge Seurat objects
      if (length(processed_list) == 1) {
        merged_obj <- processed_list[[1]]
      } else {
        common_features <- Reduce(intersect, lapply(processed_list, rownames))
        processed_list <- lapply(processed_list, function(obj) {
          subset(obj, features = common_features)
        })
        processed_list <- lapply(processed_list, function(obj) {
          if ("SCT" %in% Assays(obj)) {
            obj[["SCT"]] <- NULL
          }
          obj
        })
        merged_obj <- merge(
          x = processed_list[[1]],
          y = processed_list[-1],
          collapse = FALSE,
          add.cell.ids = add_ids,
          project = "SeuratObj"
        )
      }
      rna_assay <- merged_obj[["RNA"]]
      # Join layers (important for Seurat v5 objects)
      if (inherits(rna_assay, "Assay5")) {
        layer_names <- Layers(rna_assay)
        if (length(layer_names) > 1) {
          message("Joining layers in RNA assay...")
          merged_obj[["RNA"]] <- JoinLayers(rna_assay)
        } else {
          message("RNA assay has a single layer. No join needed.")
        }
      } else {
        message("RNA assay is v4-style Assay. JoinLayers not applicable.")
      }
      
      merged_obj@commands$scRDS_DataCompiler <- list(
        timestamp = Sys.time(),
        paths     = rv$compiler_log
      )
      
      showModal(modalDialog(
        title = "Generating New RDS file",
        paste("Saving to:", input$save_path),
        footer = NULL,
        easyClose = FALSE
      ))
      
      # Save final object
      saveRDS(merged_obj, file = input$save_path)
      gc()
      removeModal()
      showModal(modalDialog(
        title = "Success",
        paste("New RDS file saved to:", input$save_path),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("Failed to generate RDS file:\n", e$message),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })
  
}

shinyApp(ui, server)