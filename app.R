library(shiny)
library(shinyjs)
library(Seurat)
library(ggplot2)
library(tools)
# ================================
# UI for a single Seurat path module
# ================================
pathModuleUI <- function(id) {
  ns <- NS(id)
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
      uiOutput(ns("numeric_filter_ui")),  # UI output for numeric filters and violins
      uiOutput(ns("subset_ui"))
    )
  )
}

# ================================
# Server logic for a single module
# ================================
pathModuleServer <- function(id, remove_callback) {
  moduleServer(id, function(input, output, session) {
    
    local_rv <- reactiveValues(
      obj = NULL,
      path_value = "",
      obj_loaded = FALSE,
      numeric_filter_range = NULL  # Store numeric filter ranges here
    )

    # Enable/disable inputs based on load status
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
    
    # Load the Seurat object from provided RDS path
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
        
        isolate({ local_rv$numeric_filter_range <- list() })  # Reset numeric filter ranges
        
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
    
    # Reset module to initial state
    observeEvent(input$reset, {
      local_rv$obj <- NULL
      local_rv$obj_loaded <- FALSE
      updateTextInput(session, "path", value = "")
      
      isolate({ local_rv$numeric_filter_range <- list() })
    })
    
    # Remove module UI and server instance
    observeEvent(input$remove, {
      remove_callback(id)
    })
    
    # UI to select meta.data columns (exclude "CB")
    output$meta_ui <- renderUI({
      req(local_rv$obj_loaded)
      obj <- local_rv$obj
      
      cols <- colnames(obj@meta.data)
      choices <- setdiff(cols, "CB")
      
      selectInput(
        session$ns("meta_col"),
        "Select meta.data column(s)",
        choices = choices,
        multiple = TRUE
      )
    })
    
    # UI for filtering character meta columns: checkbox + rename text input
    output$subset_ui <- renderUI({
      req(local_rv$obj_loaded, input$meta_col)
      obj <- local_rv$obj
      
      meta_cols <- input$meta_col
      # Only non-numeric columns here
      non_numeric_cols <- meta_cols[!sapply(obj@meta.data[meta_cols], is.numeric)]
      if (length(non_numeric_cols) == 0) return(NULL)
      
      tagList(
        lapply(non_numeric_cols, function(col) {
          vals <- sort(unique(obj@meta.data[[col]]))
          tagList(
            tags$h5(col),
            lapply(vals, function(val) {
              safe <- make.names(paste(col, val, sep = "__"))
              fluidRow(
                column(
                  4,
                  checkboxInput(
                    session$ns(paste0("keep_", safe)),
                    label = val,
                    value = TRUE
                  )
                ),
                column(
                  8,
                  textInput(
                    session$ns(paste0("rename_", safe)),
                    label = NULL,
                    value = val
                  )
                )
              )
            })
          )
        })
      )
    })
    
    
    # UI for numeric filters: violin plot + range slider
    output$numeric_filter_ui <- renderUI({
      req(local_rv$obj_loaded, input$meta_col)
      obj <- local_rv$obj
      
      meta_cols <- input$meta_col
      numeric_cols <- meta_cols[sapply(obj@meta.data[meta_cols], is.numeric)]
      if (length(numeric_cols) == 0) return(NULL)
      
      tagList(
        lapply(numeric_cols, function(col) {
          vals <- obj@meta.data[[col]]
          min_val <- floor(min(vals, na.rm = TRUE))
          max_val <- ceiling(max(vals, na.rm = TRUE))
          
          tagList(
            tags$h5(col),
            plotOutput(session$ns(paste0("violin_", col)), height = "100px"),
            fluidRow(
              column(6,
                     numericInput(
                       session$ns(paste0("min_", col)),
                       label = paste0("Min ", col),
                       value = min_val,
                       min = min_val,
                       max = max_val,
                       step = 1
                     )
              ),
              column(6,
                     numericInput(
                       session$ns(paste0("max_", col)),
                       label = paste0("Max ", col),
                       value = max_val,
                       min = min_val,
                       max = max_val,
                       step = 1
                     )
              )
            ),
            tags$hr()
          )
        })
      )
    })
    
    observe({
      req(local_rv$obj_loaded, input$meta_col)
      obj <- local_rv$obj
      
      meta_cols <- input$meta_col
      numeric_cols <- meta_cols[sapply(obj@meta.data[meta_cols], is.numeric)]
      
      lapply(numeric_cols, function(col) {
        vals <- obj@meta.data[[col]]
        min_val <- floor(min(vals, na.rm = TRUE))
        max_val <- ceiling(max(vals, na.rm = TRUE))
        
        slider_val <- local_rv$numeric_filter_range[[col]]
        if (!is.null(slider_val)) {
          updateNumericInput(
            session,
            paste0("min_", col),
            value = slider_val[1],
            min = min_val,
            max = max_val
          )
          updateNumericInput(
            session,
            paste0("max_", col),
            value = slider_val[2],
            min = min_val,
            max = max_val
          )
        }
      })
    })
    
    # Render violin plots for numeric columns
    observe({
      req(local_rv$obj_loaded, input$meta_col)
      obj <- local_rv$obj
      meta_cols <- input$meta_col
      numeric_cols <- meta_cols[sapply(obj@meta.data[meta_cols], is.numeric)]
      
      for (col in numeric_cols) {
        local({
          col_local <- col
          output[[paste0("violin_", col_local)]] <- renderPlot({
            df <- data.frame(value = obj@meta.data[[col_local]])
            ggplot(df, aes(x = "", y = value)) +
              geom_violin(trim = FALSE, fill = "lightblue", color = "gray40") +
              geom_jitter(width = 0.3, alpha = 0.5, size = 0.1) +
              coord_flip() +
              theme_minimal() +
              labs(y = col_local, x = NULL) +
              theme(axis.text.x = element_text(), axis.ticks.x = element_line()) +
              scale_y_continuous(breaks = scales::pretty_breaks(n = 6))
          })
        })
      }
    })
    
    # Update numeric_filter_range when sliders change
    observe({
      req(local_rv$obj_loaded, input$meta_col)
      meta_cols <- input$meta_col
      obj <- local_rv$obj
      
      numeric_cols <- meta_cols[sapply(obj@meta.data[meta_cols], is.numeric)]
      
      for (col in numeric_cols) {
        local({
          col_local <- col
          min_id <- paste0("min_", col_local)
          max_id <- paste0("max_", col_local)
          
          observeEvent({
            input[[min_id]]
          }, {
            isolate({
              if (is.null(local_rv$numeric_filter_range)) local_rv$numeric_filter_range <- list()
              current_max <- local_rv$numeric_filter_range[[col_local]][2] %||% Inf
              min_val <- input[[min_id]]
              if (!is.null(min_val) && min_val <= current_max) {
                local_rv$numeric_filter_range[[col_local]][1] <- min_val
              }
            })
          }, ignoreInit = TRUE)
          
          observeEvent({
            input[[max_id]]
          }, {
            isolate({
              if (is.null(local_rv$numeric_filter_range)) local_rv$numeric_filter_range <- list()
              current_min <- local_rv$numeric_filter_range[[col_local]][1] %||% -Inf
              max_val <- input[[max_id]]
              if (!is.null(max_val) && max_val >= current_min) {
                local_rv$numeric_filter_range[[col_local]][2] <- max_val
              }
            })
          }, ignoreInit = TRUE)
        })
      }
    })
    
    
    # Return reactive values for parent server to use
    list(
      obj = reactive({ local_rv$obj }),
      path = reactive({ local_rv$path_value }),
      obj_loaded = reactive({ local_rv$obj_loaded }),
      meta_col = reactive({ input$meta_col }),
      numeric_filter_range = reactive({ local_rv$numeric_filter_range }),
      get_choices = reactive({
        req(local_rv$obj_loaded, input$meta_col)
        
        obj <- local_rv$obj
        out <- list()
        
        # Only character meta columns have choices (checkbox + rename)
        for (col in input$meta_col) {
          if (!is.numeric(obj@meta.data[[col]])) {
            vals <- sort(unique(obj@meta.data[[col]]))
            col_choices <- list()
            
            for (val in vals) {
              safe <- make.names(paste(col, val, sep = "__"))
              keep_id <- paste0("keep_", safe)
              rename_id <- paste0("rename_", safe)
              
              if (!is.null(input[[keep_id]]) && input[[keep_id]]) {
                col_choices[[val]] <- input[[rename_id]]
              }
            }
            
            if (length(col_choices) > 0) {
              out[[col]] <- col_choices
            }
          }
        }
        
        out
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
  tags$hr(),
  div(id = "path_container"),
  fluidRow(
    style = "margin-left: 0px;",
    actionButton("add_path", "➕ Add RDS file"),
    actionButton("generate_new", "Generate New RDS File", class = "btn-success")
  ),
  br(),br(),br()
)

# ================================
# Main Server
# ================================
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    path_ids = character(),
    path_modules = list(),
    counter = 0,
    initialized = FALSE,
    compiler_log = NULL
  )
  filter_summary_rv <- reactiveValues(
    compiler_log = NULL
  )
  # Add new path module dynamically
  observeEvent(input$add_path, {
    new_id <- paste0("sample_", rv$counter + 1)
    rv$counter <- rv$counter + 1
    
    isolate({
      rv$path_ids <- c(rv$path_ids, new_id)
    })
    
    insertUI(
      selector = "#path_container",
      where = "beforeEnd",
      ui = pathModuleUI(new_id)
    )
    
    remove_callback <- function(module_id) {
      isolate({
        rv$path_ids <- setdiff(rv$path_ids, module_id)
        rv$path_modules[[module_id]] <- NULL
      })
      removeUI(selector = paste0("#", module_id, "_wrapper"))
    }
    
    module <- pathModuleServer(new_id, remove_callback)
    
    isolate({
      rv$path_modules[[new_id]] <- module
    })
  })
  
  # Automatically add first module on app start
  observe({
    if (!rv$initialized) {
      rv$initialized <- TRUE
      shinyjs::click("add_path")
    }
  })
  
  # Show summary modal before generating new object
  observeEvent(input$generate_new, {
    req(length(rv$path_ids) > 0)
    
    compiler_log <- list()
    
    summary_blocks <- lapply(rv$path_ids, function(id) {
      module <- rv$path_modules[[id]]
      if (is.null(module) || !module$obj_loaded()) return(NULL)
      
      path_val <- module$path()
      meta_cols <- module$meta_col()
      choices <- module$get_choices()
      numeric_ranges <- module$numeric_filter_range()
      numeric_ranges <- if (is.null(numeric_ranges)) list() else numeric_ranges
      
      # Log for compiler
      compiler_log[[id]] <<- list(
        path = path_val,
        numeric_map = numeric_ranges,
        meta_map = choices
      )
      
      meta_blocks <- list()
      
      if (length(numeric_ranges) > 0) {
        meta_blocks <- c(meta_blocks, lapply(names(numeric_ranges), function(col) {
          rng <- numeric_ranges[[col]]
          tags$div(
            br(),
            tags$h5(paste("Meta column:", col)),
            tags$p(paste0("&nbsp;&nbsp;&nbsp;&nbsp;Range: ", rng[1], " to ", rng[2]))
          )
        }))
      }
      
      if (length(choices) > 0) {
        meta_blocks <- c(meta_blocks, lapply(names(choices), function(col) {
          rename_map <- choices[[col]]
          if (length(rename_map) == 0) {
            tags$p(strong(col), ": no values selected")
          } else {
            tags$div(
              br(),
              tags$h5(paste("Meta column:", col)),
              tags$table(
                class = "table table-condensed",
                tags$thead(
                  tags$tr(
                    tags$th("Original"),
                    tags$th(""),
                    tags$th("Renamed To")
                  )
                ),
                tags$tbody(
                  lapply(names(rename_map), function(old) {
                    tags$tr(
                      tags$td(old),
                      tags$td("→"),
                      tags$td(rename_map[[old]])
                    )
                  })
                )
              )
            )
          }
        }))
      }
      
      if (length(meta_blocks) == 0) {
        meta_blocks <- list(tags$p(em("No meta.data filters selected.")))
      }
      
      tagList(
        tags$hr(),
        tags$h4(tags$b("Input RDS File: "), path_val),
        meta_blocks
      )
    })
    
    # Clean NULLs
    summary_blocks <- Filter(Negate(is.null), summary_blocks)
    compiler_log <- Filter(Negate(is.null), compiler_log)
    
    filter_summary_rv$compiler_log <- compiler_log
    
    if (length(summary_blocks) == 0) {
      showModal(modalDialog(
        title = "No Data Loaded",
        "Please load at least one RDS file before generating.",
        easyClose = TRUE,
        footer = modalButton("Close")
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
  
  # Show modal to configure saving options after confirm
  observeEvent(input$confirm_generate, {
    removeModal()
    showModal(modalDialog(
      title = "Saving RDS file",
      size = "l",
      
      textInput("new_meta_name", "Meta.data column name in new RDS file after processing",
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
  
  # Handle final saving
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
      paste("Subsetting, Renaming and Merging Datasets..."),
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
        meta_cols <- module$meta_col()
        choices <- module$get_choices()
        numeric_ranges <- module$numeric_filter_range()
        
        if (is.null(meta_cols) || (length(choices) == 0 && length(numeric_ranges) == 0)) next
        
        cells_keep_list <- list()
        
        # Filter by character meta columns
        if (length(choices) > 0) {
          for (col in names(choices)) {
            keep_vals <- names(choices[[col]])
            cells_keep_list[[col]] <- rownames(obj@meta.data)[obj@meta.data[[col]] %in% keep_vals]
          }
        }
        
        # Filter by numeric meta columns (range)
        if (length(numeric_ranges) > 0) {
          for (col in names(numeric_ranges)) {
            rng <- numeric_ranges[[col]]
            cells_keep_list[[col]] <- rownames(obj@meta.data)[
              obj@meta.data[[col]] >= rng[1] & obj@meta.data[[col]] <= rng[2]
            ]
          }
        }
        
        # Intersect all filtering criteria
        if (length(cells_keep_list) == 0) next
        final_cells <- Reduce(intersect, cells_keep_list)
        obj_sub <- subset(obj, cells = final_cells)
        
        # Rename clusters in character meta columns
        for (col in names(choices)) {
          new_labels <- choices[[col]][obj_sub@meta.data[[col]]]
          obj_sub@meta.data[[col]] <- unname(as.character(new_labels))
        }
        
        # Merge logic based on mode
        if (merge_mode == "path_is_sample") {
          # Use path as sample identity
          if (length(add_ids) == 0) {
            obj_sub@meta.data[[new_meta_name]] <- rep(file_path_sans_ext(basename(module$path())), ncol(obj_sub))
            merged_obj <- obj_sub
          } else {
            obj_sub@meta.data[[new_meta_name]] <- rep(file_path_sans_ext(basename(module$path())), ncol(obj_sub))
            merged_obj <- merge(merged_obj, obj_sub)
          }
          add_ids <- c(add_ids, id)
          
        } else if (merge_mode == "rename_based") {
          # Rename based on cluster renaming (for all meta columns)
          if (length(add_ids) == 0) {
            new_meta_col_values <- paste0(
              do.call(paste, c(obj_sub@meta.data[names(choices)], sep = "_"))
            )
            obj_sub@meta.data[[new_meta_name]] <- new_meta_col_values
            merged_obj <- obj_sub
          } else {
            new_meta_col_values <- paste0(
              do.call(paste, c(obj_sub@meta.data[names(choices)], sep = "_"))
            )
            obj_sub@meta.data[[new_meta_name]] <- new_meta_col_values
            merged_obj <- merge(merged_obj, obj_sub)
          }
          add_ids <- c(add_ids, id)
        }
      }
      
      if (length(add_ids) == 0) {
        removeModal()
        showModal(modalDialog(
          title = "No Cells Selected",
          "No cells passed the filtering criteria in any loaded RDS file.",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        return()
      }
      
      showModal(modalDialog(
        title = "Generating New RDS file",
        paste("Saving to:", input$save_path),
        footer = NULL,
        easyClose = FALSE
      ))
      
      commands_list <- list()
      for (id in rv$path_ids) {
        module <- rv$path_modules[[id]]
        if (is.null(module) || !module$obj_loaded()) next
        
        obj <- module$obj()
        if (!is.null(obj@commands)) {
          commands_list[[id]] <- obj@commands
        }
      }

      filter_log <- filter_summary_rv$compiler_log
      
      merged_obj@commands <- filter_log
      saveRDS(merged_obj, input$save_path)
      removeModal()
      
      showModal(modalDialog(
        title = "Success",
        paste("New RDS file saved successfully to:", input$save_path),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("Error during processing:\n", e$message),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })
}

# ================================
# Run the app
# ================================
shinyApp(ui, server)
