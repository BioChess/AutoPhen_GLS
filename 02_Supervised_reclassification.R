################################################################################
# 02 - Supervised Reclassification (Shiny App)
# Adapted from 02_Supervised_reclassification.qmd
# Author: Diego Vicente-Sastre
# Adapted for Shiny interactive windows
################################################################################

# == 1. Packages ===============================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, dplyr, lubridate,
  readr,
  ggplot2, ggpubr,
  leaflet, shiny, DT
)

# == 2. Helper functions =======================================================

my_theme <- theme_bw() +
  theme(panel.border     = element_rect(color = "black", fill = NA),
        plot.title       = element_text(size = 20),
        legend.title     = element_text(size = 14, face = "bold"),
        legend.text      = element_text(size = 14),
        axis.title.x     = element_text(size = 16),
        axis.text.x      = element_text(size = 14),
        axis.title.y     = element_text(size = 16),
        axis.text.y      = element_text(size = 14))

breed_color_map <- list(
  "B"    = "orange",
  "mig1" = "black",
  "mig2" = "grey",
  "W"    = "steelblue"
)

# Prepare track_df for visualisation
prep_track <- function(track_df) {
  track_df %>%
    filter(!is.na(x) & !is.na(y)) %>%
    mutate(
      col1 = case_when(
        Breed2 == "B"    ~ "orange",
        Breed2 == "mig1" ~ "black",
        Breed2 == "mig2" ~ "grey",
        Breed2 == "W"    ~ "steelblue"
      ),
      x    = round(x, 2),
      y    = round(y, 2),
      date = as.Date(date)
    ) %>%
    dplyr::select(x, y, date, any_of("Equinox"), Breed, Breed2, col1)
}

# Prepare seg_df with standard columns
prep_seg <- function(seg_df) {
  seg_df %>%
    mutate(
      f_day   = as.Date(f_day),
      l_day   = as.Date(l_day),
      segment = row_number()
    ) %>%
    dplyr::select(
      ID, State = state, Seg = segment, Breed,
      Start = f_day, End = l_day, Eq, Agreement, Check, Check_reason
    ) %>%
    mutate(Breed2 = case_when(
      Breed == "mig" & cumsum(replace_na(Breed == "W", FALSE)) > 0 &
        rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) == 0 ~ "mig2",
      Breed == "mig" & rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) > 0 ~ "mig1",
      TRUE ~ Breed
    )) %>%
    dplyr::select(-Breed)
}

# Join track with seg to get State and Seg for each position.
join_track_seg <- function(track, seg) {
  idx <- vapply(
    seq_len(nrow(track)),
    function(i) {
      d <- track$date[i]
      j <- which(seg$Start <= d & seg$End >= d)
      if (length(j) > 0) j[1] else NA_integer_
    },
    integer(1)
  )
  
  track %>%
    mutate(
      State = as.integer(seg$State[idx]),
      Seg   = as.integer(seg$Seg[idx])
    ) %>%
    left_join(
      seg %>% dplyr::select(-Breed2) %>% mutate(State = as.integer(State),
                                                Seg   = as.integer(Seg)),
      by = join_by(State, Seg)
    )
}

# == 3. Main Shiny app =========================================================

ui <- fluidPage(
  titlePanel("GLS Supervised Reclassification"),
  
  # -- Global sidebar (individual and step selection) --
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("📂 Data"),
      textInput("tracks_path", "Tracks CSV path:",
                value = "output/Phenology_tracks.csv"),
      textInput("segs_path",   "Segments CSV path:",
                value = "output/Segmentation_results.csv"),
      actionButton("load_data", "Load data", class = "btn-primary"),
      hr(),
      
      h4("👤 Individual"),
      uiOutput("id_selector"),
      actionButton("load_ind", "Load individual", class = "btn-info"),
      hr(),
      
      h4("🔄 Step"),
      radioButtons("step", "Reclassification mode:",
                   choices = c("By segment"  = "seg",
                               "By position" = "point"),
                   selected = "seg"),
      hr(),
      
      h4("💾 Save"),
      actionButton("save_ind",  "Save current individual"),
      actionButton("save_all",  "Save all and export CSV",
                   class = "btn-success"),
      hr(),
      verbatimTextOutput("status")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(id = "tabs",
                  # Tab 1: Track overview
                  tabPanel("Overview",
                           fluidRow(
                             column(6, plotOutput("plot_state",  height = 400)),
                             column(6, plotOutput("plot_breed2", height = 400))
                           )
                  ),
                  # Tab 2: Reclassification by segment
                  tabPanel("Reclassif. by segment",
                           fluidRow(
                             column(8,  leafletOutput("map_seg",   height = 560)),
                             column(4,
                                    h4("Editing panel"),
                                    selectInput("seg_breed2", "New phenological stage:",
                                                choices = names(breed_color_map)),
                                    actionButton("seg_update", "Update segment",
                                                 class = "btn-warning"),
                                    actionButton("seg_clear",  "Clear selection"),
                                    actionButton("seg_reset",  "Reset"),
                                    hr(),
                                    h5("Segment summary"),
                                    DTOutput("seg_table")
                             )
                           )
                  ),
                  # Tab 3: Reclassification by position
                  tabPanel("Reclassif. by position",
                           fluidRow(
                             column(8,  leafletOutput("map_pt",    height = 560)),
                             column(4,
                                    h4("Editing panel"),
                                    selectInput("pt_breed2", "New phenological stage:",
                                                choices = names(breed_color_map)),
                                    actionButton("pt_update", "Update points",
                                                 class = "btn-warning"),
                                    actionButton("pt_clear",  "Clear selection"),
                                    actionButton("pt_reset",  "Reset"),
                                    hr(),
                                    h5("Selected points"),
                                    DTOutput("pt_table")
                             )
                           )
                  ),
                  # Tab 4: Accumulated results
                  tabPanel("Results",
                           h4("Processed individuals"),
                           DTOutput("results_table"),
                           hr(),
                           h4("Corrected tracks (preview)"),
                           DTOutput("tracks_preview")
                  )
      )
    )
  )
)

# ==============================================================================
server <- function(input, output, session) {
  
  # ---- Reactive state --------------------------------------------------------
  rv <- reactiveValues(
    tracks_ph   = NULL,   # full tracks table
    seg_all     = NULL,   # full segments table
    ids         = NULL,   # vector of IDs
    # active individual
    track_raw   = NULL,
    seg_raw     = NULL,
    track_prep  = NULL,   # prepared track (with col1)
    seg_prep    = NULL,   # prepared segments
    track_seg   = NULL,   # track joined with seg (for segment reclassification)
    track_pt    = NULL,   # track for point reclassification (starts from track_seg)
    # selections
    sel_seg     = vector(),
    sel_pts     = vector(),
    # accumulated results
    segs_list   = list(),
    track_list  = list(),
    status_msg  = "Waiting for data..."
  )
  
  # ---- 3.1 Data loading ------------------------------------------------------
  observeEvent(input$load_data, {
    req(input$tracks_path, input$segs_path)
    tryCatch({
      # ID always as character (there are numeric IDs, e.g. 11001001,
      # and alphanumeric ones, e.g. BD120001)
      rv$tracks_ph <- read_csv(input$tracks_path, show_col_types = FALSE,
                               col_types = cols(ID = col_character()))
      rv$seg_all   <- read_csv(input$segs_path,   show_col_types = FALSE,
                               col_types = cols(ID = col_character()))
      rv$ids       <- as.character(unique(rv$tracks_ph$ID))
      rv$status_msg <- paste0("✅ Data loaded. ",
                              length(rv$ids), " individuals found.")
    }, error = function(e) {
      rv$status_msg <- paste0("❌ Error loading data: ", e$message)
    })
  })
  
  output$id_selector <- renderUI({
    req(rv$ids)
    selectInput("sel_id", "Select ID:", choices = rv$ids)
  })
  
  # ---- 3.2 Load individual ---------------------------------------------------
  # FIXED: wrapped in tryCatch so that an error does not close the app
  observeEvent(input$load_ind, {
    req(rv$tracks_ph, rv$seg_all, input$sel_id)
    id <- input$sel_id
    
    tryCatch({
      seg_raw   <- rv$seg_all   %>% filter(as.character(ID) == id)
      track_raw <- rv$tracks_ph %>% filter(as.character(ID) == id)
      
      # Fill NAs in Breed2
      if (any(is.na(track_raw$Breed2))) {
        track_raw <- track_raw %>%
          tidyr::fill(Breed,  .direction = "down") %>%
          tidyr::fill(Breed2, .direction = "down")
      }
      
      rv$track_raw  <- track_raw
      rv$seg_raw    <- seg_raw
      rv$track_prep <- prep_track(track_raw)
      rv$seg_prep   <- prep_seg(seg_raw)
      rv$track_seg  <- join_track_seg(rv$track_prep, rv$seg_prep)
      rv$track_pt   <- rv$track_seg
      rv$sel_seg    <- vector()
      rv$sel_pts    <- vector()
      rv$status_msg <- paste0("👤 Individual loaded: ", id)
    }, error = function(e) {
      rv$status_msg <- paste0("❌ Error loading individual: ", conditionMessage(e))
    })
  })
  
  # ---- 3.3 Overview (ggplot figures) -----------------------------------------
  make_plots <- reactive({
    req(rv$track_raw)
    track_p <- rv$track_raw %>%
      filter(!is.na(x) & !is.na(y)) %>%
      mutate(
        col1  = case_when(
          Breed2 == "B"    ~ "orange",
          Breed2 == "mig1" ~ "black",
          Breed2 == "mig2" ~ "grey",
          Breed2 == "W"    ~ "steelblue"),
        state = as.factor(state),
        lab   = case_when(
          Breed2 == "B"    ~ "Breeding",
          Breed2 == "mig1" ~ "Migration post",
          Breed2 == "mig2" ~ "Migration pre",
          Breed2 == "W"    ~ "Wintering")
      )
    list(
      p1 = ggplot(track_p) +
        geom_path(aes(x, y)) +
        geom_point(aes(x, y, color = state), size = 3) +
        xlab("Longitude") + ylab("Latitude") + my_theme +
        ggtitle("By state (segmentation)"),
      p2 = ggplot(track_p) +
        geom_path(aes(x, y)) +
        geom_point(aes(x, y, fill = col1), size = 3, shape = 21) +
        scale_fill_identity(
          guide  = "legend", name = "Classification",
          breaks = unique(track_p$col1),
          labels = unique(track_p$lab)) +
        xlab("Longitude") + ylab("Latitude") + my_theme +
        ggtitle("By phenological stage (Breed2)")
    )
  })
  
  output$plot_state  <- renderPlot({ make_plots()$p1 })
  output$plot_breed2 <- renderPlot({ make_plots()$p2 })
  
  # ---- 3.4 Segment map -------------------------------------------------------
  base_map_seg <- reactive({
    req(rv$track_seg)
    d <- rv$track_seg
    leaflet(d) %>%
      addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
      addCircleMarkers(lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addPolylines(lng = ~x, lat = ~y, color = "black", weight = 2)
  })
  
  output$map_seg <- renderLeaflet({ base_map_seg() })
  
  # Segment selection on click
  observeEvent(input$map_seg_marker_click, {
    req(rv$track_seg)
    cid <- as.numeric(input$map_seg_marker_click$id)
    clicked_seg <- rv$track_seg$Seg[cid]
    req(!is.na(clicked_seg))
    rv$sel_seg <- clicked_seg
    d <- rv$track_seg
    leafletProxy("map_seg") %>%
      clearMarkers() %>% clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addCircles(data = d[which(d$Seg == clicked_seg), ],
                 lat = ~y, lng = ~x, radius = 8,
                 color = "red", weight = 2, opacity = 1,
                 fill = FALSE, group = "sel")
  })
  
  # Update segment
  observeEvent(input$seg_update, {
    req(length(rv$sel_seg) > 0, rv$track_seg)
    new_d <- rv$track_seg
    rows  <- which(new_d$Seg %in% rv$sel_seg)
    new_d[rows, "Breed2"] <- input$seg_breed2
    new_d[rows, "col1"]   <- breed_color_map[[input$seg_breed2]]
    rv$track_seg <- new_d
    leafletProxy("map_seg") %>%
      clearMarkers() %>% clearGroup("sel") %>%
      addCircleMarkers(data = new_d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(new_d))) %>%
      addPolylines(data = new_d, lng = ~x, lat = ~y,
                   color = "black", weight = 2)
    # Propagate to track_pt so it starts from the result of the previous step
    rv$track_pt <- new_d
  })
  
  # Clear segment selection
  observeEvent(input$seg_clear, {
    req(rv$track_seg)
    rv$sel_seg <- vector()
    d <- rv$track_seg
    leafletProxy("map_seg") %>%
      clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d)))
  })
  
  # Reset segments
  observeEvent(input$seg_reset, {
    req(rv$track_prep, rv$seg_prep)
    rv$track_seg <- join_track_seg(rv$track_prep, rv$seg_prep)
    rv$sel_seg   <- vector()
    d <- rv$track_seg
    leafletProxy("map_seg") %>%
      clearMarkers() %>% clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addPolylines(data = d, lng = ~x, lat = ~y,
                   color = "black", weight = 2)
  })
  
  # Segment summary table
  output$seg_table <- renderDT({
    req(rv$track_seg)
    rv$track_seg %>%
      dplyr::select(Seg, Ph.stage = Breed2, Start, End, Eq, Agreement, Check, Check_reason) %>%
      distinct() %>% drop_na(Ph.stage) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ---- 3.5 Point map ---------------------------------------------------------
  base_map_pt <- reactive({
    req(rv$track_pt)
    d <- rv$track_pt
    leaflet(d) %>%
      addProviderTiles(providers$Esri.WorldGrayCanvas) %>%
      addCircleMarkers(lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addPolylines(lng = ~x, lat = ~y, color = "black", weight = 2)
  })
  
  output$map_pt <- renderLeaflet({ base_map_pt() })
  
  # Point selection on click (each click adds or removes the point)
  observeEvent(input$map_pt_marker_click, {
    req(rv$track_pt)
    cid <- as.numeric(input$map_pt_marker_click$id)
    # Toggle: remove if already selected, add otherwise
    if (cid %in% rv$sel_pts) {
      rv$sel_pts <- rv$sel_pts[rv$sel_pts != cid]
    } else {
      rv$sel_pts <- c(rv$sel_pts, cid)
    }
    d <- rv$track_pt
    leafletProxy("map_pt") %>%
      clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addCircles(data = d[rv$sel_pts, ],
                 lat = ~y, lng = ~x, radius = 8,
                 color = "red", weight = 2, opacity = 1,
                 fill = FALSE, group = "sel")
  })
  
  # Update points
  observeEvent(input$pt_update, {
    req(length(rv$sel_pts) > 0, rv$track_pt)
    new_d <- rv$track_pt
    new_d[rv$sel_pts, "Breed2"] <- input$pt_breed2
    new_d[rv$sel_pts, "col1"]   <- breed_color_map[[input$pt_breed2]]
    rv$track_pt <- new_d
    leafletProxy("map_pt") %>%
      clearMarkers() %>% clearGroup("sel") %>%
      addCircleMarkers(data = new_d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(new_d))) %>%
      addPolylines(data = new_d, lng = ~x, lat = ~y,
                   color = "black", weight = 2)
  })
  
  # Clear point selection
  observeEvent(input$pt_clear, {
    req(rv$track_pt)
    rv$sel_pts <- vector()
    d <- rv$track_pt
    leafletProxy("map_pt") %>%
      clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d)))
  })
  
  # Reset points
  observeEvent(input$pt_reset, {
    req(rv$track_seg)
    rv$track_pt <- rv$track_seg   # go back to the result of the segment step
    rv$sel_pts  <- vector()
    d <- rv$track_pt
    leafletProxy("map_pt") %>%
      clearMarkers() %>% clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addPolylines(data = d, lng = ~x, lat = ~y,
                   color = "black", weight = 2)
  })
  
  # Table of selected points / full track
  output$pt_table <- renderDT({
    req(rv$track_pt)
    d <- rv$track_pt
    show <- if (length(rv$sel_pts) > 0) d[rv$sel_pts, ] else d
    datatable(show, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ---- 3.6 Save current individual -------------------------------------------
  observeEvent(input$save_ind, {
    req(rv$track_pt, rv$track_seg, input$sel_id)
    id <- input$sel_id
    
    # final track: fill NAs if any
    track_final <- rv$track_pt
    if (any(is.na(track_final$Breed2))) {
      track_final <- track_final %>%
        tidyr::fill(Breed2, .direction = "down")
    }
    
    # final segments derived from the reclassified track_seg
    seg_final <- rv$track_seg %>%
      dplyr::select(Breed2, State, Seg, Start, End, Eq, Agreement, Check, Check_reason) %>%
      distinct() %>% drop_na(Breed2)
    
    rv$track_list[[id]] <- track_final
    rv$segs_list[[id]]  <- seg_final
    rv$status_msg <- paste0("✅ Individual '", id, "' saved. ",
                            "Total saved: ", length(rv$track_list))
  })
  
  # ---- 3.7 Export all --------------------------------------------------------
  # CHANGED: bind_rows(.id = "ID") instead of plyr::ldply (avoids depending on plyr)
  observeEvent(input$save_all, {
    req(length(rv$track_list) > 0)
    track_all <- dplyr::bind_rows(rv$track_list, .id = "ID")
    seg_all2  <- dplyr::bind_rows(rv$segs_list,  .id = "ID")
    
    # Saved in the same folder as the input files (output/ or output/full/)
    out_dir <- dirname(input$tracks_path)
    write_csv(track_all, file.path(out_dir, "Phenology_tracks_corrected.csv"))
    write_csv(seg_all2,  file.path(out_dir, "Segmentation_results_corrected.csv"))
    
    rv$status_msg <- paste0("💾 Exported: ",
                            nrow(track_all), " positions / ",
                            nrow(seg_all2),  " segments saved.")
  })
  
  # ---- 3.8 Results tab -------------------------------------------------------
  output$results_table <- renderDT({
    req(length(rv$track_list) > 0)
    data.frame(ID = names(rv$track_list),
               n_points = sapply(rv$track_list, nrow)) %>%
      datatable(options = list(pageLength = 20))
  })
  
  output$tracks_preview <- renderDT({
    req(length(rv$track_list) > 0)
    dplyr::bind_rows(rv$track_list, .id = "ID") %>%
      head(200) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ---- Status ----------------------------------------------------------------
  output$status <- renderText({ rv$status_msg })
}

# == 4. Launch the app =========================================================
shinyApp(ui = ui, server = server)