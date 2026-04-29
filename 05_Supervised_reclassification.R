################################################################################
# 05 - Supervised Reclassification (Shiny App)
# Adaptado desde 05_Supervised_reclassification.qmd
# Author: Diego Vicente-Sastre
# Adapted for Shiny interactive windows
################################################################################

# == 1. Paquetes ===============================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, dplyr, lubridate,
  readr,
  ggplot2, ggpubr,
  leaflet, shiny, DT
)

# == 2. Funciones auxiliares ===================================================

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

# Prepara track_df para visualización
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
    dplyr::select(x, y, date, Breed, Breed2, col1)
}

# Prepara seg_df con columnas estándar
prep_seg <- function(seg_df) {
  seg_df %>%
    mutate(
      f_day   = as.Date(f_day),
      l_day   = as.Date(l_day),
      segment = row_number()
    ) %>%
    dplyr::select(
      ID, State = state, Seg = segment, Breed,
      Start = f_day, End = l_day, Eq, Confidence, Check
    ) %>%
    mutate(Breed2 = case_when(
      Breed == "mig" & cumsum(replace_na(Breed == "W", FALSE)) > 0 &
        rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) == 0 ~ "mig2",
      Breed == "mig" & rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) > 0 ~ "mig1",
      TRUE ~ Breed
    )) %>%
    dplyr::select(-Breed)
}

# Une track con seg para obtener State y Seg por punto
join_track_seg <- function(track, seg) {
  track %>%
    mutate(
      State = as.integer(map_int(list(date = date), ~ {
        idx <- which(seg$Start <= ..1 & seg$End >= ..1)
        if (length(idx) > 0) seg$State[idx[1]] else NA_character_
      })),
      Seg = as.integer(map_int(list(date = date), ~ {
        idx <- which(seg$Start <= ..1 & seg$End >= ..1)
        if (length(idx) > 0) seg$Seg[idx[1]] else NA_character_
      }))
    ) %>%
    left_join(seg, by = join_by(Breed2, State, Seg))
}

# == 3. Shiny App principal ====================================================

ui <- fluidPage(
  titlePanel("GLS Supervised Reclassification"),

  # -- Barra lateral global (selección de individuo y paso) --
  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("📂 Datos"),
      textInput("tracks_path", "Ruta tracks CSV:",
                value = "output/Phenology_tracks.csv"),
      textInput("segs_path",   "Ruta segmentos CSV:",
                value = "output/Segmentation_results.csv"),
      actionButton("load_data", "Cargar datos", class = "btn-primary"),
      hr(),

      h4("👤 Individuo"),
      uiOutput("id_selector"),
      actionButton("load_ind", "Cargar individuo", class = "btn-info"),
      hr(),

      h4("🔄 Paso"),
      radioButtons("step", "Modo de reclasificación:",
                   choices = c("Por segmento" = "seg",
                               "Por posición"  = "point"),
                   selected = "seg"),
      hr(),

      h4("💾 Guardar"),
      actionButton("save_ind",  "Guardar individuo actual"),
      actionButton("save_all",  "Guardar todos y exportar CSV",
                   class = "btn-success"),
      hr(),
      verbatimTextOutput("status")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(id = "tabs",
        # Tab 1: Vista general del track
        tabPanel("Vista general",
          fluidRow(
            column(6, plotOutput("plot_state",  height = 400)),
            column(6, plotOutput("plot_breed2", height = 400))
          )
        ),
        # Tab 2: Reclasificación por segmento
        tabPanel("Reclasif. por segmento",
          fluidRow(
            column(8,  leafletOutput("map_seg",   height = 560)),
            column(4,
              h4("Panel de edición"),
              selectInput("seg_breed2", "Nuevo estadio fenológico:",
                          choices = names(breed_color_map)),
              actionButton("seg_update", "Actualizar segmento",
                           class = "btn-warning"),
              actionButton("seg_clear",  "Limpiar selección"),
              actionButton("seg_reset",  "Reset"),
              hr(),
              h5("Resumen de segmentos"),
              DTOutput("seg_table")
            )
          )
        ),
        # Tab 3: Reclasificación por puntos
        tabPanel("Reclasif. por posición",
          fluidRow(
            column(8,  leafletOutput("map_pt",    height = 560)),
            column(4,
              h4("Panel de edición"),
              selectInput("pt_breed2", "Nuevo estadio fenológico:",
                          choices = names(breed_color_map)),
              actionButton("pt_update", "Actualizar puntos",
                           class = "btn-warning"),
              actionButton("pt_clear",  "Limpiar selección"),
              actionButton("pt_reset",  "Reset"),
              hr(),
              h5("Puntos seleccionados"),
              DTOutput("pt_table")
            )
          )
        ),
        # Tab 4: Resultados acumulados
        tabPanel("Resultados",
          h4("Individuos procesados"),
          DTOutput("results_table"),
          hr(),
          h4("Tracks corregidas (preview)"),
          DTOutput("tracks_preview")
        )
      )
    )
  )
)

# ==============================================================================
server <- function(input, output, session) {

  # ---- Estado reactivo -------------------------------------------------------
  rv <- reactiveValues(
    tracks_ph   = NULL,   # tracks completo
    seg_all     = NULL,   # segmentos completo
    ids         = NULL,   # vector de IDs
    # individuo activo
    track_raw   = NULL,
    seg_raw     = NULL,
    track_prep  = NULL,   # track preparado (con col1)
    seg_prep    = NULL,   # seg preparado
    track_seg   = NULL,   # track unido con seg (para reclasif. seg)
    track_pt    = NULL,   # track para reclasif. puntos (parte de track_seg)
    # selecciones
    sel_seg     = vector(),
    sel_pts     = vector(),
    # acumulado
    segs_list   = list(),
    track_list  = list(),
    status_msg  = "Esperando datos..."
  )

  # ---- 3.1 Carga de datos ----------------------------------------------------
  observeEvent(input$load_data, {
    req(input$tracks_path, input$segs_path)
    tryCatch({
      rv$tracks_ph <- read_csv(input$tracks_path, show_col_types = FALSE)
      rv$seg_all   <- read_csv(input$segs_path,   show_col_types = FALSE)
      rv$ids       <- unique(rv$tracks_ph$ID)
      rv$status_msg <- paste0("✅ Datos cargados. ",
                               length(rv$ids), " individuos encontrados.")
    }, error = function(e) {
      rv$status_msg <- paste0("❌ Error al cargar: ", e$message)
    })
  })

  output$id_selector <- renderUI({
    req(rv$ids)
    selectInput("sel_id", "Seleccionar ID:", choices = rv$ids)
  })

  # ---- 3.2 Cargar individuo --------------------------------------------------
  observeEvent(input$load_ind, {
    req(rv$tracks_ph, rv$seg_all, input$sel_id)
    id <- input$sel_id

    seg_raw <- rv$seg_all %>% filter(ID == id)
    track_raw <- rv$tracks_ph %>% filter(ID == id)

    # Rellenar NAs de Breed2
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
    rv$status_msg <- paste0("👤 Individuo cargado: ", id)
  })

  # ---- 3.3 Vista general (gráficos ggplot) -----------------------------------
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
        ggtitle("Por estado (segmentación)"),
      p2 = ggplot(track_p) +
        geom_path(aes(x, y)) +
        geom_point(aes(x, y, fill = col1), size = 3, shape = 21) +
        scale_fill_identity(
          guide  = "legend", name = "Clasificación",
          breaks = unique(track_p$col1),
          labels = unique(track_p$lab)) +
        xlab("Longitude") + ylab("Latitude") + my_theme +
        ggtitle("Por estadio fenológico (Breed2)")
    )
  })

  output$plot_state  <- renderPlot({ make_plots()$p1 })
  output$plot_breed2 <- renderPlot({ make_plots()$p2 })

  # ---- 3.4 Mapa de segmentos -------------------------------------------------
  base_map_seg <- reactive({
    req(rv$track_seg)
    d <- rv$track_seg
    leaflet(d) %>%
      addProviderTiles(providers$CartoDB.PositronNoLabels) %>%
      addCircleMarkers(lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addPolylines(lng = ~x, lat = ~y, color = "black", weight = 2)
  })

  output$map_seg <- renderLeaflet({ base_map_seg() })

  # Selección de segmento al clic
  observeEvent(input$map_seg_marker_click, {
    req(rv$track_seg)
    cid <- as.numeric(input$map_seg_marker_click$id)
    clicked_seg <- rv$track_seg$Seg[cid]
    rv$sel_seg <- clicked_seg
    d <- rv$track_seg
    leafletProxy("map_seg") %>%
      clearMarkers() %>% clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addCircles(data = d[d$Seg == clicked_seg, ],
                 lat = ~y, lng = ~x, radius = 8,
                 color = "red", weight = 2, opacity = 1,
                 fill = FALSE, group = "sel")
  })

  # Actualizar segmento
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
    # Propagar al track_pt para que parta del resultado del paso anterior
    rv$track_pt <- new_d
  })

  # Limpiar selección segmento
  observeEvent(input$seg_clear, {
    rv$sel_seg <- vector()
    d <- rv$track_seg
    leafletProxy("map_seg") %>%
      clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d)))
  })

  # Reset segmento
  observeEvent(input$seg_reset, {
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

  # Tabla resumen de segmentos
  output$seg_table <- renderDT({
    req(rv$track_seg)
    rv$track_seg %>%
      dplyr::select(Seg, Ph.stage = Breed2, Start, End, Eq, Confidence, Check) %>%
      distinct() %>% drop_na(Ph.stage) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- 3.5 Mapa de puntos ----------------------------------------------------
  base_map_pt <- reactive({
    req(rv$track_pt)
    d <- rv$track_pt
    leaflet(d) %>%
      addProviderTiles(providers$CartoDB.PositronNoLabels) %>%
      addCircleMarkers(lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addPolylines(lng = ~x, lat = ~y, color = "black", weight = 2)
  })

  output$map_pt <- renderLeaflet({ base_map_pt() })

  # Selección de puntos al clic (acumulativa con Ctrl simulado: cada clic agrega)
  observeEvent(input$map_pt_marker_click, {
    req(rv$track_pt)
    cid <- as.numeric(input$map_pt_marker_click$id)
    # Toggle: si ya está, quitar; si no, añadir
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

  # Actualizar puntos
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

  # Limpiar selección puntos
  observeEvent(input$pt_clear, {
    rv$sel_pts <- vector()
    d <- rv$track_pt
    leafletProxy("map_pt") %>%
      clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d)))
  })

  # Reset puntos
  observeEvent(input$pt_reset, {
    rv$track_pt <- rv$track_seg   # vuelve al resultado del paso de segmentos
    rv$sel_pts  <- vector()
    d <- rv$track_pt
    leafletProxy("map_pt") %>%
      clearMarkers() %>% clearGroup("sel") %>%
      addCircleMarkers(data = d, lat = ~y, lng = ~x, color = ~col1,
                       radius = 4, layerId = ~seq_len(nrow(d))) %>%
      addPolylines(data = d, lng = ~x, lat = ~y,
                   color = "black", weight = 2)
  })

  # Tabla de puntos seleccionados / track completo
  output$pt_table <- renderDT({
    req(rv$track_pt)
    d <- rv$track_pt
    show <- if (length(rv$sel_pts) > 0) d[rv$sel_pts, ] else d
    datatable(show, options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- 3.6 Guardar individuo actual ------------------------------------------
  observeEvent(input$save_ind, {
    req(rv$track_pt, rv$track_seg, input$sel_id)
    id <- input$sel_id

    # track final: rellena NAs si los hubiera
    track_final <- rv$track_pt
    if (any(is.na(track_final$Breed2))) {
      track_final <- track_final %>%
        tidyr::fill(Breed2, .direction = "down")
    }

    # seg_final derivado de track_seg reclasificado
    seg_final <- rv$track_seg %>%
      dplyr::select(Breed2, State, Seg, Start, End, Eq, Confidence, Check) %>%
      distinct() %>% drop_na(Breed2)

    rv$track_list[[id]] <- track_final
    rv$segs_list[[id]]  <- seg_final
    rv$status_msg <- paste0("✅ Individuo '", id, "' guardado. ",
                             "Total guardados: ", length(rv$track_list))
  })

  # ---- 3.7 Exportar todos ----------------------------------------------------
  observeEvent(input$save_all, {
    req(length(rv$track_list) > 0)
    track_all <- plyr::ldply(rv$track_list, data.frame)
    seg_all2  <- plyr::ldply(rv$segs_list,  data.frame)

    dir.create("output", showWarnings = FALSE)
    write_csv(track_all, "output/Phenology_tracks_corrected.csv")
    write_csv(seg_all2,  "output/Segmentation_results_corrected.csv")

    rv$status_msg <- paste0("💾 Exportado: ",
                             nrow(track_all), " posiciones / ",
                             nrow(seg_all2),  " segmentos guardados.")
  })

  # ---- 3.8 Tab Resultados ----------------------------------------------------
  output$results_table <- renderDT({
    req(length(rv$track_list) > 0)
    data.frame(ID = names(rv$track_list),
               n_puntos = sapply(rv$track_list, nrow)) %>%
      datatable(options = list(pageLength = 20))
  })

  output$tracks_preview <- renderDT({
    req(length(rv$track_list) > 0)
    plyr::ldply(rv$track_list, data.frame) %>%
      head(200) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- Status ----------------------------------------------------------------
  output$status <- renderText({ rv$status_msg })
}

# == 4. Lanzar la app ==========================================================
shinyApp(ui = ui, server = server)
