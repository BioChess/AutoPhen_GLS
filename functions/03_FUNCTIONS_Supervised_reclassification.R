################################################################################

reclassSegs <- function(track, seg_df) {
  require(shiny)
  require(leaflet)
  require(dplyr)
  require(DT)
  
  track2 <- track%>%
    mutate(
      State = as.integer(pmap_chr(
        list(date = date),
        ~ {
          idx <- which(seg_df$Start <= ..1 & seg_df$End >= ..1)
          if (length(idx) > 0) {
            seg_df$State[idx[1]]
          } else {
            NA_character_
          }
        }
      )),
      Seg = as.integer(pmap_chr(
        list(date = date),
        ~ {
          idx <- which(seg_df$Start <= ..1 & seg_df$End >= ..1)
          if (length(idx) > 0) {
            seg_df$Seg[idx[1]]
          } else {
            NA_character_
          }
        }
      )))
  
  track3 <- track2%>%
    left_join(., seg_df, by = join_by(Breed2, State, Seg))
  
  
  # Mapeo de Breed2 a colores
  breed_color_map <- list(
    "B" = "orange",
    "mig1" = "black",
    "mig2" = "grey",
    "W" = "steelblue"
  )
  
  # UI
  ui <- fluidPage(
    titlePanel("Selection of segments for reclassification"),
    sidebarLayout(
      sidebarPanel(
        h4("Edit panel"),
        selectInput("point_breed2", "Phenological stages/ behaviours",
                    choices = names(breed_color_map)),
        actionButton("update_breed2", "Update"),
        actionButton("clear_selection", "Clean selection"),
        actionButton("reset_process", "Reset"),
        actionButton("save_changes", "Save and close"),
        hr(),
        h4("Segments summary"),
        DTOutput("table")  # Mostrar la tabla de segmentos
      ),
      mainPanel(
        leafletOutput("map", height = 600, width = 600)
      )
    )
  )
  
  # Server
  server <- function(input, output, session) {
    original_points <- track3  # Guardar los datos originales
    points <- reactiveVal(track3)
    segments <- reactiveVal(seg_df)  # Tabla de segmentos inicial
    selected_segs <- reactiveVal(vector())  # Almacena segmentos seleccionados
    
    # Renderizar el mapa inicial
    output$map <- renderLeaflet({
      leaflet(points()) %>%
        addProviderTiles(providers$CartoDB.PositronNoLabels) %>%
        addCircleMarkers(
          data = points(),
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(points()))
        ) %>%
        addPolylines(
          data = points(),
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        )
    })
    
    # Actualizar la tabla resumida cuando cambien los puntos
    observe({
      new_segments <- points() %>%
       dplyr::select(Seg, Ph.stage =Breed2, 
                     Start, End,
                     Eq, Confidence,  Check)%>%
        dplyr::distinct()%>%
        drop_na(Ph.stage)
      segments(new_segments)
    })
    
    # Seleccionar puntos por segmento
    observeEvent(input$map_marker_click, {
      clicked_id <- as.numeric(input$map_marker_click$id)
      req(clicked_id)
      
      current_points <- points()
      clicked_seg <- current_points$Seg[clicked_id]
      
      selected_segs(clicked_seg)
      
      leafletProxy("map") %>%
        clearMarkers() %>%
        clearGroup("selection") %>%
        addCircleMarkers(
          data = current_points,
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(current_points))
        ) %>%
        addCircles(
          data = current_points[current_points$Seg == clicked_seg, ],
          lat = ~y, lng = ~x,
          radius = 8,
          color = "red",
          weight = 2,
          opacity = 1,
          fill = FALSE,
          group = "selection"
        )
    })
    
    # Actualizar la categoría Breed2 de los puntos seleccionados
    observeEvent(input$update_breed2, {
      req(selected_segs())
      
      new_points <- points()
      selected_rows <- which(new_points$Seg %in% selected_segs())
      
      new_breed <- input$point_breed2
      new_color <- breed_color_map[[new_breed]]
      
      new_points[selected_rows, "Breed2"] <- new_breed
      new_points[selected_rows, "col1"] <- new_color
      
      points(new_points)
      
      leafletProxy("map") %>%
        clearMarkers() %>%
        clearGroup("selection") %>%
        addCircleMarkers(
          data = new_points,
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(new_points))
        ) %>%
        addPolylines(
          data = new_points,
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        )
    })
    
    # Limpiar la selección
    observeEvent(input$clear_selection, {
      selected_segs(vector())
      leafletProxy("map") %>%
        clearGroup("selection") %>%
        addCircleMarkers(
          data = points(),
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(points()))
        ) %>%
        addPolylines(
          data = points(),
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        )
    })
    
    # Reiniciar el proceso
    observeEvent(input$reset_process, {
      points(original_points)  # Reiniciar los puntos al estado original
      selected_segs(vector())
      leafletProxy("map") %>%
        clearMarkers() %>%
        clearGroup("selection") %>%
        addCircleMarkers(
          data = points(),
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(points()))
        ) %>%
        addPolylines(
          data = points(),
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        )
    })
    
    # Mostrar la tabla resumida
    output$table <- renderDT({
      datatable(segments(), options = list(pageLength = nrow(seg_df)))
    })
    
    # Guardar y cerrar
    observeEvent(input$save_changes, {
      stopApp(points())
    })
  }
  
  return(runApp(shinyApp(ui, server)))
}

################################################################################
reclassPoints <- function(track) {
  require(shiny)
  require(leaflet)
  require(dplyr)
  require(DT)
  

  # Breed2-to-color mapping
  breed_color_map <- list(
    "B" = "orange",
    "mig1" = "black",
    "mig2" = "grey",
    "W" = "steelblue"
  )

  # UI
  ui <- fluidPage(
    titlePanel("Selection of points for reclassification"),
    sidebarLayout(
      sidebarPanel(
        h4("Edit panel"),
        selectInput("point_breed2", "Phenological stages/ behaviours:", 
                    choices = names(breed_color_map)),  # Now selecting Breed2, not color
        actionButton("update_breed2", "Update"),  # Renamed button
        actionButton("clear_selection", "Clean Selection"),
        actionButton("reset_process", "Reset"),
        actionButton("save_changes", "Save and Close"),
        hr(),
        h4("Track data"),
        DTOutput("table")
      ),
      mainPanel(
        leafletOutput("map", height = 600,width = 600)
      )
    )
  )
  
  # Server
  server <- function(input, output, session) {
    original_points <- track
    points <- reactiveVal(track)
    selected_points <- reactiveVal(vector())  # Store multiple selected points
    
    # Render initial map
    output$map <- renderLeaflet({
      leaflet(points()) %>%
        addProviderTiles(providers$CartoDB.PositronNoLabels) %>%
        addCircleMarkers(
          data = points(),
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(points()))
        ) %>%
        addPolylines(
          data = points(),
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        )
    })
    
    # Select multiple points on click
    observeEvent(input$map_marker_click, {
      selected_id <- as.numeric(input$map_marker_click$id)
      current_selection <- selected_points()
      
      # Toggle selection (deselect if already selected)
      if (selected_id %in% current_selection) {
        selected_points(setdiff(current_selection, selected_id))
      } else {
        selected_points(c(current_selection, selected_id))
      }
      
      # Update map with selection highlights (WITHOUT REMOVING PATHS)
      leafletProxy("map") %>%
        clearGroup("selection") %>%  # Only remove selection circles, NOT the path
        addCircles(
          data = points()[selected_points(), ],
          lat = ~y, lng = ~x,
          radius = 6, 
          color = "red", 
          weight = 2, 
          opacity = 1,
          fillOpacity = 0,
          group = "selection"  # Add selection highlights to a group
        )
    })
    
    # Update Breed2 values of selected points
    observeEvent(input$update_breed2, {
      req(selected_points())
      original_points <- track  # Guardar los datos originales
      new_points <- points()
      selected_ids <- selected_points()
      
      # Get the selected Breed2 category and its corresponding color
      new_breed <- input$point_breed2
      new_color <- breed_color_map[[new_breed]]
      
      # Apply the changes
      new_points[selected_ids, "Breed2"] <- new_breed
      new_points[selected_ids, "col1"] <- new_color
      
      points(new_points)
      # Limpiar la selección después de actualizar
      selected_points(vector())  # Resetear los puntos seleccionados
      
      # Refresh map with updated colors and keep selection highlights
      leafletProxy("map") %>%
        clearMarkers() %>%
        addCircleMarkers(
          data = new_points,
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(new_points))
        ) %>%
        addPolylines(
          data = new_points,
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        ) %>%
        clearGroup("selection")  # Only clear selection circles, keep the path
   
    })
    
    # Reiniciar el proceso
    observeEvent(input$reset_process, {
      points(original_points)  # Reiniciar los puntos al estado original
      selected_points(vector())
      leafletProxy("map") %>%
        clearMarkers() %>%
        clearGroup("selection") %>%
        addCircleMarkers(
          data = points(),
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(points()))
        ) %>%
        addPolylines(
          data = points(),
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        )
    })
    # Limpiar la selección
    observeEvent(input$clear_selection, {
      selected_points(vector())
      leafletProxy("map") %>%
        clearGroup("selection") %>%
        addCircleMarkers(
          data = points(),
          lat = ~y, lng = ~x,
          color = ~col1,
          radius = 4,
          layerId = ~seq_len(nrow(points()))
        ) %>%
        addPolylines(
          data = points(),
          lng = ~x, lat = ~y,
          color = "black",
          weight = 2
        )
    })
    # Show updated data table
    # Mostrar la tabla filtrada por puntos seleccionados
    output$table <- renderDT({
      # Obtener los datos y los puntos seleccionados
      current_points <- points()
      selected_ids <- selected_points()
      
      # Filtrar si hay puntos seleccionados
      if (length(selected_ids) > 0) {
        datatable(current_points[selected_ids, ], options = list(pageLength = 10))
      } else {
        datatable(current_points, options = list(pageLength = 10))
      }
    })
    
    # Save and return modified dataframe
    observeEvent(input$save_changes, {
      stopApp(points())  # Stops the app and returns the modified dataframe
    })
  }
  
  return(runApp(shinyApp(ui, server)))
}
################################################################################
my_theme <- theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 20),
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) 

################################################################################
supervisedClass <- function(seg_df, track_df, byseg = T, x11 = T) {
  
  # Check if the input data frames are valid
  if (nrow(seg_df) == 0) {
    stop("Segments data frame is empty. Please provide valid data.")
  }
  
  if (nrow(track_df) == 0) {
    stop("Track data frame is empty. Please provide valid data.")
  }
 
  # Formatear los datos iniciales
  track_df2 <- track_df %>%
    filter(!is.na(x) & !is.na(y)) %>%
    mutate(
      col1 = case_when(
        Breed2 == 'B' ~ 'orange',
        Breed2 == 'mig1' ~ 'black',
        Breed2 == 'mig2' ~ 'grey',
        Breed2 == 'W' ~ 'steelblue'
      )
    ) %>%
    dplyr::select(x, y, date, Breed, Breed2, col1) %>%
    mutate(
      x = round(x, digits = 2),
      y = round(y, digits = 2),
      date = as.Date(date)
    )
  seg_df2 <- seg_df %>%
    mutate(
      f_day = as.Date(f_day),
      l_day = as.Date(l_day),
      segment = row_number()
    ) %>%
    dplyr::select(ID, State = state, Seg = segment, Breed, Start = f_day, End = l_day, Eq, Confidence,Check) %>%
    mutate(Breed2 = case_when(
      Breed == 'mig' & cumsum(replace_na(Breed == 'W', FALSE)) > 0 & 
        rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) == 0 ~ 'mig2',
      Breed == "mig" & rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) > 0 ~ "mig1", 
      TRUE ~ Breed
    )) %>%
    dplyr::select(-Breed)
  

  
  print(seg_df2)
  
  # Plot the track classification and segmentation overview
  track_p <- track_df %>%
    filter(!is.na(x) & !is.na(y)) %>%
    mutate(
      col1 = case_when(
        Breed2 == 'B' ~ 'orange',
        Breed2 == 'mig1' ~ 'black',
        Breed2 == 'mig2' ~ 'grey',
        Breed2 == 'W' ~ 'steelblue'
      ),
      state = as.factor(state),
      lab = case_when(Breed2 == 'B' ~ 'Breeding',
                      Breed2 == 'mig1' ~ 'Migration post',
                      Breed2 == 'mig2' ~ 'Migration pre',
                      Breed2 == 'W' ~ 'Wintering')
    )
  breaks <- unique(track_p$col1)
  labels <- unique(track_p$lab)
  p1 <- ggplot(data = track_p) +
    geom_path(aes(x, y)) +
    geom_point(aes(x, y, color = state), size = 3) +    
    xlab('Longitude') + ylab('Latitude') + my_theme
  
  p2 <- ggplot(data = track_p) +
    geom_path(aes(x, y)) +
    geom_point(aes(x, y, fill = col1), size = 3, shape = 21) + 
    scale_fill_identity(guide = 'legend',
                        name = 'Track classification',
                        breaks = breaks,
                        labels = labels) +
    xlab('Longitude') + ylab('Latitude') + my_theme
  
  if(x11 == T){x11(width = 10, height = 8)} 
  plot(ggarrange(p1, p2))
  
  # Prompt for specific track if there are multiple
  if(byseg == T){
    repeat {
      flush.console()  # Force console to update before showing prompt
      
      # cat("Do you want to reclass any segment? (y/n): ")  # Print the question
      answer <- suppressWarnings(readline("Do you want to reclass any segment? (y/n): "))  # Read input without tryCatch
      
      flush.console()  # Ensure output is printed immediately
      
      # Convert input to lowercase and trim spaces
      answer <- trimws(tolower(answer))
      
      # Check if input is valid
      if (answer == "y") {
        cat("Proceeding with reclassification...\n")
        utrack <- reclassSegs(track_df2, seg_df2)
        return(utrack)
        break  # Exit loop after reclassification
      } else if (answer == "n") {
        utrack <- track_df2%>%
          mutate(
            State = as.integer(pmap_chr(
              list(date = date),
              ~ {
                idx <- which(seg_df2$Start <= ..1 & seg_df2$End >= ..1)
                if (length(idx) > 0) {
                  seg_df2$State[idx[1]]
                } else {
                  NA_character_
                }
              }
            )),
            Seg = as.integer(pmap_chr(
              list(date = date),
              ~ {
                idx <- which(seg_df2$Start <= ..1 & seg_df2$End >= ..1)
                if (length(idx) > 0) {
                  seg_df2$Seg[idx[1]]
                } else {
                  NA_character_
                }
              }
            )))
        
         utrack <- utrack%>%
          left_join(., seg_df2, by = join_by(Breed2, State, Seg))
        
        cat("Supervised segment reclassification is not required.\n")
        return(utrack)
        break  # Exit loop if answer is 'n'
      } else {
        cat("Invalid input. Please enter 'y' or 'n'.\n")
      }
    }
  } else if (byseg == F){
    repeat {
      flush.console()  # Force console to update before showing prompt
      
      # cat("Do you want to reclass any segment? (y/n): ")  # Print the question
      answer <- suppressWarnings(readline("Do you want to reclass any concrete point? (y/n): "))  # Read input without tryCatch
      
      flush.console()  # Ensure output is printed immediately
      
      # Convert input to lowercase and trim spaces
      answer <- trimws(tolower(answer))
      
      # Check if input is valid
      if (answer == "y") {
        cat("Proceeding with reclassification...\n")
        utrack <- reclassPoints(track_df2)
        return(utrack)
        break  # Exit loop after reclassification
      } else if (answer == "n") {
        utrack <- track_df2
        cat("Supervised point reclassification is not required.\n")
        return(utrack)
        break  # Exit loop if answer is 'n'
      } else {
        cat("Invalid input. Please enter 'y' or 'n'.\n")
      }
    }
  }
}


