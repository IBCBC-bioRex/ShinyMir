# R/plot_module.R

# ---------------------------------------------------------------------------
# plotVisualizationUI
# mode = "both"    : original layout with view_type radio (Graph | Heatmap)
# mode = "network" : shows only network section (no view_type radio)
# mode = "heatmap" : shows only heatmap section (no view_type radio)
# ---------------------------------------------------------------------------
plotVisualizationUI <- function(id,
                                show_button = TRUE,
                                title_x = "miRNA", title_y = "target",
                                between_ui = NULL,
                                mode       = "both") {
  ns <- NS(id)

  # ── Network section ───────────────────────────────────────────────────────
  network_section <- tagList(
    if (!is.null(between_ui)) {
      fluidRow(
        column(9, visNetworkOutput(ns("network"), width = "100%", height = "80vh") %>%
                    withSpinner(type = 6)),
        column(3, between_ui)
      )
    } else {
      visNetworkOutput(ns("network"), width = "100%", height = "80vh") %>%
        withSpinner(type = 6)
    },
    br(),
    actionButton(ns("snapshot_png"),
                 label = tagList(icon("camera"), " Save PNG"),
                 class = "btn-sm btn-default",
                 title = "Saves exactly what you see at current layout")
  )

  # ── Heatmap section ───────────────────────────────────────────────────────
  heatmap_section <- tagList(
    checkboxInput(ns("cluster_rows"), paste0("Cluster ", title_x, " (rows)"),  value = FALSE),
    checkboxInput(ns("cluster_cols"), paste0("Cluster ", title_y, " (columns)"), value = FALSE),
    plotlyOutput(ns("heatmap"), height = "1000px") %>% withSpinner(type = 6)
  )

  # ── Compose based on mode ─────────────────────────────────────────────────
  if (mode == "network") {
    tagList(br(), network_section)
  } else if (mode == "heatmap") {
    tagList(br(), heatmap_section)
  } else {
    # "both" – view_type radio
    tagList(
      if (show_button)
        actionButton(ns("obtain_plot"), "Visualize associations",
                     icon = icon("eye"), class = "btn-primary"),
      br(),
      radioButtons(ns("view_type"), "Visualization type:",
                   choices = c("Graph" = "graph", "Heatmap" = "heatmap"),
                   selected = "graph", inline = TRUE),
      conditionalPanel(
        condition = paste0("input['", ns("view_type"), "'] == 'graph'"),
        network_section
      ),
      conditionalPanel(
        condition = paste0("input['", ns("view_type"), "'] == 'heatmap'"),
        heatmap_section
      )
    )
  }
}


# ---------------------------------------------------------------------------
# plotVisualizationServer
# ---------------------------------------------------------------------------
plotVisualizationServer <- function(id, get_data_fn, type_network,
                                    x1, x2,
                                    extra_graph_args   = reactive(list()),
                                    extra_heatmap_args = reactive(list()),
                                    max_nodes          = maximum_nodo_graph,
                                    mode               = "both",
                                    external_threshold = NULL,
                                    external_trigger   = NULL,
                                    physics_input      = NULL) {
  moduleServer(id, function(input, output, session) {

    # ── Snapshot PNG – composite canvas (network + legend) ──────────────────
    observeEvent(input$snapshot_png, {
      if (mode == "heatmap") return()
      net_id   <- session$ns("network")
      filename <- paste0(id, "_snapshot_", format(Sys.Date(), "%Y%m%d"), ".png")
      shinyjs::runjs(sprintf(
        "(function(){
           var container = document.getElementById('%s');
           if (!container) return;
           var netCanvas = container.querySelector('canvas');
           if (!netCanvas) { alert('Render the network first.'); return; }

           // Parse legend items from sibling column (if present)
           var items = [];
           var row = container.closest('.row');
           if (row) {
             var legendInner = row.querySelector('[id*=\"legend\"] > div');
             if (legendInner) {
               legendInner.querySelectorAll('div').forEach(function(div) {
                 var st = div.getAttribute('style') || '';
                 if (!st.includes('background:') && !st.includes('background-color:')) return;
                 var cm = st.match(/background(?:-color)?:\\s*([^;]+)/);
                 if (!cm) return;
                 var color  = cm[1].trim();
                 var circle = st.includes('border-radius:50%%') || st.includes('border-radius: 50%%');
                 var label  = (div.parentElement || div).textContent.trim().replace(/^\\s+/, '');
                 // label is sibling text — get parent text minus the color-box text
                 var parent = div.parentElement;
                 if (parent) {
                   var txt = '';
                   parent.childNodes.forEach(function(n){
                     if (n.nodeType === 3) txt += n.textContent;
                   });
                   label = txt.trim();
                 }
                 if (label) items.push({ color: color, label: label, circle: circle });
               });
             }
           }

           // Build composite canvas
           var PAD = 12, LEGEND_W = 220, ROW_H = 24, HEADER_H = 36;
           var hasLegend = items.length > 0;
           var legendH   = hasLegend ? HEADER_H + items.length * ROW_H + PAD : 0;
           var W = netCanvas.width + (hasLegend ? LEGEND_W : 0);
           var H = Math.max(netCanvas.height, legendH + PAD * 2);

           var out = document.createElement('canvas');
           out.width = W; out.height = H;
           var ctx = out.getContext('2d');

           ctx.fillStyle = '#ffffff';
           ctx.fillRect(0, 0, W, H);
           ctx.drawImage(netCanvas, 0, 0);

           if (hasLegend) {
             var lx = netCanvas.width + PAD;
             ctx.strokeStyle = '#dddddd';
             ctx.lineWidth = 1;
             ctx.strokeRect(lx, PAD, LEGEND_W - PAD * 2, legendH);

             ctx.font = 'bold 13px sans-serif';
             ctx.fillStyle = '#222222';
             ctx.fillText('Legend', lx + 8, PAD + 20);

             ctx.font = '12px sans-serif';
             items.forEach(function(item, i) {
               var y = PAD + HEADER_H + i * ROW_H + ROW_H / 2;
               ctx.fillStyle = item.color;
               if (item.circle) {
                 ctx.beginPath();
                 ctx.arc(lx + 14, y, 7, 0, Math.PI * 2);
                 ctx.fill();
               } else {
                 ctx.fillRect(lx + 7, y - 7, 14, 14);
               }
               ctx.fillStyle = '#222222';
               ctx.fillText(item.label, lx + 28, y + 4);
             });
           }

           var link = document.createElement('a');
           link.download = '%s';
           link.href = out.toDataURL('image/png');
           document.body.appendChild(link);
           link.click();
           document.body.removeChild(link);
         })();",
        net_id, filename
      ))
    })

    # ── Unified trigger ──────────────────────────────────────────────────────
    .plot_trigger <- reactive({
      (input$obtain_plot %||% 0L) +
        if (!is.null(external_trigger)) external_trigger() else 0L
    })

    # ── Data state ───────────────────────────────────────────────────────────
    data_state <- eventReactive(.plot_trigger(), {
      get_data_fn()
    })

    # ── Display state ────────────────────────────────────────────────────────
    display_state <- reactive({
      list(
        cluster_rows    = isTRUE(input$cluster_rows),
        cluster_cols    = isTRUE(input$cluster_cols),
        count_threshold = if (!is.null(external_threshold)) external_threshold() else 1L,
        extra_graph     = extra_graph_args(),
        extra_heatmap   = extra_heatmap_args()
      )
    })

    data_for_plot <- reactive(data_state())

    # ── Network output ───────────────────────────────────────────────────────
    output$network <- renderVisNetwork({
      df   <- data_state()
      req(nrow(df) > 0)
      n_nodes <- length(unique(c(df[[x1]], df[[x2]])))
      shiny::validate(shiny::need(
        n_nodes <= 500L,
        paste0("Network too large (", n_nodes, " nodes). Apply stricter filters to go below 500 nodes.")
      ))
      disp <- display_state()
      layout_type <- if (!is.null(physics_input)) physics_input() else "forceAtlas2Based"
      args <- c(
        list(df                  = df,
             maximum_nodo_graph  = max_nodes,
             physics_type        = layout_type,
             filter_edge         = disp$count_threshold,
             type_network        = type_network,
             node_click_input_id = session$ns("network_node_click")),
        disp$extra_graph
      )
      net <- do.call(mirna_target_graph_fun, args)
      # Auto-capture composite canvas (network + legend) after stabilization
      net_el  <- session$ns("network")
      png_inp <- session$ns("net_png")
      net %>% visNetwork::visEvents(
        stabilizationIterationsDone = sprintf(
          "(function() {
            // Delay to ensure canvas is fully painted after stabilization
            setTimeout(function() {
              var container = document.getElementById('%s');
              if (!container) return;
              var netCanvas = container.querySelector('canvas');
              if (!netCanvas || netCanvas.width === 0 || netCanvas.height === 0) return;
              // Build composite with legend (mirrors snapshot_png logic)
              var items = [];
              var row = container.closest('.row');
              if (row) {
                var legendInner = row.querySelector('[id*=\"legend\"] > div');
                if (legendInner) {
                  legendInner.querySelectorAll('div').forEach(function(div) {
                    var st = div.getAttribute('style') || '';
                    if (!st.includes('background:') && !st.includes('background-color:')) return;
                    var cm = st.match(/background(?:-color)?:\\s*([^;]+)/);
                    if (!cm) return;
                    var color = cm[1].trim();
                    var parent = div.parentElement;
                    var label = '';
                    if (parent) {
                      parent.childNodes.forEach(function(n){ if (n.nodeType===3) label += n.textContent; });
                      label = label.trim();
                    }
                    if (label) items.push({color: color, label: label});
                  });
                }
              }
              var PAD=12, LEGEND_W=220, ROW_H=24, HEADER_H=36;
              var hasLegend = items.length > 0;
              var legendH = hasLegend ? HEADER_H + items.length * ROW_H + PAD : 0;
              var W = netCanvas.width + (hasLegend ? LEGEND_W : 0);
              var H = Math.max(netCanvas.height, legendH + PAD*2);
              var out = document.createElement('canvas');
              out.width = W; out.height = H;
              var ctx = out.getContext('2d');
              ctx.fillStyle = '#ffffff'; ctx.fillRect(0,0,W,H);
              ctx.drawImage(netCanvas, 0, 0);
              if (hasLegend) {
                var lx = netCanvas.width + PAD;
                ctx.strokeStyle='#dddddd'; ctx.lineWidth=1;
                ctx.strokeRect(lx, PAD, LEGEND_W-PAD*2, legendH);
                ctx.font='bold 13px sans-serif'; ctx.fillStyle='#222222';
                ctx.fillText('Legend', lx+8, PAD+20);
                ctx.font='12px sans-serif';
                items.forEach(function(item,i){
                  var y = PAD+HEADER_H+i*ROW_H+ROW_H/2;
                  ctx.fillStyle = item.color;
                  ctx.beginPath(); ctx.arc(lx+14, y, 7, 0, Math.PI*2); ctx.fill();
                  ctx.fillStyle='#222222'; ctx.fillText(item.label, lx+28, y+4);
                });
              }
              Shiny.setInputValue('%s', out.toDataURL('image/png'), {priority:'event'});
            }, 400);
          })()",
          net_el, png_inp
        )
      )
    })

    # ── Heatmap ──────────────────────────────────────────────────────────────
    heatmap_data <- reactive({
      ds <- data_state()
      req(nrow(ds) > 0)
      prepare_heatmap_data(ds, display_state()$count_threshold, x1, x2)
    })

    output$heatmap <- renderPlotly({
      mat  <- heatmap_data()
      req(!is.null(mat))
      shiny::validate(shiny::need(
        nrow(mat) * ncol(mat) <= 5000L,
        paste0("Heatmap too large (", nrow(mat), "×", ncol(mat), "). Apply stricter filters to reduce the result set.")
      ))
      disp <- display_state()
      render_mirna_target_heatmap(
        mat,
        disp$cluster_rows,
        disp$cluster_cols,
        x1, x2,
        col_annotation = disp$extra_heatmap$col_annotation
      )
    })

    return(list(
      data        = data_for_plot,
      filter_edge = reactive(display_state()$count_threshold),
      node_click  = reactive(input$network_node_click),
      net_png     = reactive(input$net_png)
    ))
  })
}
