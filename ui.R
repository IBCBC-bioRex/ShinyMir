# ==============================
# UI
# ==============================

ui <- dashboardPage(
  skin = "blue",

  # ==============================
  # HEADER
  # ==============================

  dashboardHeader(title = "ShinyMIR"),
  
  # ==============================
  # SIDEBAR
  # ==============================
  dashboardSidebar(
    useShinyjs(),
    
    sidebarMenu(
      id = "sidebar",
      menuItem("miRNA - Disease",   tabName = "mirna_disease",      icon = icon("dna")),
      menuItem("miRNA - Gene",      tabName = "mirna_gene",         icon = icon("vial")),
      menuItem("Pathway Annotation",tabName = "over-representation",icon = icon("sitemap")),
      menuItem("miRNA - Reaction",  tabName = "mirna_reaction",     icon = icon("flask")),
      menuItem("Data upload",       tabName = "data_upload",        icon = icon("database"))
          ),
    hr(),
    div(style = "padding: 10px; padding-top: 0;",
      actionButton("download_report_btn", "Generate report", icon = icon("file-alt"),
                   style = "width:100%; background-color:#27ae60; color:white; border:none;")
    ),
    div(style = "padding: 10px; padding-top: 0;",
      actionButton("reset_all", "Reset all results", icon = icon("redo"),
                   style = "width:100%; background-color:#c0392b; color:white; border:none;")
    )
  ),
  
  # ==============================
  # BODY
  # ==============================
  dashboardBody(
    # Styles + scripts hoisted to <head> by Shiny (safe inside body)
    tags$head(
      tags$style(customStyles()),
      tags$script(scripts),
      tags$script(HTML('
Shiny.addCustomMessageHandler("trigger_zip_download", function(msg) {
  var bytes = atob(msg.b64);
  var arr = new Uint8Array(bytes.length);
  for (var i = 0; i < bytes.length; i++) arr[i] = bytes.charCodeAt(i);
  var blob = new Blob([arr], {type: "application/zip"});
  var url = URL.createObjectURL(blob);
  var a = document.createElement("a");
  a.href = url; a.download = msg.filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(function() { URL.revokeObjectURL(url); }, 3000);
});
      '))
    ),
    tabItems(

      # =========================================================
      # Data upload
      # =========================================================
      tabItem(
        tabName = "data_upload",
        h4(icon("database"), "Data upload"),
        p(icon("info-circle"), style = "color:#555;",
          "Upload custom tables to replace the built-in database. Original data can always be restored.",
          tags$br(),
          tags$em("Workflow: upload here and run Search in the analysis tabs to apply the new data.")),

        tabsetPanel(
          id = "data_upload_tabs",
          
          database_upload_tab_fixed_ids(
            title = "miRNA-disease source",
            source_id = "db_source",
            upload_id = "upload_database",
            update_id = "update_mda_database",
            restore_id = "restore_mda_database",
            apply_id = "apply_mda_database",
            help_id = "documentation_miRNA_disease_table_upload",
            preview_output_id = "mda_database_preview",
            label = "miRNA-disease",
            show_3p5p_button =TRUE
          ),
          
          database_upload_tab_fixed_ids(
            title = "metabolism source",
            source_id = "db_metabolism_source",
            upload_id = "upload_metabolism_database",
            update_id = "update_metabolism_database",
            restore_id = "restore_metabolism_database",
            apply_id = "apply_metabolism_database",
            help_id = "documentation_metabolism_table_upload",
            preview_output_id = "metabolism_database_preview",
            label = "metabolism",
            show_3p5p_button = FALSE
          )
        )
      )
      ,      

      # =========================================================
      # miRNA - Disease
      # =========================================================
      tabItem(
        tabName = "mirna_disease",
        
        h4(icon("dna"), "Analysis of miRNA-Disease Interactions"),
        p("Search curated miRNA-disease associations."),
        box(
          title = tagList(icon("filter"), "Filters"),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          
          fluidRow(
            column(
              6,
              box(
                title = tagList(icon("virus"), "miRNA selection"),
                width = 12,
                upload_list_ui("ul_mirna_select",
                tagList(
                  selectInput("mirna_select", "Choose one or more miRNAs:", choices = NULL, multiple = TRUE),
                  checkboxInput("filter_mirna_all_checkbox", "Use all miRNAs", value = TRUE),
                  actionButton("reset_mirna", "Reset miRNA selection", icon = icon("undo"))
                )
                )
              )
            ),
            column(
              6,
              box(
                title = tagList(icon("notes-medical"), "Disease selection"),
                width = 12,
                upload_list_ui("ul_disease_select",
                  tagList(
                    selectInput("do_category_filter", "Disease category (DO):",
                      choices = c("All categories" = ""), multiple = FALSE, selected = ""),
                    selectInput("disease_select", "Choose one or more diseases:", choices = NULL, multiple = TRUE),
                    checkboxInput("filter_disease_checkbox", "Use all available diseases", value = FALSE),
                    checkboxInput("filter_human_mda", "Human studies only", value = TRUE),
                    checkboxInput("filter_diseaseancestor_checkbox", "Include disease descendants", value = FALSE),
                    conditionalPanel(
                      condition = "input.filter_diseaseancestor_checkbox == true",
                      radioButtons("disease_descendant_depth", NULL,
                        choices  = c("All levels (recursive)" = "all",
                                     "First level only"       = "first"),
                        selected = "all", inline = TRUE)
                    ),
                    tags$small(style = "color:#888; display:block; margin-top:-8px; margin-bottom:6px;",
                      "Expands each disease to child terms in the DO hierarchy."),
                    actionButton("reset_disease", "Reset disease selection", icon = icon("undo"))
                  )
                )
              )
            )
          ),
          tags$details(
            tags$summary(style = "cursor:pointer; color:#3a8fc7; font-size:12px; margin-top:8px; margin-bottom:4px;",
              icon("sliders-h"), " Advanced filters"),
            fluidRow(style = "margin-top:8px;",
             column(3,
                    numericInput("min_assoc_mda", "Min. miRNA-disease publications:",
                                 value = 3L, min = 1L, max = 9999L, step = 1L)
             ),
              column(3,
                numericInput("min_degree_mirna_mda", "Min. miRNA degree (# diseases):",
                             value = 1L, min = 1L, max = 9999L, step = 1L)
              ),
              column(3,
                numericInput("min_degree_disease_mda", "Min. disease degree (# miRNAs):",
                             value = 1L, min = 1L, max = 9999L, step = 1L)
              ),
              column(3,
                checkboxGroupInput("sample_type_mda", "Sample type:",
                  choices = c("Circulating" = "circulating", "Exosome" = "exosome",
                              "Tissue" = "tissue", "Other" = "other"),
                  selected = character(0))
              )
            )
          )
        ),

        fluidRow(
          column(6,
            actionButton("obtain_mda", "Search associations", icon = icon("search"), class = "btn-primary")
          )
        ),
        hr(),
        tabsetPanel(
          id = "mirna_disease_tabs",
          tabPanel(
          title = tagList(icon("table"), "Results"),
          br(),
          fluidRow(
              column(6,
                     radioButtons("mda_view", NULL,
                                  choices = c("Table" = "table", "Network" = "net"),
                                  selected = "table", inline = TRUE)
              )
          ),          
          
          conditionalPanel(
            condition = "input.mda_view == 'table'",
            fluidRow(
              column(4,
                     selectInput("group_mda", "Group by:", choices = group_mdp_choices, selected = "miRNA and disease")
              )
            ),
            br(),
            htmlOutput("count_mda"),
            DTOutput("result_table_mda") %>% withSpinner(type = 6),
            br(),
            checkboxInput("download_abstract_mda", "Include article abstracts in download", value = FALSE),
            downloadButton("download_mda_csv", "Download CSV", class = "btn-sm btn-success")
          ),
          conditionalPanel(
            condition = "input.mda_view == 'net'",
            fluidRow(
              column(4, selectInput("mda_network_layout", "Layout:",
                          choices = .network_layout_choices,
                          selected = "forceAtlas2Based", width = "100%"))
            ),
            plotVisualizationUI("plot_md_net",
                                mode       = "network",
                                between_ui = uiOutput("mda_legend"))
          )
          ),
          # Co-targeting
          tabPanel(            title = tagList(icon("share-alt"), "Co-targeting"),
            radioButtons("mda_cotarget_node", "Node type:",
              choices = c("miRNA-miRNA" = "mirna_mirna", "Disease-Disease" = "disease_disease"),
              selected = "mirna_mirna", inline = TRUE),
            .cotarget_controls_ui("cotarget_mda"),
            plotlyOutput("cotarget_mda_heatmap", height = "700px") %>% withSpinner(type = 6),
            br(),
            h5(icon("table"), "Co-targeting pairs"),
            htmlOutput("cotarget_mda_count"),
            DTOutput("cotarget_mda_table") %>% withSpinner(type = 6)
          ),

          # Network Analysis: Robustness + Statistical Validation (merged)
          tabPanel(
            title = tagList(icon("chart-line"), "Network Analysis"),
            br(),
            p(icon("info-circle"), style = "color:#888; font-size:12px;",
              tags$em("Run Table search first. All analyses use the current filter selection.")),
            radioButtons("mda_network_analysis_mode", NULL,
              choices = c("Network robustness"      = "network",
                          "Co-targeting robustness" = "cotarget"),
              selected = "network", inline = TRUE),
            hr(),

            
            conditionalPanel(
              condition = "input.mda_network_analysis_mode == 'network'",
              fluidRow(
                column(3, numericInput("rob_mda_min_thresh", "Min threshold:",
                                       value = 1L, min = 1L, max = 100L, step = 1L)),
                column(3, numericInput("rob_mda_max_thresh", "Max threshold:",
                                       value = 20L, min = 1L, max = 100L, step = 1L))
              ),
              actionButton("run_rob_mda", "Run robustness", icon = icon("play"), class = "btn-primary btn-sm"),
              br(),
              plotlyOutput("rob_mda_plot", height = "420px") %>% withSpinner(type = 6),
              br(),
              h5(icon("table"), "Summary table"),
              DTOutput("rob_mda_table") %>% withSpinner(type = 6)
            ),

            
            conditionalPanel(
              condition = "input.mda_network_analysis_mode == 'cotarget'",
              fluidRow(
                column(3, numericInput("rob_mda_min_thresh_ct", "Min threshold:",
                                       value = 1L, min = 1L, max = 100L, step = 1L)),
                column(3, numericInput("rob_mda_max_thresh_ct", "Max threshold:",
                                       value = 20L, min = 1L, max = 100L, step = 1L)),
                column(3, sliderInput("rob_mda_jaccard_cutoff", "Jaccard cutoff:",
                                      min = 0, max = 1, value = 0.1, step = 0.05)),
                column(3, radioButtons("rob_mda_cotarget_node", "Node type:",
                  choices = c("miRNA-miRNA" = "mirna_mirna", "Disease-Disease" = "disease_disease"),
                  selected = "mirna_mirna", inline = TRUE))
              ),
              actionButton("run_rob_mda_cotarget", "Run co-targeting robustness",
                           icon = icon("play"), class = "btn-primary btn-sm"),
              br(),
              plotlyOutput("rob_mda_cotarget_plot", height = "420px") %>% withSpinner(type = 6),
              br(),
              h5(icon("table"), "Pair stability scores"),
              DTOutput("rob_mda_stability_table") %>% withSpinner(type = 6)
            )
          )
        )
      ),

      # =========================================================
      # miRNA - Gene
      # =========================================================
      tabItem(
        tabName = "mirna_gene",
        
        h4(icon("vial"), "Analysis of miRNA-Gene Interactions"),
        p("Search experimentally validated miRNA-gene target associations.",
          "miRNAs and diseases selected in the", tags$b("miRNA-Disease"), "tab propagate here automatically."),

        fluidRow(
          column(
            12,
            box(
              title = tagList(icon("filter"), "Filters"),
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              collapsible = TRUE,
              fluidRow(
                  column(
                    6,
                    box(title = tagList(icon("virus"), "miRNA selection"), width = 12,
                    radioButtons("mga_mirna_source", NULL,
                      choices = c( "From miRNA-Disease" = "disease","Search from scratch" = "scratch"),
                      selected = "disease", inline = TRUE),
                    conditionalPanel(
                      condition = "input.mga_mirna_source == 'scratch'",
                      checkboxInput("mga_use_3p5p_scratch", "Use 3p/5p information", value = FALSE),
                      upload_list_ui("ul_mirna_gene_scratch",
                        tagList(
                          selectizeInput("mirna_gene_select_scratch", "Select miRNAs:", choices = NULL, multiple = TRUE),
                          checkboxInput("mga_use_all_scratch", "Use all miRNAs", value = FALSE)
                        )
                      )
                    ),
                    conditionalPanel(
                      condition = "input.mga_mirna_source == 'disease'",
                      selectizeInput("mirna_gene_select", "miRNAs from miRNA\u2013Disease:", choices = NULL, multiple = TRUE),
                      checkboxInput("filter_mirna_common_checkbox", "Use common miRNAs between diseases", value = FALSE),
                      checkboxInput("use_all_mirnas_from_mda", "Use all miRNAs from Disease table", value = FALSE),
                      radioButtons("mga_mda_arm", "miRNA arm:",
                        choices  = c("Both (5p + 3p)" = "both","5p only" = "5p", "3p only" = "3p",
                                      "Most associated" = "max"),
                        selected = "both", inline = TRUE)
                    ),
                    actionButton("reset_mirna_gene", "Reset miRNA selection", icon = icon("undo"))
                    )
                  ),
                  column(
                    6,
                    box(
                      title = tagList(icon("dna"), "Gene selection"),
                      width = 12,
                      upload_list_ui("ul_gene_search",
                        tagList(
                          selectizeInput("gene_search_select", "Select one or more genes:", choices = NULL, multiple = TRUE),
                          checkboxInput("mga_use_all_genes",          "Use all genes",             value = FALSE),
                          checkboxInput("mga_only_metabolic_genes",   "Only metabolic genes",      value = FALSE),
                          tags$small(style = "color:#888; display:block; margin-top:-8px;",
                            "Genes with at least one reaction in Human-GEM.")
                        )
                      ),
                      actionButton("reset_gene_search", "Reset gene selection", icon = icon("undo"))
                    )
                  )
              ),
              checkboxInput("filter_human_mga", "Human studies only", value = TRUE),
              hr(),
              fluidRow(
                column(6, evidenceFiltersUI(suffix = "", reset_label = "Reset filters")),
                column(6, expressionFilterUI("_mga"))
              ),
              tags$details(
                tags$summary(style = "cursor:pointer; color:#3a8fc7; font-size:12px; margin-top:8px; margin-bottom:4px;",
                  icon("sliders-h"), " Advanced filters"),
                fluidRow(style = "margin-top:8px;",
                   column(3,
                          numericInput("min_assoc_mga", "Min. publications per miRNA-gene pair:",
                                       value = 2L, min = 1L, max = 100L, step = 1L)
                          
                   ),
                  column(3,
                    numericInput("min_degree_mirna_mga", "Min. miRNA degree (# genes):",
                                 value = 1L, min = 1L, max = 500L, step = 1L)
                  ),
                  column(3,
                    numericInput("min_degree_gene_mga", "Min. gene degree (# miRNAs):",
                                 value = 1L, min = 1L, max = 500L, step = 1L)
                  )
                )
              )
            )
          )
        ),

        fluidRow(
          column(6,
            actionButton("obtain_mga", "Search associations", icon = icon("search"), class = "btn-primary")
          )
        ),
        br(),
        tabsetPanel(
          # Co-targeting (without cotarget_min_assoc slider)
          tabPanel(
            title = tagList(icon("table"), "Results"),
            id = "mirna_gene_tabs",
            fluidRow(
              column(6,
                     radioButtons("mga_view", NULL,
                                  choices = c("Table" = "table", "Network" = "net"),
                                  selected = "table", inline = TRUE)
              )
            ),
            
            conditionalPanel(
              condition = "input.mga_view == 'table'",
              fluidRow(
                column(5, selectInput("group_mga", "Group by:", choices = group_mga_choices, selected = "miRNA and gene"))
              ),
              br(),
              htmlOutput("count_mga"),
              DTOutput("result_table_mga") %>% withSpinner(type = 6),
              br(),
              checkboxInput("download_abstract_mga", "Include article abstracts in download", value = FALSE),
              downloadButton("download_mga_csv", "Download CSV", class = "btn-sm btn-success")
            ),
            conditionalPanel(
              condition = "input.mga_view == 'net'",
              fluidRow(
                column(6, selectInput("mga_network_layout", "Layout:",
                            choices = .network_layout_choices,
                            selected = "forceAtlas2Based", width = "100%")),
                column(6, uiOutput("mga_pathway_filter_ui"))
              ),
              plotVisualizationUI("plot_mg_net",
                                  mode       = "network",
                                  title_x    = "miRNA", title_y = "Gene",
                                  between_ui = uiOutput("mga_legend"))
            )
            ),
          # Co-targeting (without cotarget_min_assoc slider)
          tabPanel(
            title = tagList(icon("share-alt"), "Co-targeting"),
            br(),
            radioButtons("mga_cotarget_node", "Node type:",
              choices = c("miRNA-miRNA" = "mirna_mirna", "Gene-Gene" = "gene_gene"),
              selected = "mirna_mirna", inline = TRUE),
            .cotarget_controls_ui("cotarget"),
            conditionalPanel(
              condition = "input.mga_cotarget_node === 'mirna_mirna' && !input.cotarget_cluster"
            ),
            plotlyOutput("cotarget_heatmap", height = "700px") %>% withSpinner(type = 6),
            br(),
            h5(icon("table"), "Co-targeting pairs"),
            htmlOutput("cotarget_count"),
            DTOutput("cotarget_table") %>% withSpinner(type = 6)
          ),

          # Network Analysis (robustness + validation merged)
          tabPanel(
            title = tagList(icon("chart-line"), "Network Analysis"),
            br(),
            p(icon("info-circle"), style = "color:#888; font-size:12px;",
              tags$em("Run Table search first. All analyses use the current filter selection.")),
            radioButtons("mga_network_analysis_mode", NULL,
              choices = c("Network robustness"      = "thresh",
                          "Co-targeting robustness" = "cotarget"),
              selected = "thresh", inline = TRUE),
            hr(),

            
            conditionalPanel(
              condition = "input.mga_network_analysis_mode == 'thresh'",
              fluidRow(
                column(3, numericInput("rob_mga_min_thresh", "Min threshold:",
                                       value = 1L, min = 1L, max = 100L, step = 1L)),
                column(3, numericInput("rob_mga_max_thresh", "Max threshold:",
                                       value = 20L, min = 1L, max = 100L, step = 1L))
              ),
              actionButton("run_rob_mga", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
              br(),
              plotlyOutput("rob_mga_plot", height = "420px") %>% withSpinner(type = 6),
              br(),
              h5(icon("table"), "Summary table"),
              DTOutput("rob_mga_table") %>% withSpinner(type = 6)
            ),
            conditionalPanel(
              condition = "input.mga_network_analysis_mode == 'cotarget'",
              fluidRow(
                column(3, numericInput("rob_mga_min_thresh_ct", "Min threshold:",
                                       value = 1L, min = 1L, max = 100L, step = 1L)),
                column(3, numericInput("rob_mga_max_thresh_ct", "Max threshold:",
                                       value = 20L, min = 1L, max = 100L, step = 1L)),
                column(3, sliderInput("rob_mga_jaccard_cutoff", "Jaccard cutoff:",
                                      min = 0, max = 1, value = 0.1, step = 0.05)),
                column(3, radioButtons("rob_mga_cotarget_node", "Node type:",
                  choices = c("miRNA-miRNA" = "mirna_mirna", "Gene-Gene" = "gene_gene"),
                  selected = "mirna_mirna", inline = TRUE))
              ),
              actionButton("run_rob_mga_cotarget", "Run", icon = icon("play"), class = "btn-primary btn-sm"),
              br(),
              plotlyOutput("rob_mga_cotarget_plot", height = "420px") %>% withSpinner(type = 6),
              br(),
              h5(icon("table"), "Pair stability scores"),
              DTOutput("rob_mga_stability_table") %>% withSpinner(type = 6)
            )
          )
        )
      ),
      
      
      # =========================================================
      # Pathway Annotation / Over-representation
      # =========================================================
      tabItem(
        tabName = "over-representation",

        h4(icon("sitemap"), "Pathway Annotation"),
        p("Tests which Reactome pathways are enriched in the gene set targeted by your miRNAs.",
          "Run a miRNA-Gene search first"),

        box(
          title = tagList(icon("filter"), "miRNA & Gene selection"),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          tags$small(style = "color:#888; display:block; margin-bottom:6px;",
            "Uses the miRNAs and genes from the last Search in the miRNA-Gene tab."),
          hr(),
          fluidRow(
            column(4,
              selectizeInput("ontology_mirna_select", "Select miRNAs:", choices = NULL, multiple = TRUE),
              checkboxInput("use_ontology_all_mirnas_from_mga", "Use all miRNAs", value = FALSE)
            ),
            column(4,
              selectizeInput("ontology_gene_select", "Select genes:", choices = NULL, multiple = TRUE),
              checkboxInput("filter_ontology_genes_reaction", "Use all genes", value = FALSE)
            )
          )
        ),

        tags$div(style = "display:none;",
          selectizeInput("ovr_pathway_upload_names", NULL, choices = NULL, multiple = TRUE)
        ),
        box(
          title = tagList(icon("filter"), "Pathway filters"),
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          width = 12,
          fluidRow(
            column(3,
              radioButtons("ovr_pathway_source", "Pathway database:",
                choices = c("Reactome (built-in)" = "reactome", "Custom gene sets (GMT)" = "gmt"),
                selected = "reactome", inline = FALSE)
            ),
            column(2, numericInput("ontology_top_n", "Top N pathways:", min = 5, max = 500, value = 20, step = 5))
          ),
          conditionalPanel(
            condition = "input.expression_filter_mode_mga != 'none' && input.expression_filter_mode_mga != null && input.ovr_pathway_source == 'reactome'",
            hr(),
            uiOutput("ovr_expr_context_banner"),
            fluidRow(
              column(4,
                sliderInput("ovr_expr_min_frac",
                  "Min % pathway genes expressed (universe filter):",
                  min = 0, max = 100, value = 20, step = 5, post = "%")
              ),
              column(4,
                sliderInput("ovr_tissue_coverage_min",
                  "Min tissue coverage in results:",
                  min = 0, max = 100, value = 0, step = 5, post = "%")
              ),
              column(4, br(), uiOutput("ovr_expr_filter_info"))
            )
          ),
          conditionalPanel(
            condition = "input.ovr_pathway_source == 'gmt'",
            hr(),
            fluidRow(
              column(6,
                h5(icon("upload"), "Upload GMT file (MSigDB, Hallmarks, KEGG, â€¦)"),
                tags$small(style = "color:#888; display:block; margin-bottom:6px;",
                  "GMT format: one gene set per line â€” name, URL/description, then tab-separated gene symbols."),
                upload_list_ui("ul_ovr_gmt",
                  tagList(
                    uiOutput("ovr_gmt_status")
                  )
                )
              )
            )
          ),
          conditionalPanel(
            condition = "input.ovr_pathway_source == 'reactome'",
            hr(),
            uiOutput("ovr_hierarchy_ui")
          )
        ),

        fluidRow(
          column(3,
            actionButton("run_overpresentation", "Run annotation", icon = icon("play"), class = "btn-primary")
          ),
          column(9, htmlOutput("count_mgoa"))
        ),
        br(),

        tabsetPanel(
          id = "ovr_results_tabs",

          tabPanel(
            title = tagList(icon("table"), "Table"),
            br(),
            DTOutput("pathway_overrepresentation_table") %>% withSpinner(type = 6),
            br(),
            downloadButton("download_overrep_csv", "Download CSV", class = "btn-sm btn-success")
          ),

          tabPanel(
            title = tagList(icon("bar-chart"), "Overview"),
            br(),
            conditionalPanel(
              condition = "input.ovr_pathway_source == 'reactome'",
              checkboxInput("ovr_barplot_group_l1", "Group bars by L1 macro-category (colour by L1)", value = FALSE)
            ),
            plotlyOutput("ovr_barplot_combined", height = "700px") %>% withSpinner(type = 6),
            conditionalPanel(
              condition = "input.ovr_pathway_source == 'reactome'",
              br(),
              h5("Summary by L1 pathway category", style = "color:#6b7a8d; font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:.7px; margin-bottom:6px;"),
              DTOutput("ovr_l1_summary_table") %>% withSpinner(type = 6)
            )
          ),

          tabPanel(
            title = tagList(icon("th"), "Heatmap"),
            br(),
            fluidRow(
              column(4,
                checkboxInput("ovr_heatmap_cluster_row",  "Cluster miRNAs (rows)",   value = FALSE),
                checkboxInput("ovr_heatmap_cluster_col",  "Cluster pathways (cols)", value = FALSE),
                checkboxInput("ovr_heatmap_transpose",    "Transpose",               value = FALSE)
              ),
              column(4,
                selectInput("ovr_heatmap_palette", "Color scale:",
                  choices = .heatmap_palette_choices,
                  selected = "YlGnBu")
              )
            ),
            hr(style = "margin-top:6px; margin-bottom:12px;"),
            plotlyOutput("mirna_pathway_heatmap", height = "800px") %>% withSpinner(type = 6)
          )
        )
      ),

      # =========================================================
      # miRNA - Reaction
      # =========================================================
      tabItem(
        tabName = "mirna_reaction",
        
        h4(icon("flask"), "Analysis of miRNA-Reaction Interactions"),
        p("Maps miRNA target genes to metabolic reactions via gene-protein rules (GPR) from Human-GEM.",
          tags$br(),
          tags$em("'By miRNA-Gene' starts from miRNAs and finds their metabolic reactions.",
                  " 'By reaction/subsystem' starts from the metabolic network and finds regulating miRNAs.")),

        box(
          title = tagList(icon("filter"), "Filters"),
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          fluidRow(
            column(4,
              h5(icon("dna"), "miRNA selection"),
              selectizeInput("mirna_reaction_select", "miRNAs from Gene table:", choices = NULL, multiple = TRUE),
              checkboxInput("filter_rmirna_common_checkbox", "Use common miRNAs", value = FALSE),
              checkboxInput("use_all_mirnas_from_mga", "Use all miRNAs", value = FALSE),
              actionButton("reset_mirna_reaction", "Reset", icon = icon("undo"), class = "btn-sm")
            ),
            column(4,
              h5(icon("dna"), "Gene selection"),
              selectizeInput("gene_reaction_select", "Genes from Gene table:", choices = NULL, multiple = TRUE),
              checkboxInput("filter_genes_reaction", "Use all genes from table", value = FALSE)
            )
          ),
            fluidRow(style = "margin-top:8px;",
              column(6,
                h6(icon("filter"), "Restrict to category/subsystem"),
                selectizeInput("mr_filter_l1", "Category:",
                  choices = category_choices, selected = NULL, multiple = TRUE,
                  options = list(placeholder = "All categories...")),
                selectizeInput("mr_filter_l2", "Subsystem:",
                  choices = subsystem_choices, selected = NULL, multiple = TRUE,
                  options = list(placeholder = "All subsystems...", maxOptions = 500L)),
                actionButton("reset_mr_pathway_filter", "Reset", icon = icon("undo"), class = "btn-sm")
              ),
              column(6,
                h6(icon("times-circle"), "Exclude subsystems/categories"),
                selectizeInput("mr_filter_exclude",
                  label = NULL,
                  choices = sort(union(category_choices, subsystem_choices)),
                  selected = NULL, multiple = TRUE,
                  options = list(placeholder = "Exclude...", maxOptions = 500L))
              ),
            )
        ),
        fluidRow(
          column(3,
            actionButton("obtain_mr", "Search associations", icon = icon("search"), class = "btn-primary")
          ),
          column(9, htmlOutput("count_mr"))
        ),
        br(),
        tabsetPanel(
          id = "mirna_reaction_tabs",

          tabPanel(
            title = tagList(icon("search"), "Results"),
            br(),
            fluidRow(
              column(3,
                radioButtons("mr_view", "View:",
                  choices = c("Table" = "table", "Network" = "network",
                              "Overview" = "barplot", "Heatmap" = "heatmap"),
                  selected = "table", inline = FALSE)
              ),
              conditionalPanel(
                condition = "input.mr_view == 'table'",
                column(4,
                  selectInput("group_mr", "Group by:", choices = group_mr_choices, selected = "None")
                )
              )
            ),

            # Table
            conditionalPanel(
              condition = "input.mr_view == 'table'",
              DTOutput("result_table_mr") %>% withSpinner(type = 6),
              br(),
              downloadButton("download_mr_csv", "Download CSV", class = "btn-sm btn-success")
            ),

            # Network
            conditionalPanel(
              condition = "input.mr_view == 'network'",
              fluidRow(
                column(6,
                  radioButtons("mr_net_color_by", "Colour reactions by:",
                    choices = c("Subsystem" = "subsystem", "Category" = "category"),
                    selected = "subsystem", inline = TRUE)
                ),
                column(6,
                  selectInput("mr_network_layout", "Layout:",
                    choices = .network_layout_choices,
                    selected = "forceAtlas2Based", width = "100%")
                )
              ),
              plotVisualizationUI("plot_mr_net",
                mode       = "network",
                title_x    = "miRNA", title_y = "reaction",
                between_ui = uiOutput("subsystem_legend_mr"))
            ),

            # Barplot (Overview)
            conditionalPanel(
              condition = "input.mr_view == 'barplot'",
              br(),
              fluidRow(
                column(4,
                  radioButtons("mr_barplot_group_by", "Group by:",
                    choices = c("Subsystem" = "subsystem", "Category" = "category"),
                    selected = "subsystem", inline = TRUE)
                )
              ),
              br(),
              plotlyOutput("mr_barplot", height = "700px") %>% withSpinner(type = 6)
            ),

            # Heatmap (metabolic miRNA - subsystem/reaction matrix)
            conditionalPanel(
              condition = "input.mr_view == 'heatmap'",
              br(),
              fluidRow(
                column(3,
                  radioButtons("metab_group_by", "Group by:",
                    choices = c("Category" = "category", "Subsystem" = "subsystem"),
                    selected = "subsystem")
                ),
                column(3),
                column(3, br(),
                  checkboxInput("metab_cluster", "Cluster rows/cols", value = FALSE)
                ),
                column(3,
                  selectInput("metab_heatmap_palette", "Color scale:",
                    choices = .heatmap_palette_choices,
                    selected = "YlGnBu")
                )
              ),
              br(),
              plotlyOutput("metab_heatmap", height = "800px") %>% withSpinner(type = 6),
              br(),
              tags$details(
                tags$summary(style = "cursor:pointer; font-weight:600; color:#555;",
                             icon("table"), " Summary table"),
                tags$div(style = "margin-top:8px;",
                  downloadButton("download_metab_heatmap_table", "Download CSV",
                                 class = "btn-sm btn-success", style = "margin-bottom:8px;"),
                  DTOutput("metab_heatmap_table") %>% withSpinner(type = 6)
                )
              )
            )
          ),

        tabPanel(
            title = tagList(icon("share-alt"), "Co-targeting"),
            radioButtons("mr_cotarget_node", "Node type:",
              choices = c("miRNA-miRNA" = "mirna_mirna", "Subsystem-Subsystem" = "subsystem_subsystem"),
              selected = "mirna_mirna", inline = TRUE),
            .cotarget_controls_ui("cotarget_mr"),
            plotlyOutput("cotarget_mr_heatmap", height = "700px") %>% withSpinner(type = 6),
            br(),
            h5(icon("table"), "Co-targeting pairs"),
            htmlOutput("cotarget_mr_count"),
            DTOutput("cotarget_mr_table") %>% withSpinner(type = 6)
          ),

          tabPanel(
            title = tagList(icon("chart-bar"), "Essentiality Analysis"),
            br(),
            p(icon("info-circle"),
              tags$em("Per-miRNA essentiality on the metabolic network.",
                      "For each miRNA: how many targeted reactions are fully controlled (ESS_FRAC >= threshold) vs partially.")),
            fluidRow(
              column(3,
                numericInput("mr_essentiality_threshold",
                             "Classify as essential if ESS_FRAC >=:",
                             value = 1.0, min = 0.0, max = 1.0, step = 0.05),
                tags$small(style = "color:#888;",
                  "Classifies reactions for this plot only.",
                  "To filter data upstream use the ESS_FRAC slider in the search panel.")
              ),
              column(3, br(),
                actionButton("run_mr_essentiality", "Compute essentiality",
                             icon = icon("play"), class = "btn-primary")
              )
            ),
            br(),
            plotlyOutput("mr_essentiality_plot", height = "500px") %>% withSpinner(type = 6),
            br(),
            h5(icon("table"), "Essentiality summary"),
            DTOutput("mr_essentiality_table") %>% withSpinner(type = 6),
            br(),
            downloadButton("download_mr_essentiality_csv", "Download CSV", class = "btn-sm btn-success")
          )
        )
      )
    )
    ) #end dashplot
  ) #end dashboardPage
