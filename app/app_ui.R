# Web UI for shiny

#v0.7

#dashboard version

# HELPER FUNCTIONS ------------------------------

#adds id parameter to a menuItem. Allows calling the generated list item by selector reference
menuItemAddId = function( menuitem, id){
  #menuitem$attribs$id <- id;
  menuitem$children[[1]]$attribs$id <- id;
  return (menuitem);
}

# UI COMPONENTS -------------------------------
#header is the top stripe panel
header <- dashboardHeader(title = paste0("DreamTK v",app.version.major,".",app.version.minor, ".", app.version.revision));

#side bar is the expandable menu bar on the left side
sidebar <- dashboardSidebar(
  # Custom CSS to hide the default logout panel
  tags$head(tags$style(HTML('.shiny-server-account { display: none; }'))),
  
  # The dynamically-generated user panel
  uiOutput("userpanel"),
  
  div(id = "sidebardiv", 
  sidebarMenu(id = "sidebar",
    menuItemAddId( menuItem("Home", tabName = "hometab", icon = icon("home")), id="li_home" ),
    menuItemAddId( menuItem("Search", tabName = "searchtab", icon = icon("search")), id="li_search" ),
    menuItemAddId( menuItem("Analysis", tabName = "analysistab", icon = icon("pie-chart")), id="li_analysis" ),
    menuItemAddId( menuItem("MFA", tabName = "mfatab", icon = icon("line-chart")), id="li_mfa" ),
	menuItemAddId( menuItem("BER Analysis", tabName = "bertab", icon = icon("user", lib = "glyphicon")), id="li_ber" ),
    menuItemAddId( menuItem("Save/Load", tabName = "savetab", icon = icon("floppy-o")), id="li_save" ),
    menuItemAddId( menuItem("Help", icon = icon("question-circle"), startExpanded = FALSE,
                            actionLink(inputId = "link_help_overview", label = "Interface Help"),
                            actionLink(inputId = "link_help_search", label = "Search Help"),
                            actionLink(inputId = "link_help_save", label = "Data Backup Help"),
                            actionLink(inputId = "link_help_analysis", label = "Analysis Help")
                            ),
                   id = "li_help" ),
    menuItemAddId( menuItem("Quick Options", icon = icon("cog"), startExpanded = FALSE,
             radioButtons(inputId = "radio_listtype", label = "Chemical display options:", selected = "name",
                          choices = list("By CASN" = "casn",
                                         "By Name" = "name") ),
             #radioButtons(inputID = "radio_plotorientation", label = "Plot orientation options:", selected = "vertical",
              #            choices = list("Horizontal" ="horizontal",
              #                           "Vertical" = "vertical")),
             radioButtons(inputId = "radio_dbtype", label = "Database options:", selected = "rdata", width = "100%",
                          choices = list("MySQL: dreamtk_db (v0.7, based on httk-1.7, tcpl-v2)" = "mysql",
                                         "RData: DreamTKv0.8.RData (based on httk-1.7, tcpl-v2)" = "rdata") ),
             uiOutput(outputId="ui_database_status")
             ), 
             id="li_quickoptions" )
  ))
);

#body is the main dashboard are with all the required tabs
body <- dashboardBody(
  
  shinyjs::useShinyjs(),
  includeCSS( "./www/progressbar.css"),
  introjsUI(),
  
  tabItems(
    tabItem(tabName = "hometab",
            h2("Home"),
            fluidRow( div(id="welcomeboxdiv",
              box(status = "primary", title = "Welcome", collapsible = TRUE, width = 6,
                  h5("Welcome to DreamTK, an R application which facilitates toxicokinetic analysis of a variety of chemicals."),
                  h5("ToxCast Pipeline for High-Throughput Screening Data (tcpl v2.0) (Filer et al., 2016, US EPA) is the  primary chemical and assay database used by this application."),
                  h5("High-Throughput Toxicokinetics (httk v1.8) (Pearce et al., 2018) database is used to obtain the necessary toxicokinetic constants."),
				  h5("Only Assays that hit their cutoffs are considered for analysis.")
              ))
            ),
            fluidRow(
              box(status = "primary", title = "How To", collapsible = TRUE, width = 6,
                  h5("The following tutorials can help you familiarize with the application. These are also available from the Help tab on the Navigation menu."),
                  actionButton(inputId = "button_tour_overview", label = "Interface Tutorial"),
                  actionButton(inputId = "button_tour_search", label = "Search Tutorial"),
                  actionButton(inputId = "button_tour_save", label = "Data Backup Tutorial"),
                  actionButton(inputId = "button_tour_analysis", label = "Analysis Tutorial")
              )
            ),
            fluidRow(
              box(status = "primary", title = "News and Announcements", collapsible = TRUE, width = 6,
                  includeHTML("./www/news.html")
              )
            )
    ),
    
    tabItem(tabName = "searchtab",
            h3("Search Chemical Database"),
            fluidRow(
              
              tabBox(title = tagList(shiny::icon("compass")), id = "searchtabset",
                     
                     #search panel
                     tabPanel("Search", icon = icon("search"), fluid = TRUE,
                              fluidRow(
                                  column(12,
                                         div(
                                           style = "display: inline-block;vertical-align:baseline; width: 90%;",
                                           selectizeInput(inputId = "field_search", label = "", width = "100%",
                                                          options = list(placeholder = "Enter one or more chemicals separated by spaces",
                                                                         create = TRUE, createOnBlur = TRUE, createFilter = "^[\\w-,()\\[\\]]+$", persist = FALSE,
                                                                         maxOptions = 1000, loadThrottle = 800, openOnFocus = FALSE,
                                                                         delimiter = " ", hideSelected = TRUE, closeAfterSelect = TRUE,
                                                                         plugins = list("restore_on_backspace", "remove_button")),
                                                          multiple = TRUE, choices = NULL)
                                         ),
                                         div(
                                           style = "display: inline-block;vertical-align:70%; width: 9%;",
                                           actionButton(inputId = "button_search", label = "", icon = icon("search"), style="color:#0e76b7"),
                                           bsTooltip(id="button_search", "Search", placement = "bottom", trigger = "hover", options = NULL)
                                         ),
                                         radioButtons(inputId = "radio_searchtype", label = "Search options:", inline = TRUE, selected = "casn",
                                                      choices = list("By CASN" = "casn",
                                                                     "By Name" = "name") )
                                  )
                              )
                     ),
                     
                     #custom file panel
                     tabPanel("Custom File", icon = icon("table"), fluid = TRUE,
                              fluidRow(
                                  column(10, 
                                         fileInput(inputId = "file_customchem", label = "Choose CSV File", accept = c("csv", ".csv"), width = "100%"),
                                         uiOutput(outputId="ui_customchems_status")
                                  ),
                                  column(2,
                                         br(),
                                         actionButton(inputId = "button_load_customchems", label = "Parse file", icon = icon("folder-open-o"), 
                                                      style="padding-left:10px; padding-right:10px; padding-top:10px; padding-bottom:10px; white-space: normal;")
                                  )
                              ),
                              fluidRow(
                                column(12,
                                       br(),
                                       actionLink(inputId = "link_customchem_hint", label = "File requirements hint", icon = icon("hand-o-up")),
                                       hidden(fluidRow(id="panelCustomchemHint",
                                                       column(12,
                                                              h4("Template:"),
                                                              downloadLink(outputId = "button_customchem_template",label = tagList(shiny::icon("table"), "Download CSV Template"), 
                                                                           style="padding-left:10px; padding-right:10px; width:100%; white-space: normal; font-size:20px"),
                                                              HTML( 
                                                                str_c( "<p>Comma-separated value (.csv) file with first row containing column names. ",
                                                                       "<strong>'casn'</strong> column is mandatory. Other columns are optional, and if missing, will be looked up in the database.</p>",
                                                                       "<strong>Table columns:</strong>",
                                                                       "<blockquote style = 'font-size:14px; border-left: 10px solid #fff;'><strong>casn:</strong> CAS Registry Number",
                                                                       "<br><strong>name:</strong> Chemical name",
                                                                       "<br><strong>cytotoxicity_um:</strong> Cytotoxic concentration (<strong>Units:</strong> uM)",
                                                                       "<br><strong>cytotoxic:</strong> is the chemical considered cytotoxic at above concentration? Y/N flag",
                                                                       "<br><strong>mw:</strong> Molecular weight (<strong>Units:</strong> g/mol)",
                                                                       "<br><strong>human_funbound_plasma:</strong> Unbound fraction of chemical in blood (<strong>Units:</strong> none)",
                                                                       "<br><strong>human_clint:</strong> Intrinsic invitro hepatic clearance (<strong>Units:</strong> uL/min/10<sup>6</sup> cells)",
                                                                       "<br><strong>human_rblood2plasma:</strong> Blood to plasma concentration ratio (<strong>Units:</strong> none)",
                                                                       "<br><strong>log_kow:</strong> Octanol:Water partition coefficient at 25*C (<strong>Units:</strong> LOG10 value)",
                                                                       "<br><strong>pka:</strong> Equivalent chemical ionization constant (<strong>Units:</strong> - LOG10 value)</blockquote>"
                                                                )
                                                              )
                                                       )
                                       ))
                                )
                              )
                     )
                     
              )
              
            ),
            fluidRow(
              
              box(status = "primary", title = "Select and examine chemicals and chemical assays", collapsible = TRUE, width = 12, id = "chemlistbox",
                fluidRow(
					column(2,
					#if for some reason we change the colours primary won't work and we will have to figure a way to put colours in those switchtes. Style is not useable. https://rdrr.io/cran/shinyWidgets/src/R/input-pretty.R Gabriel
					prettySwitch(inputId = "select_hit", label = "Only Active Assays", value = TRUE, status = "primary",
									  fill = TRUE, bigger = TRUE, inline = TRUE,
									  width = NULL)
					),
					column(2,
					prettySwitch(inputId = "select_background", label = "Include Background Measurements", value = TRUE, status = "primary",
									  fill = TRUE, bigger = TRUE, inline = TRUE,
									  width = NULL)
					)
					
				),
				fluidRow(
                  column(12,
                         div(
                           style = "display: inline-block;vertical-align:bottom; width: 33%;",
                           selectInput(inputId = "select_chemical", label = "Searched chemical list", selectize = FALSE, size = 10,
                                       choices = NULL )
                         ),	
                         div(
                           style = "display: inline-block;vertical-align:bottom; width: 33%;",
                           selectInput(inputId = "select_assay", label = "Chemical assay list", selectize = FALSE, size = 10,
                                       choices = NULL )
                         ),
                         div(
                           style = "display: inline-block;vertical-align:bottom; width: 33%;",
                           selectInput(inputId = "select_assay_comp", label = "Assay component list", selectize = FALSE, size = 10,
                                       choices = NULL )
                         )
                  )

                ),
                fluidRow(
                  column(1,
                    actionButton(inputId = "button_delete_chem", label = "", icon = icon("trash-o"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_delete_chem", "Delete selected chemical", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(1,
                    actionButton(inputId = "button_delete_missing", label = "*", icon = icon("trash"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_delete_missing", "Delete missing chemicals", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(1,
                    actionButton(inputId = "button_list_missing", label = "*", icon = icon("list-ul"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_list_missing", "List missing chemicals", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(1,
                    actionButton(inputId = "button_clear_list", label = "", icon = icon("times"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;"),
                    bsTooltip(id="button_clear_list", "Clear chemical list", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(2,
                    downloadButton(outputId = "button_savechems", label = "Chem CSV", icon = icon("table"), 
                                   style="padding-left:10px; padding-right:10px; width:95%; white-space: normal;"),
                    bsTooltip(id="button_savechems", "Save information for all listed chemicals as a CSV file", placement = "bottom", trigger = "hover", options = NULL)
                  ),
                  column(2,
                    downloadButton(outputId = "button_saveassays", label = "Assays CSV", icon = icon("table"), style="padding-left:10px; padding-right:10px; width:95%; white-space: normal;"),
                    bsTooltip(id="button_saveassays", "Save assay information for the selected chemical as a CSV file", placement = "bottom", trigger = "hover", options = NULL)
                  )
                )							   
              )
              
          ),
          fluidRow(
            #chemical and assay information panels
            box(status = "primary", title = "Chemical info", collapsible = TRUE, width = 6,
                htmlOutput(outputId = "html_chemicalinfo")
            ),
            box(status = "primary", title = "Assay info", collapsible = TRUE, width = 6,
                htmlOutput(outputId = "html_assayinfo")
            )
          )
    ),
    
    tabItem(tabName = "analysistab",
            h3("Analyse chemicals"),
            #chemical selection panel
            fluidRow(
              box(status = "primary", title = "Select chemicals and desired analysis", collapsible = TRUE, width = 12,
                column(6,
                       wellPanel(id = "stats_selectcontrol",
                         selectizeInput(inputId = "select_chemical_stats", label = "Selected chemicals", 
                                        options = list(placeholder = "Click me to select chemicals",
                                                       maxOptions = 10000, loadThrottle = 800,
                                                       delimiter = " ", hideSelected = TRUE,
                                                       plugins = list("restore_on_backspace", "remove_button")),
                                        multiple = TRUE, choices = NULL),
                         fluidRow(
                           column(3,
                                  actionButton(inputId = "button_stats_selectall", label = "Select All", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;overflow:hidden")
                           ),
                           column(3,
                                  actionButton(inputId = "button_stats_deselectall", label = "Deselect All", icon = icon("ban"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;overflow:hidden")
                           )
                         )
                       ),
					   prettySwitch(inputId = "analyse_background", label = "Include Background Measurements", value = TRUE, status = "primary",
									  fill = TRUE, bigger = TRUE, inline = TRUE,
									  width = NULL)
                ),
                column(6,
                       wellPanel(id = "stats_optionscontrol",
                         checkboxGroupInput(inputId = "checkbox_stats", label = "Select statistics", 
                                            choices = c("Target Family Counts and ac50 values" = "tfcounts",
                                                        "Hierarchical cluster heatmap of Target Subfamily activities" = "tfhm",
                                                        "Hierarchical cluster heatmap of Assay Endpoint activities" = "assayhm",
                                                        "Ac50 vs ScalarTop" = "scalartop_ac50",
                                                        "OED vs ScalarTop" = "scalartop_oed",
														"Burst Assay vs Not Burst Assay" = "ac50_box",
                                                        "Chemical ToxPI Plots (Individual)" = "toxpi",
                                                        "Chemical ToxPI Plots (Cytotoxicity)" = "toxpi2",
                                                        "Chemical ToxPI Plots (Grouped)" = "toxpigroup"),
                                            selected = NULL),
						#those are magic numbers. Do not change them.					
						fluidRow(
                           column(4,
								actionButton(inputId = "button_select_all_stats", label = "Select/Deselect all", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; white-space: normal; width:100%; display:inline-block; overflow:hidden")
                                  
                           ),
                           column(3, offset = 0.5,
                                actionButton(inputId = "button_stats_run", label = "Run Stats", icon = icon("bar-chart"), style="padding-left:10px; padding-right:10px; white-space: normal; display:inline-block; overflow:hidden")  
                           )
                         ),
					    busyIndicator(text="Working...")
                       )
                )
              )
            ),
            fluidRow(
              box(status = "primary", title = "Analysis results", collapsible = TRUE, width = 12,
                uiOutput(outputId = "ui_stats")
              )
            )

    ),
    
    tabItem(tabName = "mfatab",
            h3("Multiple Factor Analysis"),
            #chemical selection panel
            fluidRow(
              box(status = "primary", title = "Select chemicals and desired analysis", collapsible = TRUE, width = 6,
                  wellPanel(id = "mfa_selectcontrol",
                    selectizeInput(inputId = "select_chemical_mfa", label = "Selected chemicals", 
                                   options = list(placeholder = "Click me to select chemicals",
                                                  maxOptions = 10000, loadThrottle = 800,
                                                  delimiter = " ", hideSelected = TRUE,
                                                  plugins = list("restore_on_backspace", "remove_button")),
                                   multiple = TRUE, choices = NULL),
                    fluidRow(
                      column(3,
                             actionButton(inputId = "button_mfa_selectall", label = "Select All", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_mfa_deselectall", label = "Deselect All", icon = icon("ban"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_mfa_run", label = "Run MFA", icon = icon("bar-chart"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(1,
                             busyIndicator(text="Working...")
                      )
                    )
                  )
              )
            ),
            fluidRow(
              box(status = "primary", title = "Analysis results", collapsible = TRUE, width = 12,
                  uiOutput(outputId = "ui_mfa")
              )
            )
    ),
	
	tabItem(tabName = "bertab",
            h3("Biological Exposure Ratio Analysis"),
            #chemical selection panel
            fluidRow(
              box(status = "primary", title = "Select chemicals and desired analysis", collapsible = TRUE, width = 6,
                  wellPanel(id = "ber_selectcontrol",
                    selectizeInput(inputId = "select_chemical_ber", label = "Selected chemicals", 
                                   options = list(placeholder = "Click me to select chemicals",
                                                  maxOptions = 10000, loadThrottle = 800,
                                                  delimiter = " ", hideSelected = TRUE,
                                                  plugins = list("restore_on_backspace", "remove_button")),
                                   multiple = TRUE, choices = NULL),
                    fluidRow(
                      column(3,
                             actionButton(inputId = "button_ber_selectall", label = "Select All", icon = icon("mouse-pointer"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_ber_deselectall", label = "Deselect All", icon = icon("ban"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(3,
                             actionButton(inputId = "button_ber_run", label = "Run BER analysis", icon = icon("bar-chart"), style="padding-left:10px; padding-right:10px; width:100%; white-space: normal;")
                      ),
                      column(1,
                             busyIndicator(text="Working...")
                      )
                    )
                  )
				  
              ),
			  box(status = "primary", title = "Extra Information and Assumptions", collapsible = TRUE, width = 6, height = 250,
					p("The calculations are based on the", a( href = "https://pubs.acs.org/doi/10.1021/es502513w", "SHEDS-HT", style = "color: blue;", target = "_blank", rel = "noopener noreferrer"), " exposure model."),
					h4("Assumptions"),
					tags$ul(tags$li("Physical activity index  =  1.75"),
						tags$li("Basal Alveolar Ventilation Rate  =  15.7 m",tags$sup("3"),"/day"),
						tags$li("Vapor pressure  =  0.876 Pa")
						),
					
					h4("Warnings"),
					tags$ul(tags$li("Ac50 values are used as a surrogate for OED values as the vast majority of chemicals do not have any values for OED.")))
            ),
            fluidRow(
              box(status = "primary", title = "Analysis results", collapsible = TRUE, width = 12,
                  uiOutput(outputId = "ui_ber")
              )
            )
    ),
    
    tabItem(tabName = "savetab",
            h3("Save/Load workspace"),
            box(status = "primary", title = "Save", width = 3,
                       radioButtons(inputId = "radio_savetype", label = "Save options:", selected = "cm",
                                    choices = list("R Chemical object list" = "cm") ),
                       downloadButton(outputId = "button_savefile", label = "Save", style="padding-left:10px; padding-right:10px; white-space: normal;")
            ),
            box(status = "primary", title = "Load", width = 4,
                       radioButtons(inputId = "radio_loadtype", label = "Load options:", selected = "cm",
                                    choices = list("R Chemical object list" = "cm") ),
                       fileInput(inputId = "file_load", label = "Choose File"),
                       actionButton(inputId = "button_loadfile", label = "Load", style="padding-left:10px; padding-right:10px; white-space: normal;"),
                       busyIndicator(text="Working..."),
                       uiOutput(outputId="ui_load_status")
            )
    )
	)
);




# main ui function
ui <- dashboardPage( header, sidebar, body, skin = "blue" );