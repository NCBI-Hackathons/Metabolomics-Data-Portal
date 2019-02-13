require(shiny)
require(DT)
require(shinydashboard)
require(ggplot2)
require(grid)
require(Hmisc)
require(plotly)
require(igraph)
require(networkD3)
require(R.utils)
require(MolPhenoMatch)

# Make sure you are in the current directory of your local installation of the MetDataPortal package
setwd("/Users/lillian.rosa/Downloads/MetDataPortal/R/")
source("metDataPortal_appFns.r")

kmx = 15

ui = dashboardPage(
  dashboardHeader(title = "MetDataPortal"),
  dashboardSidebar(sidebarMenu(id = "tab",
                               menuItem("Pathway Knowledge", tabName = "pathway", icon = icon("user-circle-o")),
                               menuItem("Disease-Specific Networks", tabName="networks", icon= icon("search")))),
  dashboardBody(height="100%",
                tabItems(
                  tabItem(tabName="pathway",
                          h2("Interpret Using Pathway Knowledgebases"),
                          fluidRow(box(title="Select Patient(s)", status="warning", solidHeader = TRUE,
                                       splitLayout(cellWidths=c("20%", "80%"),
                                                   selectInput(inputId = "diagClass", label = "Select diagnosis.", choices = c("ARG", "ASLD", "OTC", "ZSD"), selected = "ZSD"),
                                                   checkboxGroupInput(inputId = "ptIDs", label = "Select patients.", choices = "")),
                                       align="left", width=12)
                                   ),
                          fluidRow(box(title="Pathway Map", status="primary", solidHeader = TRUE,
                                       splitLayout(cellWidths=c("33%", "33%", "33%"),
                                                   selectInput(inputId = "pathwayMapId", label = "Pathway Map", choices = ""),
                                                   sliderInput(inputId = "scalingFactor", label="Node Scaling Factor", min=1, max=5, step=1, value=1),
                                                   plotOutput("colorbar")),
                                       imageOutput("pathwayMap"),
                                       align="left", width=12, collapsible=TRUE)
                                   ),
                          fluidRow(box(title = "Patient Report", status="info", solidHeader = TRUE,
                                       downloadButton("downloadPatientReport", "Download Patient Report"),
                                       splitLayout(cellWidths=c("60%", "40%"), dataTableOutput("patientReport"), dataTableOutput("missingMets")),
                                       align="left", width=12, collapsible=TRUE)
                                   ),
                          fluidRow(box(title="Set-based Enrichment Analyses", status="info", solidHeader=TRUE, width=4, collapsible=TRUE#,
                                       #tableOutput("setbased")
                                       ), 
                                   box(title="Topological-based Enrichment Analyses", status="info", solidHeader=TRUE, width=8, collapsible=TRUE#,
                                       #dataTableOutput("topological")
                                       )
                                   )
                          ),
                  tabItem(tabName="networks",
                          h2("Interpret Using Data-Derived Networks"),
                          fluidRow(tabBox(title="", id="diagOrExtract",
                                          tabPanel("Modular Feature Extraction",
                                                   splitLayout(cellWidths=c("25%", "25%", "25%", "25%"),
                                                               #selectInput(inputId = "diagnosis", label="Select Disease-Specific Knowledge Graph", choices = names(cohorts), selected = names(cohorts)[1]),
                                                               selectInput(inputId = "ptID", label="Select Patient 1", choices=NULL),
                                                               selectInput(inputId = "ptID2", label="Select Patient 2", choices=NULL),
                                                               selectInput(inputId = "kmax", label="Select k", choices = c(2:kmx))),
                                                   splitLayout(cellWidths=c("60%", "40%"), dataTableOutput("sim"), forceNetworkOutput("pt_cmp"))), width=12)
                                   )
                          )
                )
              )
)

server = function(input, output, session) {
  observe({
    print(sprintf("%s tab is selected.", input$tab))
  })

  observeEvent(input$tab, {
    if (input$tab == "pathway") {
      # Pathway Analysis code
      observeEvent(input$diagClass, priority=1, {
        print(input$diagClass)
        #updateCheckboxGroupInput(session, "ptIDs", choices = cohorts[[input$diagClass]], selected = cohorts[[input$diagClass]])

        report = reactive(getPatientReport(input, .GlobalEnv$all_raw_data, .GlobalEnv$all_norm_data, .GlobalEnv$all_data))
        output$patientReport = DT::renderDataTable({
          d = report()$patientReport
          DT::datatable(d, rownames=FALSE, options=list(scrollX=TRUE))
        })
        output$missingMets = DT::renderDataTable(report()$missingMets, rownames = FALSE)
        output$downloadPatientReport <- downloadHandler(
          filename = function() { paste(input$biofluid, "-", input$patientID, ".txt", sep="") },
          content = function(file) { write.table(report()$patientReport, file, sep="\t", col.names = TRUE, row.names = FALSE) }
        )
        #output$pathwayEnrichment = renderTable(getPathwayEnrichment2(input))

        observeEvent(input$pathwayMapId, priority=0, {
          print(input$pathwayMapId)
          updateSelectInput(session, "pathwayMapId", choices = c("All", "Arginine Metabolism", "Ascorbate Metabolism", "Asp-Glu Metabolism",
                                        "BCAA Metabolism", "Benzoate Metabolism", "Beta-Oxidation", "Bile-Acid Metabolism",
                                        "Carnitine Biosynthesis", "Cholesterol Synthesis", "Creatine Metabolism", "Dicoarboxylic Acid Metabolism",
                                        "Eicosanoids", "Endocannabinoid Synthesis", "Fatty Acid Metabolism", "Fibrinogen Cleavage Peptides",
                                        "GABA Shunt", "Galactose Metabolism", "Glutathione Metabolism", "Gly-Ser-Thr Metabaolism", "Glycogen Metabolism",
                                        "Glycolysis", "Glycosylation", "Hemoglobin-Porphyrin Metabolism", "Histidine Metabolism", "Inositol Metabolism",
                                        "Ketone Bodies", "Lysine Catabolism", "Met-Cys Metabolism", "Mevalonate Metabolism", "Nicotinate-Nicotinamide Metabolism",
                                        "Pantothenate Metabolism", "Pentose-Phosphate Metabolism", "Phe-Tyr Metabolism", "Phospholipid Metabolism", "Polyamine Metabolism",
                                        "Proline Metabolism", "Protein Degradation", "Purine Metabolism", "Pyridoxal Metabolism", "Pyrimidine Metabolism",
                                        "Riboflavin Metabolism", "Secondary-Bile-Acids", "Sorbitol-Glycerol Metabolism", "Sphingolipid-Metabolism",
                                        "Steroid-Hormone Biosynthesis", "TCA Cycle", "Thyroid Hormone Synthesis", "Tryptophan Metabolism"),
                            selected="Arginine Metabolism")
          observeEvent(input$scalingFactor, priority=-1, {
            pmap = reactive(isolate(getPathwayMap(input, .GlobalEnv$all_data)))
            output$pathwayMap = renderImage({pmap()$pmap})
            output$colorbar = renderPlot({
              grid.newpage()
              grid.layout(nrow = 1, ncol = 1, just = c("right", "top"))
              grid.draw(pmap()$colorbar)
            }, height = 20)
          })
        })
      })

    } else if (input$tab == "networks") {
      data = eval(parse(text=sprintf("graphs$%s$data", input$diagnosis)))
      #updateSelectInput(session, "ptID", label="Select a patient to diagnose", choices = colnames(data), selected=colnames(data)[1])
      #updateSelectInput(session, "ptID2", label="Select a patient to diagnose", choices = colnames(data), selected=colnames(data)[2])

      #output$print_diagnosis = renderText(input$diagnosis)
      #output$algSig_dscore = renderDataTable(comparePatientModPerts(input), rownames=FALSE)
      #output$mds = renderPlotly(getMDS(input))
      #res = reactive(extractModPerts(input))
      #output$sim = renderDataTable(res()$sim, rownames=FALSE)
      #output$pt_cmp = renderForceNetwork(res()$pt_ig)
    } else {
      print("No tab selected")
    }
  })



}

shinyApp(ui, server)
