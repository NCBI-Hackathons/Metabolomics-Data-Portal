require(shiny)
require(DT)
require(shinydashboard)
require(ggplot2)
require(igraph)
require(grid)
require(MetabolomicsDataPortal)
require(Hmisc)

# Make sure you are in the current directory of your local installation of the MetDataPortal package
pkgdir <- find.package("MetabolomicsDataPortal")
miller_data = readRDS(file.path(pkgdir, "data/Miller2015_Heparin.rds"))
miller_data = miller_data[,-which(colnames(miller_data) %in% c("FillRate", "Mean.NormPop", "STD.NormPop"))]
print(colnames(miller_data))
print(dim(miller_data))
miller_data = apply(miller_data, c(1,2), as.numeric)

wangler_data = readRDS(file.path(pkgdir, "data/Wangler2017_EDTA.rds"))
wangler_data = wangler_data[,-which(colnames(wangler_data) %in% c("Average", "Standard.dev"))]
print(colnames(wangler_data))
print(dim(wangler_data))
wangler_data = apply(wangler_data, c(1,2), as.numeric)

cohorts = list()
cohorts[["zsd"]] = colnames(wangler_data)
cohorts[["arg"]] = colnames(miller_data)[grep("Argininemia", colnames(miller_data))]
cohorts[["asld"]] = colnames(miller_data)[grep("AL", colnames(miller_data))]
cohorts[["btd"]] = colnames(miller_data)[grep("biotinidase", colnames(miller_data))]
cohorts[["cbl"]] = colnames(miller_data)[grep("CBL", colnames(miller_data))]
cohorts[["cit"]] = colnames(miller_data)[grep("Citrullinemia", colnames(miller_data))]
cohorts[["cys"]] = colnames(miller_data)[grep("Cystinosis", colnames(miller_data))]
cohorts[["ga"]] = colnames(miller_data)[grep("GA", colnames(miller_data))]
cohorts[["gamt"]] = colnames(miller_data)[grep("GAMT", colnames(miller_data))]
cohorts[["hmg"]] = colnames(miller_data)[grep("HMG", colnames(miller_data))]
cohorts[["holo"]] = colnames(miller_data)[grep("holocarboxylase", colnames(miller_data))]
cohorts[["hcys"]] = colnames(miller_data)[grep("Homocystinuria", colnames(miller_data))]
cohorts[["ia"]] = colnames(miller_data)[grep("IA", colnames(miller_data))]
cohorts[["lpi"]] = colnames(miller_data)[grep("LPI", colnames(miller_data))]
cohorts[["mcad"]] = colnames(miller_data)[grep("MCAD", colnames(miller_data))]
cohorts[["mcc"]] = colnames(miller_data)[grep("MCC", colnames(miller_data))]
cohorts[["mma"]] = colnames(miller_data)[grep("MMA", colnames(miller_data))]
cohorts[["mngie"]] = colnames(miller_data)[grep("MNGIE", colnames(miller_data))]
cohorts[["msud"]] = colnames(miller_data)[grep("MSUD", colnames(miller_data))]
cohorts[["otc"]] = colnames(miller_data)[grep("OTC", colnames(miller_data))]
cohorts[["paa"]] = colnames(miller_data)[grep("PAA", colnames(miller_data))]
cohorts[["pku"]] = colnames(miller_data)[grep("PKU", colnames(miller_data))]
cohorts[["tmhle"]] = colnames(miller_data)[grep("TMLHE", colnames(miller_data))]
cohorts[["vlcad"]] = colnames(miller_data)[grep("VLCAD", colnames(miller_data))]
cohorts[["xcre"]] = colnames(miller_data)[grep("Xcreatine", colnames(miller_data))]
cohorts[["none"]] = colnames(miller_data)[grep("None", colnames(miller_data))]
print(names(cohorts))

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
                                                   selectInput(inputId = "diagClass", label = "Select diagnosis.", choices = names(cohorts), selected = "zsd"),
                                                   checkboxGroupInput(inputId = "ptIDs", label = "Select patients.", choices = "")),
                                       align="left", width=12)
                                   ),
                          fluidRow(box(title="Set-based Enrichment Analyses", status="info", solidHeader=TRUE, width=4, collapsible=TRUE#,
                                       #tableOutput("msea")
                                      ), 
                                    box(title="Topological-based Enrichment Analyses", status="info", solidHeader=TRUE, width=8, collapsible=TRUE,
                                     dataTableOutput("cepa")
                                    )
                                  ),
                          fluidRow(box(title="Pathway Map", status="primary", solidHeader = TRUE,
                                       splitLayout(cellWidths=c("33%", "33%", "33%"),
                                                   selectInput(inputId = "pathwayMapId", label = "Pathway Map", 
                                                               choices =  c("All", "Arginine-Metabolism", "Ascorbate-Metabolism", "Asp-Glu-Metabolism", 
                                                                            "BCAA-Metabolism", "Benzoate-Metabolism", "Beta-Oxidation", "Bile-Acid-Metabolism", 
                                                                            "Carnitine-Biosynthesis", "Cholesterol-Synthesis", "Creatine-Metabolism", 
                                                                            "DicarboxylicAcid-Metabolism", "Eicosanoids", "Endocannabinoid-Synthesis", 
                                                                            "FattyAcid-Metabolism", "Fibrinogen-Cleavage-Peptides", "GABA-Shunt", 
                                                                            "Galactose-Metabolism", "Glutathione-Metabolism", "Gly-Ser-Thr-Metabolism", 
                                                                            "Glycogen-Metabolism", "Glycolysis", "Glycosylation", "Hemoglobin-Porphyrin-Metabolism", 
                                                                            "Histidine-Metabolism", "Inositol-Metabolism", "Ketone-Bodies", 
                                                                            "Lysine-Catabolism", "Met-Cys-Metabolism", "Mevalonate-Metabolism", 
                                                                            "Nicotinate-Nicotinamide-Metabolism", "Pantothenate-Metabolism", 
                                                                            "Pentose-Phosphate-Metabolism", "Phe-Tyr-Metabolism", "Phospholipid-Metabolism", 
                                                                            "Polyamine-Metabolism", "Proline-Metabolism", "Protein-Degradation", 
                                                                            "Purine-Metabolism", "Pyridoxal-Metabolism", "Pyrimidine-Metabolism", 
                                                                            "Riboflavin-Metabolism", "Secondary-Bile-Acids", "Sorbitol-Glycerol-Metabolism", 
                                                                            "Sphingolipid-Metabolism", "Steroid-Hormone-Biosynthesis", "TCA-Cycle", 
                                                                            "Thyroid-Hormone-Synthesis", "Tryptophan-Metabolism"),
                                                               selected="Arginine-Metabolism"),
                                                   sliderInput(inputId = "scalingFactor", label="Node Scaling Factor", min=1, max=5, step=1, value=1),
                                                   plotOutput("colorbar")),
                                       imageOutput("pathwayMap"),
                                       align="left", width=12, collapsible=TRUE)
                                   ),
                          fluidRow(box(title = "Patient Report", status="info", solidHeader = TRUE,
                                       downloadButton("downloadPatientReport", "Download Patient Report"),
                                       splitLayout(cellWidths=c("60%", "40%"), dataTableOutput("patientReport"), dataTableOutput("missingMets")),
                                       align="left", width=12, collapsible=TRUE)
                                   )
                          ),
                  tabItem(tabName="networks", h2("Interpret Using Data-Derived Networks"))
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
        updateCheckboxGroupInput(session, "ptIDs", choices = cohorts[[input$diagClass]], selected = cohorts[[input$diagClass]][1], inline = TRUE)
        
        if (input$diagClass=="zsd") {
          .GlobalEnv$all_data = wangler_data
        } else {
          .GlobalEnv$all_data = miller_data
        }
        print(head(all_data))
        report = reactive(getPatientReport(input, all_data))
        output$patientReport = DT::renderDataTable({
          d = report()$patientReport
          DT::datatable(d, rownames=FALSE, options=list(scrollX=TRUE))
        })
        output$downloadPatientReport <- downloadHandler(
          filename = function() { paste(input$diagClass, "-", input$patientID, ".txt", sep="") },
          content = function(file) { write.table(report()$patientReport, file, sep="\t", col.names = TRUE, row.names = FALSE) }
        )
        #output$msea = renderTable(getMSEA(input))
        output$cepa = renderDataTable({
            d = report()$patientReport
            nms = d[,1]
            vec = d[,2]
            names(vec) = nms
            return (shiny.get.cepa(z_vec=vec, pathway.name=input$pathwayMapId, pmap.path=file.path(pkgdir, "extdata")))
          })
        #output$spia = renderTable(getSPIA(input))
        
        observeEvent(input$pathwayMapId, priority=0, {
          print(input$pathwayMapId)
          observeEvent(input$scalingFactor, priority=-1, {
            pmap = reactive(isolate(getPathwayMap(input, all_data)))
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
      #output$print_diagnosis = renderText(input$diagnosis)
      #output$sim = renderDataTable(res()$sim, rownames=FALSE)
      #output$pt_cmp = renderForceNetwork(res()$pt_ig)
    } else {
      print("Sanity check. Should never print this statement.")
    }
  })



}

shinyApp(ui, server)
