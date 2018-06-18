source(file.path("utils", "helpers.R"),  local = TRUE)

load(file.path("data", "pstat.rda"))

tables <- pstat.extraction(pstat)
#tables <- files.extraction()
taxlevs <- colnames(tables$TAX)
covariates <- colnames(tables$SAM)

tabPanel("Relative Abundance",
  tabsetPanel(
    tabPanel("Sample Relative Abundance",
      tags$br(),
      sidebarLayout(
        sidebarPanel(

          # Sort the samples by a condition
          conditionalPanel(
            condition = "input.group_samples == false",
            selectizeInput('sra_select_conditions', 'Color Samples by Condition', choices=covariates, multiple=TRUE)
          ),
          conditionalPanel(
            condition = "input.group_samples == true",
            selectizeInput('gra_select_conditions', 'Color Samples by Condition', choices=covariates)
          ),

          # Sample aggregation
          checkboxInput("group_samples", "Group Samples by Condition"),

          # Select taxon level
          selectInput("sra_taxlev", "Tax Level", choices=taxlevs, selected="family"),

          # Dynamically generate based on tax level
          uiOutput("sra_order_organisms"),

          # Sort the bars
          radioButtons("sra_sort_by", "Sort By",
          c("No Sorting" = "nosort", "Conditions" = "conditions", "Organisms" = "organisms"), selected="nosort"),

          # Legend toggle
          checkboxInput("sra_show_legend", "Show Legend", value=TRUE),
          width=3
        ),
        mainPanel(
          plotlyOutput("ra_plot", width="800px", height="600px"),
          width=9        
        )
      )
    ),
    tabPanel("Heatmap",
      tags$br(),
      sidebarLayout(
        sidebarPanel(
          width=3
        ),
        mainPanel(
          width=9        
        )
      )
    )
  )
)
