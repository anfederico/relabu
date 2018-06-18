library(shiny)

#runApp("/Users/anthonyfederico/Work/johnson/sampleviz/relabu")

ui <- navbarPage(
  title = "Testing",
  theme = "simplex",
  source(file.path("ui", "tab1.R"),  local = TRUE)$value
)

server <- function(input, output, session) {
  source(file.path("server", "tab1.R"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)

