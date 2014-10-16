

library(shiny)

# Define UI for dataset viewer application
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Kinase Inhibitor Search"),
  
  # Sidebar with controls to select a dataset and specify the
  # number of observations to view
  sidebarLayout(
    sidebarPanel(
      # Input query
      selectInput("query", 
                  label = "Choose a kinase for inhibition:", 
                  choices = unique(data$kinase)),
      # Input query cut-off
      sliderInput("query_cut", 
                  label = "Residual activity:",
                  min = 0, max = 100, value = 10), 
      # Input exclusion
      selectInput("exclusion", 
                  label = "Exclude off-target kinases:", 
                  choices = unique(data$kinase), 
                  multiple = TRUE), 
      # Input query cut-off
      sliderInput("exclusion_cut", 
                  label = "Residual activity:",
                  min = 0, max = 100, value = 50)
    ),

    # Show a summary of the dataset and an HTML table with the 
    # requested number of observations
    mainPanel(
      tableOutput(outputId = "table")
    )
  )
))