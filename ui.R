# Define UI for dataset viewer application
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Kinase Inhibitor Search"),
  
  # Sidebar with controls to select a dataset and specify the
  # number of observations to view
  sidebarLayout(
    sidebarPanel(
      # Select data set
      selectInput("dataset", "Choose a dataset:", 
                  choices = c("anastassiadis2011", "gao2013", "kid2014")),
      radioButtons("selection", label = "Choose query type: ",
                   choices = list("Kinase" = "kinase", "Inhibitor" = "inhibitor"), 
                   selected = "kinase"),
      # Input query 
      selectInput("query", label = "Choose query:", choices = ""),
      # Input query cut-off
      sliderInput("query_cut", label = "Residual activity:", min = 0, max = 100, value = 10), 
      # Input exclusion
      selectInput("exclusion", label = "Exclude:", choices = "", multiple = TRUE), 
      # Input query cut-off
      sliderInput("exclusion_cut", label = "Residual activity:", min = 0, max = 100, value = 50)
    ),

    # Show a summary of the dataset and an HTML table with the 
    # requested number of observations
    mainPanel(
      #tableOutput(outputId = "table")
      dataTableOutput('table')
    )
  )
))