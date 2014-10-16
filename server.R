# import libraries
library(shiny)
library(dplyr)

# load data
load("data.RData")

# Define server logic required to summarize and view the selected
# dataset
shinyServer(function(input, output) {
   
  # Retrieve the kinase for query
  queryInput <- reactive({
    input$query
  })
  
  # Retrieve the kinase for exclusion
  exclusionInput <- reactive({
    input$exclusion
  })
  
  # Show the first "n" observations
  output$table <- renderTable({
    if (is.null(exclusionInput())) {
      query <- data %>%
        filter(kinase == queryInput() & percent_activity <= input$query_cut) %>%
        select(compound, kinase, percent_activity) %>%
        arrange(percent_activity)
    } else {
      offtarget <- data %>% 
        filter(kinase %in% exclusionInput()) %>%
        #filter(kinase %in% c("Aurora A", "Aurora B")) %>%
        dcast(compound ~ kinase, value.var = "percent_activity")
      
      query <- data %>%
        filter(kinase == queryInput() & percent_activity <= input$query_cut) %>%
        #filter(kinase == "Aurora C" & percent_activity <= 10) %>%
        left_join(offtarget, by = "compound")
      
      query$offtarget <- apply(select(query, -c(compound, kinase, percent_activity)), 1, min)
      
      query %>% filter(offtarget > input$exclusion_cut) %>% 
        select(-offtarget) %>%
        arrange(percent_activity) 
    }
  })
})