shinyServer(function(input, output, session) {
  
  ### switch query type between kinase and inhibitor
  querySelection <- reactive({
    if (input$selection == "kinase") {
      sort(unique(kid$kinase))
    } else {
      sort(unique(kid$compound))
    }
  })
  
  ### update query selection based on query type
  observe({
    updateSelectInput(session, "query", choices = querySelection() )
  })
  
  ### update exclusion selection based on query type
  observe({
    updateSelectInput(session, "exclusion", choices = querySelection() )
  })
  
  ### render output as DataTable
  output$table = renderDataTable({
    if (input$selection == "kinase") {
      queryKinase(input, kid)
    } else {
      queryInhibitor(input, kid)
    }
  }, 
  options = list(pageLength = 20, lengthMenu = c(10, 20, 50))
  )
})
