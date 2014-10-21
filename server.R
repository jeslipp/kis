shinyServer(function(input, output, session) {
  dataInput <- reactive({
    get(input$dataset)
  })
  querySelection <- reactive({
    if (input$selection == "kinase") {
      unique(dataInput()$kinase)
    } else {
      unique(dataInput()$compound)
    }
  })
  observe({
    updateSelectInput(session, "query", choices = querySelection() )
  })
  observe({
    updateSelectInput(session, "exclusion", choices = querySelection() )
  })
  output$table = renderDataTable({
    if (input$selection == "kinase") {
      queryKinase(input, dataInput())
    } else {
      queryInhibitor(input, dataInput())
    }
  }, 
  options = list(pageLength = 20, lengthMenu = c(10, 20, 50))
  )
#   output$table <- renderTable({
#     if (input$selection == "kinase") {
#       queryKinase(input, dataInput())
#     } else {
#       queryInhibitor(input, dataInput())
#     }
#   })
})
