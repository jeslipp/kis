# import libraries
library(shiny)
library(dplyr)
library(reshape2)

# load data
load("anastassiadis2011.RData")
load("gao2013.RData")
load("kid2014.RData")

queryKinase <- function(input, data) {

  if (is.null(input$exclusion)) {
    query <- data %>%
      filter(kinase == input$query & percent_activity <= input$query_cut) %>%
      select(compound, kinase, percent_activity) %>%
      arrange(percent_activity)
  } else {
    offtarget <- data %>%
      filter(kinase %in% input$exclusion) %>%
      dcast(compound ~ kinase, value.var = "percent_activity", fun.aggregate = mean)
    
    query <- data %>%
      filter(kinase == input$query & percent_activity <= input$query_cut) %>%
      left_join(offtarget, by = "compound")
    
    query$offtarget <- apply(select(query, -c(compound, kinase, percent_activity)), 1, min)
    
    query <- query %>% 
      filter(offtarget > input$exclusion_cut) %>%
      select(-offtarget) %>%
      arrange(percent_activity) 
  }
  # return table
  if (nrow(query) == 0) {
    query <- data.frame("compound" = "N/A", "kinase" = "N/A", "percent_activity" = "N/A")
  } else {
    query
  }
}

queryInhibitor <- function(input, data) {

  if (is.null(input$exclusion)) {
    query <- data %>%
      filter(compound == input$query & percent_activity <= input$query_cut) %>%
      select(kinase, compound, percent_activity) %>%
      arrange(percent_activity)
  } else {
    offtarget <- data %>%
      filter(compound %in% input$exclusion) %>%
      dcast(kinase ~ compound, value.var = "percent_activity", fun.aggregate = mean)
    
    query <- data %>%
      filter(compound == input$query & percent_activity <= input$query_cut) %>%
      left_join(offtarget, by = "kinase")
    
    query$offtarget <- apply(select(query, -c(compound, kinase, percent_activity)), 1, min)
    
    query <- query %>% 
      filter(offtarget > input$exclusion_cut) %>%
      select(-offtarget) %>%
      arrange(percent_activity) 
  }
  # return table
  if (nrow(query) == 0) {
    query <- data.frame("kinase" = "N/A", "compound" = "N/A", "percent_activity" = "N/A")
  } else {
    query
  }

}