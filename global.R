# import libraries
library(shiny)
library(dplyr)
library(reshape2)

# load data
setwd("~/home/github/local/kis/")
load("anastassiadis2011.RData")
load("gao2013.RData")
load("kid2014.RData")

queryKinase <- function(input, data) {
  offtarget_rate <- data %>%
    group_by(compound) %>%
    summarise(percent_offtarget = round(sum(percent_activity <= input$query_cut, na.rm = TRUE) / n() * 100, 0))
    
  if (is.null(input$exclusion)) {
    query <- data %>%
      filter(kinase == input$query & percent_activity <= input$query_cut) %>%
      select(compound, kinase, percent_activity) %>%
      left_join(offtarget_rate, by = "compound") %>%
      arrange(percent_activity)
    
  } else {
    offtarget <- data %>%
      filter(kinase %in% input$exclusion) %>%
      dcast(compound ~ kinase, value.var = "percent_activity", fun.aggregate = mean)
    
    query <- data %>%
      filter(kinase == input$query & percent_activity <= input$query_cut) %>%
      left_join(offtarget_rate, by = "compound") %>%
      left_join(offtarget, by = "compound")
    
    query$offtarget <- apply(select(query, -c(compound, kinase, percent_activity, percent_offtarget)), 1, min, na.rm = TRUE)
    
    query <- query %>% 
      filter(offtarget > input$exclusion_cut) %>%
      select(-offtarget) %>%
      arrange(percent_activity) 
  }
  # return table
  if (nrow(query) == 0) {
    query <- data.frame("compound" = "N/A", 
                        "kinase" = "N/A", 
                        "percent_activity" = "N/A", 
                        "percent_offtarget" = "N/A")
  } else {
    query
  }
}

queryInhibitor <- function(input, data) {
  offtarget_rate <- data %>%
    group_by(kinase) %>%
    summarise(percent_offtarget = round(sum(percent_activity <= input$query_cut, na.rm = TRUE) / n() * 100, 0))
  
  if (is.null(input$exclusion)) {
    query <- data %>%
      filter(compound == input$query & percent_activity <= input$query_cut) %>%
      select(kinase, compound, percent_activity) %>%
      left_join(offtarget_rate, by = "kinase") %>%
      arrange(percent_activity)
  } else {
    offtarget <- data %>%
      filter(compound %in% input$exclusion) %>%
      dcast(kinase ~ compound, value.var = "percent_activity", fun.aggregate = mean)
    
    query <- data %>%
      filter(compound == input$query & percent_activity <= input$query_cut) %>%
      left_join(offtarget_rate, by = "kinase") %>%
      left_join(offtarget, by = "kinase")
    
    query$offtarget <- apply(select(query, -c(compound, kinase, percent_activity, percent_offtarget)), 1, min)
    
    query <- query %>% 
      filter(offtarget > input$exclusion_cut) %>%
      select(-offtarget) %>%
      arrange(percent_activity) 
  }
  # return table
  if (nrow(query) == 0) {
    query <- data.frame("compound" = "N/A", 
                        "kinase" = "N/A", 
                        "percent_activity" = "N/A", 
                        "percent_offtarget" = "N/A")
  } else {
    query
  }

}
