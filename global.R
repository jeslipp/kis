# import libraries
library(shiny)
library(dplyr)
library(reshape2)

# load data
load("anastassiadis2011.RData")
load("gao2013.RData")
load("kid2014.RData")

queryKinase <- function(input, data) {
  
  # calculate offtarget rate of compounds at given cut-off activity
  percent_offtarget <- data %>%
    group_by(compound) %>%
    summarise(percent_offtarget = as.integer(sum(percent_activity <= input$query_cut, na.rm = TRUE) / n() * 100))
  
  # determine activity of specified offtarget kinases at given exclusion cut-off activity
  if (is.null(input$exclusion)) {
    query <- data %>%
      select(-gini_kinase) %>%
      filter(kinase == input$query & percent_activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "compound") %>%
      arrange(percent_activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             "Activity (%)" = percent_activity, 
             "Gini (Inhibitor)" = gini_compound, 
             "Off-target (% Kinome)" = percent_offtarget)
  } else {
    offtarget <- data %>%
      filter(kinase %in% input$exclusion) %>%
      group_by(compound) %>%
      filter(min(percent_activity) >= input$exclusion_cut) %>%
      select(kinase, compound, percent_activity) %>%
      dcast(compound ~ kinase, value.var = "percent_activity")
    
    # determine compounds that inhibit query kinase at given cut-off activity
    query <- data %>%
      select(-gini_kinase) %>%
      filter(kinase == input$query & percent_activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "compound") %>%
      inner_join(offtarget, by = "compound") %>%
      arrange(percent_activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             "Activity (%)" = percent_activity, 
             "Gini (Inhibitor)" = gini_compound, 
             "Off-target (% Kinome)" = percent_offtarget)
  }
  
  query
}

queryInhibitor <- function(input, data) {
  # calculate offtarget rate of kinases at given cut-off activity
  percent_offtarget <- data %>%
    group_by(kinase) %>%
    summarise(percent_offtarget = as.integer(sum(percent_activity <= input$query_cut, na.rm = TRUE) / n() * 100))
  
  # determine activity of specified offtarget compounds at given exclusion cut-off activity
  if (is.null(input$exclusion)) {
    query <- data %>%
      select(-gini_compound) %>%
      filter(compound == input$query & percent_activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "kinase") %>%
      arrange(percent_activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             "Activity (%)" = percent_activity, 
             "Gini (Kinase)" = gini_kinase, 
             "Off-target (% Inhibitors)" = percent_offtarget)
  } else {
    offtarget <- data %>%
      filter(compound %in% input$exclusion) %>%
      group_by(kinase) %>%
      filter(min(percent_activity) >= input$exclusion_cut) %>%
      select(kinase, compound, percent_activity) %>%
      dcast(kinase ~ compound, value.var = "percent_activity")
    
    # determine compounds that inhibit query kinase at given cut-off activity
    query <- data %>%
      select(-gini_compound) %>%
      filter(compound == input$query & percent_activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "kinase") %>%
      inner_join(offtarget, by = "kinase") %>%
      arrange(percent_activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             "Activity (%)" = percent_activity, 
             "Gini (Kinase)" = gini_kinase, 
             "Off-target (% Inhibitors)" = percent_offtarget)
  }
  
  query
}

# queryKinase <- function(input, data) {
#   offtarget_rate <- data %>%
#     group_by(compound) %>%
#     summarise(percent_offtarget = round(sum(percent_activity <= input$query_cut, na.rm = TRUE) / n() * 100, 0))
#     
#   if (is.null(input$exclusion)) {
#     query <- data %>%
#       filter(kinase == input$query & percent_activity <= input$query_cut) %>%
#       select(compound, kinase, percent_activity) %>%
#       left_join(offtarget_rate, by = "compound") %>%
#       arrange(percent_activity)
#     
#   } else {
#     offtarget <- data %>%
#       filter(kinase %in% input$exclusion) %>%
#       dcast(compound ~ kinase, value.var = "percent_activity", fun.aggregate = mean)
#     
#     query <- data %>%
#       filter(kinase == input$query & percent_activity <= input$query_cut) %>%
#       left_join(offtarget_rate, by = "compound") %>%
#       left_join(offtarget, by = "compound")
#     
#     query$offtarget <- apply(select(query, -c(compound, kinase, percent_activity, percent_offtarget)), 1, min, na.rm = TRUE)
#     
#     query <- query %>% 
#       filter(offtarget > input$exclusion_cut) %>%
#       select(-offtarget) %>%
#       arrange(percent_activity) 
#   }
#   # return table
#   if (nrow(query) == 0) {
#     query <- data.frame("compound" = "N/A", 
#                         "kinase" = "N/A", 
#                         "percent_activity" = "N/A", 
#                         "percent_offtarget" = "N/A")
#   } else {
#     query
#   }
# }
# 
# queryInhibitor <- function(input, data) {
#   offtarget_rate <- data %>%
#     group_by(kinase) %>%
#     summarise(percent_offtarget = round(sum(percent_activity <= input$query_cut, na.rm = TRUE) / n() * 100, 0))
#   
#   if (is.null(input$exclusion)) {
#     query <- data %>%
#       filter(compound == input$query & percent_activity <= input$query_cut) %>%
#       select(kinase, compound, percent_activity) %>%
#       left_join(offtarget_rate, by = "kinase") %>%
#       arrange(percent_activity)
#   } else {
#     offtarget <- data %>%
#       filter(compound %in% input$exclusion) %>%
#       dcast(kinase ~ compound, value.var = "percent_activity", fun.aggregate = mean)
#     
#     query <- data %>%
#       filter(compound == input$query & percent_activity <= input$query_cut) %>%
#       left_join(offtarget_rate, by = "kinase") %>%
#       left_join(offtarget, by = "kinase")
#     
#     query$offtarget <- apply(select(query, -c(compound, kinase, percent_activity, percent_offtarget)), 1, min)
#     
#     query <- query %>% 
#       filter(offtarget > input$exclusion_cut) %>%
#       select(-offtarget) %>%
#       arrange(percent_activity) 
#   }
#   # return table
#   if (nrow(query) == 0) {
#     query <- data.frame("compound" = "N/A", 
#                         "kinase" = "N/A", 
#                         "percent_activity" = "N/A", 
#                         "percent_offtarget" = "N/A")
#   } else {
#     query
#   }
# 
# }
