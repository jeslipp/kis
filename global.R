# import libraries
library(shiny)
library(dplyr)
library(reshape2)

# load data
load("./data/kid.RData")

queryKinase <- function(input, data) {
  
  # calculate offtarget rate of compounds at given cut-off activity
  percent_offtarget <- data %>%
    group_by(compound) %>%
    summarise(percent_offtarget = as.integer(sum(activity <= input$query_cut, na.rm = TRUE) / n() * 100))
  
  # determine activity of specified offtarget kinases at given exclusion cut-off activity
  if (is.null(input$exclusion)) {
    query <- data %>%
      select(-gini_kinase) %>%
      filter(kinase == input$query & activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "compound") %>%
      arrange(activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             Source = source, 
             "Activity (%)" = activity, 
             "Gini (Inhibitor)" = gini_compound, 
             "Off-target (% Kinome)" = percent_offtarget)
  } else {
    offtarget <- data %>%
      filter(kinase %in% input$exclusion) %>%
      group_by(compound) %>%
      filter(min(activity) >= input$exclusion_cut) %>%
      select(kinase, compound, activity) %>%
      dcast(compound ~ kinase, value.var = "activity")
    
    # determine compounds that inhibit query kinase at given cut-off activity
    query <- data %>%
      select(-gini_kinase) %>%
      filter(kinase == input$query & activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "compound") %>%
      inner_join(offtarget, by = "compound") %>%
      arrange(activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             Source = source, 
             "Activity (%)" = activity, 
             "Gini (Inhibitor)" = gini_compound, 
             "Off-target (% Kinome)" = percent_offtarget)
  }
  
  query
}

queryInhibitor <- function(input, data) {
  # calculate offtarget rate of kinases at given cut-off activity
  percent_offtarget <- data %>%
    group_by(kinase) %>%
    summarise(percent_offtarget = as.integer(sum(activity <= input$query_cut, na.rm = TRUE) / n() * 100))
  
  # determine activity of specified offtarget compounds at given exclusion cut-off activity
  if (is.null(input$exclusion)) {
    query <- data %>%
      select(-gini_compound) %>%
      filter(compound == input$query & activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "kinase") %>%
      arrange(activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             Source = source, 
             "Activity (%)" = activity, 
             "Gini (Kinase)" = gini_kinase, 
             "Off-target (% Inhibitors)" = percent_offtarget)
  } else {
    offtarget <- data %>%
      filter(compound %in% input$exclusion) %>%
      group_by(kinase) %>%
      filter(min(activity) >= input$exclusion_cut) %>%
      select(kinase, compound, activity) %>%
      dcast(kinase ~ compound, value.var = "activity")
    
    # determine compounds that inhibit query kinase at given cut-off activity
    query <- data %>%
      select(-gini_compound) %>%
      filter(compound == input$query & activity <= input$query_cut) %>%
      inner_join(percent_offtarget, by = "kinase") %>%
      inner_join(offtarget, by = "kinase") %>%
      arrange(activity) %>%
      rename(Inhibitor = compound, 
             Kinase = kinase, 
             "Activity (%)" = activity, 
             "Gini (Kinase)" = gini_kinase, 
             "Off-target (% Inhibitors)" = percent_offtarget)
  }
  
  query
}