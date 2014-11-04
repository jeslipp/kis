library(stringr)
library(reshape2)
library(tidyr)
library(dplyr)
#library(ineq)

### Anastassiadis 2011 data
# read data
# all at 0.5uM
path <- "./data/anastassiadis2011.csv"
raw <- read.csv(path, 
                skip = 1, 
                stringsAsFactors = FALSE, 
                check.names = TRUE, 
                header = FALSE)

# clean data
data <- raw[3:302, 1:179]
inhibitor_names <- raw[1, 1:179]
inhibitor_names <- str_replace_all(unlist(inhibitor_names), "[ ]+", " ")
inhibitor_names <- str_replace_all(unlist(inhibitor_names), "\x9a", "o")
#inhibitors <- raw[2, 1:179]
colnames(data) <- paste(inhibitor_names, "0.5", sep = "@")
colnames(data)[1] <- "kinase"


# convert data to long format
data <- melt(data, 
             id.vars = "kinase", 
             variable.name = "compound", 
             value.name = "percent_activity")
data$compound <- as.character(data$compound)
data$percent_activity <- as.integer(data$percent_activity)
data$percent_activity <- ifelse(data$percent_activity < 0, 0, ifelse(data$percent_activity > 100, 100, data$percent_activity))
#data$percent_activity <- as.numeric(data$percent_activity)
data <- data %>%
  group_by(compound) %>%
  mutate(gini_compound = round(ineq(100 - percent_activity, type = "Gini"), 3)) %>%
  group_by(kinase) %>%
  mutate(gini_kinase = round(ineq(100 - percent_activity, type = "Gini"), 3))
  
data <- data[complete.cases(data), ]
data <- select(data, kinase, compound, percent_activity, gini_kinase, gini_compound)
#data <- select(data, kinase, compound, percent_activity)
anastassiadis2011 <- data.frame(data)
save(anastassiadis2011, file = "./anastassiadis2011.RData")


### Gao 2013 data

# read data
path <- "./data/gao2013.csv"
raw <- read.csv(path, 
                skip = 1, 
                stringsAsFactors = FALSE, 
                check.names = TRUE, 
                header = FALSE)

# clean data
data <- raw[2:235, ]
inhibitor_names <- unlist(raw[1, -1])
Encoding(inhibitor_names) <- "UTF-8"
inhibitor_names <- str_split(inhibitor_names, " @ ")
inhibitors <- sapply(inhibitor_names, "[", 1)
inhibitors <- str_replace(inhibitors, "\x9a", "o")
concentration <- sapply(inhibitor_names, "[", 2)
concentration <- str_extract(concentration, "^[10]{1,2}")
inhibitors <- paste(inhibitors, concentration, sep = "@")
inhibitors <- c("kinase", inhibitors)
colnames(data) <- inhibitors

data <- data %>%
  melt(id.vars = "kinase", variable.name = "compound", value.name = "percent_activity")

data$compound <- as.character(data$compound)
data$percent_activity <- as.integer(data$percent_activity)
data$percent_activity <- ifelse(data$percent_activity < 0, 0, ifelse(data$percent_activity > 100, 100, data$percent_activity))
# 
data <- data %>%
  group_by(compound) %>%
  mutate(gini_compound = round(ineq(100 - percent_activity, type = "Gini"), 3)) %>%
  group_by(kinase) %>%
  mutate(gini_kinase = round(ineq(100 - percent_activity, type = "Gini"), 3))

data <- select(data, kinase, compound, percent_activity, gini_kinase, gini_compound)
#data <- select(data, kinase, compound, percent_activity)
data <- data[complete.cases(data), ]
gao2013 <- data.frame(data)
save(gao2013, file = "./gao2013.RData")


### Kinase screen mrc data

path_screen <- "./data/kinase_inhibitor_results_2014-10-18T04-26-51.csv"
raw <- read.csv(path_screen, stringsAsFactors = FALSE)
data <- raw %>% 
  select(Kinase, CNumber, Screen.Conc, Inhibition) %>%
  rename(kinase = Kinase, 
         compound = CNumber, 
         concentration = Screen.Conc, 
         percent_activity = Inhibition) %>%
  arrange(kinase)
#names(data) <- c("kinase", "compound", "concentration", "percent_activity")

path_names <- "./kinase_inhibitor_list_2014-10-18T04-00-06.csv"
inhibitor_names <- read.csv(path_names, stringsAsFactors = FALSE)
inhibitors <- inhibitor_names$Inhibitor
names(inhibitors) <- inhibitor_names$CNumber

data$compound <- inhibitors[as.character(data$compound)]
#data$compound <- as.character(data$compound)
data$concentration <- str_replace(data$concentration, " ", "")
data$concentration <- as.numeric(data$concentration)

data <- data %>%
  group_by(kinase, compound, concentration) %>%
  summarise(percent_activity = mean(percent_activity)) %>%
  ungroup()

data$percent_activity <- as.integer(data$percent_activity)
data$percent_activity <- ifelse(data$percent_activity < 0, 0, ifelse(data$percent_activity > 100, 100, data$percent_activity))

data$concentration <- as.character(data$concentration)
data$concentration[is.na(data$concentration)] <- "unknown"
data <- data %>%
  unite(inhibitor, compound, concentration, sep = "@") %>%
  rename(compound = inhibitor) %>%
  arrange(kinase)

#data$compound <- as.character(data$compound)
#data$percent_activity <- as.integer(data$percent_activity)



data <- data %>%
  group_by(compound) %>%
  mutate(gini_compound = round(ineq(100 - percent_activity, type = "Gini"), 3)) %>%
  group_by(kinase) %>%
  mutate(gini_kinase = round(ineq(100 - percent_activity, type = "Gini"), 3))

data <- select(data, kinase, compound, percent_activity, gini_kinase, gini_compound)
#data <- select(data, kinase, compound, percent_activity)

#data <- data.frame(data)
# compound 299 is duplicated
# data <- data[!duplicated(data), ]
data <- data[complete.cases(data), ]

kid2014 <- data.frame(data)
save(kid2014, file = "./kid2014.RData")
