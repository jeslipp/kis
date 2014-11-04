rm(list = ls())

library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)
library(xlsx)
library(ineq)

setwd("~/home/github/kis")

################################################################
### create mapping file for human kinases
### takes a couple of minutes and requires connection to the internet
################################################################

if (!file.exists("./data/kinases_human.RData")) {
  source("./clean_pkinfam.R")
}
load("./data/kinases_human.RData")

################################################################
### Mapping of screening names to standard names
################################################################

map_kinome <- function(query, reference) {
  reference_set <- str_split(reference, " ")
  map <- array(NA, dim = length(screen))
  
  for (i in seq_along(screen)) {
    success = FALSE
    # string matching using grep
    if (success == FALSE) {
      candidates <- sapply(screen[[i]], grep, reference)
      candidates <- candidates[sapply(candidates, length) == 1]
      idx <- unique(unlist(candidates))
      if (length(idx) == 1) {
        map[i] <- kinases_human$kinase[idx]
        success = TRUE
      }
    }
    # string matching using set intersection
    if (success == FALSE) {
      candidates <- lapply(reference_set, intersect, screen[[i]])
      idx <- which(sapply(candidates, length) > 0)
      if (length(idx) == 1) {
        map[i] <- kinases_human$kinase[idx]
        success = TRUE
      }
    }
  }
  map
}

################################################################
### Anastassiadis 2011 
### Link to paper: http://www.ncbi.nlm.nih.gov/pubmed/22037377
################################################################

path <- "./raw/nbt.2017-S2.xls"
raw <- read.xlsx2(path, 
                  sheetIndex = 1, 
                  startRow = 2,
                  stringsAsFactors = FALSE,
                  check.names = TRUE,
                  header = FALSE)

ana11 <- data.matrix(raw[3:302, 2:179])

screen <- raw[3:302, 1]
screen <- str_replace(screen, "c-", "")
screen <- str_split(screen, "/| ")
screen <- lapply(screen, toupper)
reference <- kinases_human$gene
mapping <- map_kinome(screen, reference)
kinases <- data.frame(screen = raw[3:302, 1], mapping = mapping)
kinases$mapping[15] <- "AURKA"
kinases$mapping[16] <- "AURKB"
kinases$mapping[17] <- "AURKC"
kinases$mapping[30] <- "CAMK1"
kinases$mapping[31] <- "PNCK"
kinases$mapping[57] <- "CSNK1A1"
kinases$mapping[58] <- "CSNK1D"
kinases$mapping[59] <- "CSNK1E"
kinases$mapping[60] <- "CSNK1G1"
kinases$mapping[62] <- "CSNK1G3"
kinases$mapping[63] <- "CSNK2A1"
kinases$mapping[123] <- "GCG2"
kinases$mapping[134] <- "INSR"
kinases$mapping[168] <- "MYLK"
kinases$mapping[200] <- "MAPK12"
kinases$mapping[203] <- "PAK1"
kinases$mapping[207] <- "PAK7"
kinases$mapping[219] <- "PRKACA"
kinases$mapping[220] <- "PRKACB"
kinases$mapping[222] <- "PRKCB"
kinases$mapping[223] <- "PRKCB"
kinases$mapping[224] <- "PRKCD"
kinases$mapping[225] <- "PRKCE"
kinases$mapping[226] <- "PRKCH"
kinases$mapping[228] <- "PRKCI"
kinases$mapping[231] <- "PRKCQ"
kinases$mapping[232] <- "PRKCZ"
kinases$mapping[234] <- "PRKG1"
kinases$mapping[235] <- "PRKG1"
rownames(ana11) <- kinases$mapping

inhib <- raw[1, 2:179]
inhib_cas <- raw[2, 2:179]
inhib <- str_replace_all(unlist(inhib), "[ ]+", " ")
inhib <- paste(inhib, "0.5", sep = " @ ") # all inhibitors were used at 0.5uM
colnames(ana11) <- inhib

ana11 <- melt(ana11,
              id.vars = rownames(ana11),
              varnames = c("kinase", "compound"),
              value.name = "activity")
ana11$kinase <- as.character(ana11$kinase)
ana11$compound <- as.character(ana11$compound)
ana11["source"] <- "ana11"

################################################################
### Gao 2013
### Link to paper: http://www.ncbi.nlm.nih.gov/pubmed/23398362
################################################################

path <- "./raw/bj4510313add03.xls"
raw <- read.xlsx2(path, 
                  sheetIndex = 1, 
                  startRow = 1,
                  stringsAsFactors = FALSE,
                  check.names = TRUE,
                  header = FALSE)

gao13 <- data.matrix(raw[3:236, 2:256])

screen <- raw[3:236, 1]
screen <- str_replace_all(screen, "\\(h\\)|_|-", "")
screen <- str_split(screen, "/| ")
screen <- lapply(screen, toupper)
reference <- kinases_human$gene
mapping <- map_kinome(screen, reference)
kinases <- data.frame(screen = raw[3:236, 1], mapping = mapping)
kinases$mapping[5] <- "PRKAA1"
kinases$mapping[6] <- "PRKAA2"
kinases$mapping[10] <- "AURKA"
kinases$mapping[11] <- "AURKB"
kinases$mapping[12] <- "AURKC"
kinases$mapping[20] <- "CAMK1"
kinases$mapping[21] <- "CAMK2B"
kinases$mapping[22] <- "CAMK2G"
kinases$mapping[23] <- "CAMK2D"
kinases$mapping[25] <- "CAMK1D"
kinases$mapping[37] <- "CSNK1G1"
kinases$mapping[38] <- "CSNK1G2"
kinases$mapping[39] <- "CSNK1G3"
kinases$mapping[40] <- "CSNK1D"
kinases$mapping[41] <- "CSNK2A1"
kinases$mapping[42] <- "CSNK2A2"
kinases$mapping[43] <- "KIT"
kinases$mapping[46] <- "RAF1"
kinases$mapping[48] <- "SRC"
kinases$mapping[87] <- "GSK3A"
kinases$mapping[88] <- "GSK3B"
kinases$mapping[89] <- "GCG2"
kinases$mapping[96] <- "IGF1R"
kinases$mapping[97] <- "CHUK"
kinases$mapping[98] <- "IKBKB"
kinases$mapping[99] <- "INSR"
kinases$mapping[100] <- "INSR"
kinases$mapping[107] <- "MAPK8"
kinases$mapping[108] <- "MAPK9"
kinases$mapping[118] <- "MAPK1"
kinases$mapping[128] <- "MAP2K7"
kinases$mapping[129] <- "MLCK"
kinases$mapping[132] <- "CDC42BPA"
kinases$mapping[133] <- "CDC42BPB"
kinases$mapping[149] <- "RPS6KB1"
kinases$mapping[153] <- "PAK7"
kinases$mapping[155] <- "MARK2"
kinases$mapping[157] <- "PDGFRA"
kinases$mapping[158] <- "PDGFRB"
kinases$mapping[159] <- "PDPK1"
kinases$mapping[160] <- "NUAK2"
kinases$mapping[164] <- "PRKACA"
kinases$mapping[165] <- "AKT1"
kinases$mapping[166] <- "AKT2"
kinases$mapping[167] <- "AKT3"
kinases$mapping[168] <- "PRKCA"
kinases$mapping[169] <- "PRKCB"
kinases$mapping[170] <- "PRKCB"
kinases$mapping[171] <- "PRKCG"
kinases$mapping[172] <- "PRKCD"
kinases$mapping[173] <- "PRKCE"
kinases$mapping[174] <- "PRKCZ"
kinases$mapping[175] <- "PRKD3"
kinases$mapping[176] <- "PRKCQ"
kinases$mapping[177] <- "PRKCI"
kinases$mapping[178] <- "PRKD1"
kinases$mapping[180] <- "PRKG1"
kinases$mapping[181] <- "PRKG1"
kinases$mapping[191] <- "ROCK1"
kinases$mapping[192] <- "ROCK2"
kinases$mapping[214] <- "TAOK1"
kinases$mapping[215] <- "TAOK2"
kinases$mapping[216] <- "TAOK3"
rownames(gao13) <- kinases$mapping

inhib <- unlist(raw[2, 2:256])
inhib <- str_replace(inhib, " µM", "")
colnames(gao13) <- inhib

gao13 <- melt(gao13, 
              id.vars = rownames(gao13),
              varnames = c("kinase", "compound"),
              value.name = "activity")
gao13$kinase <- as.character(gao13$kinase)
gao13$compound <- as.character(gao13$compound)
gao13["source"] <- "gao13"


################################################################
### Kinase inhibitor database
### Link to resource: http://www.kinase-screen.mrc.ac.uk/kinase-inhibitors
################################################################

path_screen <- "./raw/kinase_inhibitor_results_2014-10-18T04-26-51.csv"
raw <- read.csv(path_screen, stringsAsFactors = FALSE)
kid14 <- raw %>% 
  select(Kinase, CNumber, Screen.Conc, Inhibition) %>%
  rename(kinase = Kinase, 
         compound = CNumber, 
         concentration = Screen.Conc, 
         activity = Inhibition)

path_names <- "./raw/kinase_inhibitor_list_2014-10-18T04-00-06.csv"
inhibitor_names <- read.csv(path_names, stringsAsFactors = FALSE)
inhibitors <- inhibitor_names$Inhibitor
names(inhibitors) <- inhibitor_names$CNumber
kid14$compound <- inhibitors[as.character(kid14$compound)]
kid14$concentration <- str_replace(kid14$concentration, " ", "")
kid14$concentration <- as.numeric(kid14$concentration)
kid14$concentration[is.na(kid14$concentration)] <- "unknown"
kid14 <- kid14 %>%
  unite(inhibitor, compound, concentration, sep = " @ ") %>%
  rename(compound = inhibitor) %>%
  arrange(kinase)

screen <- unique(kid14$kinase)
screen <- str_replace_all(screen, "-| (hum)", "")
screen <- str_split(screen, " ")
screen <- lapply(screen, toupper)
reference <- kinases_human$gene
mapping <- map_kinome(screen, reference)
kinases <- data.frame(screen = unique(kid14$kinase), mapping = mapping)
kinases$mapping[5] <- "AURKA"
kinases$mapping[6] <- "AURKB"
kinases$mapping[7] <- "AURKC"
kinases$mapping[13] <- "CAMKK1"
kinases$mapping[14] <- "CAMKK2"
kinases$mapping[15] <- "CDK2"
kinases$mapping[16] <- "CDK9"
kinases$mapping[19] <- "CSNK1D"
kinases$mapping[20] <- "CSNK1G2"
kinases$mapping[21] <- "CSNK2A1"
kinases$mapping[43] <- "GSK3B"
kinases$mapping[49] <- "IKBKB"
kinases$mapping[50] <- "IKBKE"
kinases$mapping[51] <- "INSR"
kinases$mapping[91] <- "PAK7"
kinases$mapping[94] <- "PDPK1"
kinases$mapping[95] <- "PHKG1"
kinases$mapping[100] <- "PRKACA"
kinases$mapping[101] <- "AKT1"
kinases$mapping[102] <- "AKT2"
kinases$mapping[103] <- "PRKCA"
kinases$mapping[104] <- "PRKCG"
kinases$mapping[105] <- "PRKCZ"
kinases$mapping[112] <- "ROCK2"
kinases$mapping[115] <- "RPS6KC1"
kinases$mapping[122] <- "MLCK"
kinases$mapping[125] <- "TAOK1"
kinases$mapping[142] <- "MAPK14"
kinases$mapping[143] <- "MAPK11"
kinases$mapping[144] <- "MAPK13"
kinases$mapping[145] <- "MAPK12"
rownames(kinases) <- kinases$screen
kid14$kinase <- as.character(kinases[kid14$kinase, 2])
kid14["source"] <- "kid14"

### Unite data sets

kid <- rbind_list(ana11, gao13, kid14)
kid <- kid %>%
  dcast(kinase + source ~ compound, fun.aggregate = mean, value.var = "activity") %>%
  melt(id.vars = c("kinase", "source"), value.name = "activity", variable.name = "compound") %>%
  filter(!(is.nan(activity) | is.na(activity))) %>%
  mutate(activity = as.integer(ifelse(activity < 0, 0, ifelse(activity > 100, 100, activity)))) %>%
  group_by(compound) %>%
  mutate(gini_compound = round(ineq(100 - activity, type = "Gini"), 3)) %>%
  group_by(kinase) %>%
  mutate(gini_kinase = round(ineq(100 - activity, type = "Gini"), 3)) %>%
  ungroup() %>%
  mutate(compound = as.character(compound)) %>%
  data.frame()

Encoding(kid$compound) <- "unknown"

save(kid, file = "./data/kid.RData")
