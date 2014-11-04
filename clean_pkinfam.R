rm(list = ls())

library(dplyr)
library(stringr)

if (!file.exists("./raw/pkinfam.txt")) {
  download.file(url = "http://www.uniprot.org/docs/pkinfam.txt", 
                destfile = "./raw/pkinfam.txt", 
                method = "curl")
}

read_kinase_family <- function(start, end, tag, file = "./raw/pkinfam.txt") {
  family <- read.table("./raw/pkinfam.txt", skip = start, nrow = end, 
                       stringsAsFactors = FALSE, fill = TRUE)
  family %>%
    select(kinase = V1, accession = V3) %>%
    filter(!(grepl(")", accession))) %>%
    mutate(accession = str_sub(accession, start = 2)) %>% 
    mutate(family = tag)
}

agc <- read_kinase_family(83, 141-83, "agc")
camk <- read_kinase_family(149, 230-149, "camk")
ck1 <- read_kinase_family(237, 249-237, "ck1")
cmgc <- read_kinase_family(256, 318-256, "cmgc")
nek <- read_kinase_family(325, 336-325, "nek")
ste <- read_kinase_family(343, 399-343, "ste")
tkl <- read_kinase_family(406, 440-406, "tkl")
tk <- read_kinase_family(447, 540-447, "tk")
other <- read_kinase_family(547, 614-547, "other")
adck <- read_kinase_family(621, 626-621, "adck")
alpha <- read_kinase_family(633, 639-633, "alpha")
fast <- read_kinase_family(646, 647-646, "fast")
pdk <- read_kinase_family(654, 659-654, "pdk")
pi <- read_kinase_family(666, 673-666, "pi")
rio <- read_kinase_family(680, 683-680, "rio")

kinases <- rbind_list(agc, camk, ck1, cmgc, nek, ste, tkl, tk, other, adck, alpha, fast, pdk, pi, rio)

uniprot <- list()
for (i in seq(nrow(kinases))) {
  accession <- kinases$accession[i]
  url <- paste0("http://www.uniprot.org/uniprot/?query=accession:",
                accession, 
                "&format=tab&columns=id,genes,protein%20names")
  record <- read.table(url, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  uniprot[[i]] <- record
  Sys.sleep(0.3)
}

uniprot[[106]] <- data.frame(Entry = "Q13131", 
                             Gene.names = "PRKAA1 AMPK1", 
                             Protein.names = "5'-AMP-activated protein kinase catalytic subunit alpha-1 (AMPK subunit alpha-1) (EC 2.7.11.1) (Acetyl-CoA carboxylase kinase) (ACACA kinase) (EC 2.7.11.27) (Hydroxymethylglutaryl-CoA reductase kinase) (HMGCR kinase) (EC 2.7.11.31) (Tau-protein kinase PRKAA1) (EC 2.7.11.26)", 
                             stringsAsFactors = FALSE)
uniprot[[107]] <- data.frame(Entry = "P54646", 
                             Gene.names = "PRKAA2 AMPK AMPK2", 
                             Protein.names = "5'-AMP-activated protein kinase catalytic subunit alpha-2 (AMPK subunit alpha-2) (EC 2.7.11.1) (Acetyl-CoA carboxylase kinase) (ACACA kinase) (EC 2.7.11.27) (Hydroxymethylglutaryl-CoA reductase kinase) (HMGCR kinase) (EC 2.7.11.31)", 
                             stringsAsFactors = FALSE)

uniprot <- rbind_all(uniprot)
kinases_human <- kinases %>%
  left_join(uniprot, by = c("accession" = "Entry")) %>%
  rename(gene = Gene.names, protein = Protein.names) %>%
  mutate(protein = str_replace(protein, " \\(.*$", ""))

dir.create(file.path("./data"), showWarnings = FALSE)
save(kinases_human, file = "./data/kinases_human.RData")