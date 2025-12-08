library(tidyverse)
# Reads in real food web data, then calls the rf_nulls function to generate nulls
rm(list=ls())
pld <- TRUE  # do we run bipartite?
ran <- TRUE  # do we run null models?

# List of Food Webs
fwlist <- c("Plowmanetal", "Krkonose", "CameroonW", "CameroonD")  # names
fwlen <- c (4, 4, 4, 4)  # gradient lengths

###################################
# Code specific to each dataset----
## plowman----
{#wmt<-read.csv("../simfood/dataIN/Antplantdata.csv")
  wmt <- read.csv("realfood/Antplantdata.csv")
  wmt$species_lo <- wmt$Species
  wmt$species_hi <- ifelse(wmt$Species.code == "UNC", "UNOCCUPIED", wmt$Species.code)
  wmt$elevation <- wmt$Elevation
  wmt$elevation<-as.factor(wmt$elevation)
  wmt$elevation<-as.numeric(wmt$elevation)
  wmt$elevation<-ifelse (wmt$elevation %in% 1:2, 1, 
                         ifelse (wmt$elevation %in% 3:4, 2,
                                 ifelse (wmt$elevation %in% 5:6,3,4)))
  table(wmt$elevation)
  
  wmt <- wmt[, ncol(wmt) - c(0, 1, 2)]
  wmt <- wmt[wmt$species_hi != "UNOCCUPIED", ]
  onfw <- "Plowmanetal"
  mygr1 <- mygr <- length(unique(wmt$elevation))
  head(wmt)
  plowman <- wmt
  
  # sanity checks
  length(unique(wmt$species_hi))
  length(unique(wmt$species_lo))
  }
plowman %>% 
  pivot_longer(species_hi:species_lo) %>% 
  group_by(name, value) %>% 
  summarise(n_elev = n_distinct(elevation)) %>% 
  ggplot(aes(x = n_elev, group = name)) +
  geom_histogram() +
  facet_wrap(.~name) +
  theme_bw()

# run null models
tictoc::tic("Running null models")
source("rf_nulls.R")
tictoc::toc()
#GOTO LINE 200 and continue
###################################
# robert
##cameroon & krkonose 4 elev
##Kohler  & Arroyo 3 elev
## Krkonose----
onfw <- "Krkonose"
mygr = 4
mygr1 <- mygr
#wmtx <- read.csv("dataset/Robert/Krkonose_data.csv")
krpre <- "realfood/Krkonose_data"
w1 <- read.csv(paste0(krpre, "1.csv"))
w2 <- read.csv(paste0(krpre, "2.csv"))
w2 <- w2[-c(nrow(w2)), ]
w3 <- read.csv(paste0(krpre, "3.csv"))
w4 <- read.csv(paste0(krpre, "4.csv"))
#wmtx$Species.name <- ifelse( wmtx$Species.Code..plants. != "",wmtx$Species.Code..plants., wmtx$Species.name )
#table(wmtx$LEVEL)
#wmt<-data.frame(species_lo=NA,species_hi=NA,elevation=NA)
#for (i in grep("LOWER",wmtx$LEVEL) ) for (y in 1:4) 
#  if ( (!is.na(wmtx[i,grep("site",names(wmtx))[y]])) & (wmtx[i,grep("site",names(wmtx))[y]] >0)  ) wmt<-rbind(wmt, c(wmtx[i,]$Species.name,NA,y))
#for (i in grep("UPPER",wmtx$LEVEL) ) for (y in 1:4) 
#  if ( (!is.na(wmtx[i,grep("site",names(wmtx))[y]])) & (wmtx[i,grep("site",names(wmtx))[y]] >0)  ) wmt<-rbind(wmt, c(NA,wmtx[i,]$Species.name,y))
library(data.table)
wmt <- melt(setDT(w1), variable.name = "species_hi")
names(wmt)[1] <- "species_lo"
wmt$elevation <- 1
head(wmt)
wmt2 <- melt(setDT(w2), variable.name = "species_hi")
names(wmt2)[1] <- "species_lo"
wmt2$elevation <- 2
wmt3 <- melt(setDT(w3), variable.name = "species_hi")
names(wmt3)[1] <- "species_lo"
wmt3$elevation <- 3
wmt4 <- melt(setDT(w4), variable.name = "species_hi")
names(wmt4)[1] <- "species_lo"
wmt4$elevation <- 4
wmt <- rbind(wmt, wmt2, wmt3, wmt4)
wmt <- wmt[!is.na(wmt$value), ]
wmt$elevation <- as.numeric(wmt$elevation)
krkonose <- wmt
# sanity check
length(unique(wmt$species_hi))
length(unique(wmt$species_lo))
krkonose %>%  
  as.data.frame() %>% 
  pivot_longer(species_hi:species_lo, values_to = "spec") %>% 
  group_by(name, spec) %>% 
  summarise(n_elev = n_distinct(elevation)) %>% 
  ggplot(aes(x = n_elev, group = name)) +
  geom_histogram() +
  facet_wrap(.~name) +
  theme_bw()
}
# run null models
tictoc::tic("Running null models")
source("rf_nulls.R")
tictoc::toc()
#GOTO LINE 200 and continue
###################################
## cameroon wet ----
mygr = 4
mygr1 <- mygr
onele <- c("CL", "DG", "MS", "PC")
on = "W"
onfw <- paste0("Cameroon", on)
robpre <- "realfood/c"
w1 <- read.csv(paste0(robpre, onele[1], on, ".csv"))
w2 <- read.csv(paste0(robpre, onele[2], on, ".csv"))
w3 <- read.csv(paste0(robpre, onele[3], on, ".csv"))
w4 <- read.csv(paste0(robpre, onele[4], on, ".csv"))

library(data.table)
wmt <- melt(setDT(w1), variable.name = "species_hi")
names(wmt)[1] <- "species_lo"
wmt$elevation <- 1
head(wmt)
wmt2 <- melt(setDT(w2), variable.name = "species_hi")
names(wmt2)[1] <- "species_lo"
wmt2$elevation <- 2
wmt3 <- melt(setDT(w3), variable.name = "species_hi")
names(wmt3)[1] <- "species_lo"
wmt3$elevation <- 3
wmt4 <- melt(setDT(w4), variable.name = "species_hi")
names(wmt4)[1] <- "species_lo"
wmt4$elevation <- 4
wmt <- rbind(wmt, wmt2, wmt3, wmt4)
wmt <- wmt[!is.na(wmt$value), ]
wmt$elevation <- as.numeric(wmt$elevation)
head(wmt)
cameroonW <- wmt
# sanity check
length(unique(wmt$species_hi))
length(unique(wmt$species_lo))
cameroonW %>% as.data.frame() %>% 
  pivot_longer(species_hi:species_lo, values_to = "spec") %>% 
  group_by(name, spec) %>% 
  summarise(n_elev = n_distinct(elevation)) %>% 
  ggplot(aes(x = n_elev, group = name)) +
  geom_histogram() +
  facet_wrap(.~name) +
  theme_bw()
}

# run null models
tictoc::tic("Running null models")
source("rf_nulls.R")
tictoc::toc()


#GOTO LINE 200 and continue
## cameroon dry----
{onele <- c("CL", "DG", "MS", "PC")
on = "D"
mygr1 <- mygr<-4
onfw <- paste0("Cameroon", on)
w1 <- read.csv(paste0(robpre, onele[1], on, ".csv"))
w2 <- read.csv(paste0(robpre, onele[2], on, ".csv"))
w3 <- read.csv(paste0(robpre, onele[3], on, ".csv"))
w4 <- read.csv(paste0(robpre, onele[4], on, ".csv"))

library(data.table)
wmt <- melt(setDT(w1), variable.name = "species_hi")
names(wmt)[1] <- "species_lo"
wmt$elevation <- 1
head(wmt)
wmt2 <- melt(setDT(w2), variable.name = "species_hi")
names(wmt2)[1] <- "species_lo"
wmt2$elevation <- 2
wmt3 <- melt(setDT(w3), variable.name = "species_hi")
names(wmt3)[1] <- "species_lo"
wmt3$elevation <- 3
wmt4 <- melt(setDT(w4), variable.name = "species_hi")
names(wmt4)[1] <- "species_lo"
wmt4$elevation <- 4
wmt <- rbind(wmt, wmt2, wmt3, wmt4)
wmt <- wmt[!is.na(wmt$value), ]
wmt$elevation <- as.numeric(wmt$elevation)
head(wmt)
cameroonD <- wmt
# sanity checks
length(unique(wmt$species_hi))
length(unique(wmt$species_lo))
cameroonD %>% as.data.frame() %>% 
  pivot_longer(species_hi:species_lo, values_to = "spec") %>% 
  group_by(name, spec) %>% 
  summarise(n_elev = n_distinct(elevation)) %>% 
  ggplot(aes(x = n_elev, group = name)) +
  geom_histogram() +
  facet_wrap(.~name) +
  theme_bw()

}
# run null models
tictoc::tic("Running null models")
source("rf_nulls.R")
tictoc::toc()