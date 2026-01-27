rm(list=ls())
# Load libraries----
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggh4x)
library(glue)
library(tidyverse)
# Data setup ----
fwlist <- c("Plowmanetal", "Krkonose", "CameroonW", "CameroonD")
fwlist2 <- c("Myrmecophytic\nant-plant dataset",
             "Temperate\nplant-pollinator dataset",
             "Tropical wet-season\nplant-pollinator dataset",
             "Tropical dry-season\nplant-pollinator dataset")

fwlen <- c(4, 4, 4, 4)

elevations_plowman <- c(700, 935, 1170, 1400)
elevations_krkonose <- c(450, 600, 800, 1000)
elevations_cameroon <- c(650, 1100, 1450, 2200)
names(elevations_cameroon) = names(elevations_krkonose) = names(elevations_plowman) = c(1,2,3,4)

elevations_list <- list(
  elevations_plowman,
  elevations_krkonose,
  elevations_cameroon,
  elevations_cameroon
)
names(elevations_list) <- fwlist2

# Read in datasets----
datasets <- 
  fwlist %>% 
  map(~ read_csv(paste0("realfood/netsgrad1000noi_", ., "x1.csv")))
# Set names
names(datasets) <- fwlist2

# bind rows, pivot longer
datasets <-
  datasets %>%
  bind_rows(.id = "dataset_name") %>% 
  mutate(
    gr = map2_dbl(dataset_name, gr, ~ elevations_list[[.x]][as.character(.y)])
  )  %>%
  mutate(dataset_name = 
           factor(dataset_name, 
                  levels = fwlist2)) %>% 
  pivot_longer(NOconnections:vulnerability.LL,
               names_to = "index",
               values_to = "index_value") %>% 
  mutate(index = case_when(index == "number.of.species.LL" ~ "Species~richness~LO",
                           index == "number.of.species.HL" ~ "Species~richness~HI",
                            index == "NOconnections" ~ "Total~realised~links",
                           #index == "SA" ~ "Specialisation~asymmetry",
                           index == "generality.HL" ~ "Generality",
                           index == "vulnerability.LL" ~ "Vulnerability",
                           index == "H2" ~ "H[2]'",
                           index == "nestedness" ~ "NestednessOld",
                           index == "NODF" ~ "Nestedness",
                           index == "connectance" ~ "Connectance",
                           .default = index)) %>% 
  mutate(index = factor(index, 
                        levels = 
                          c("Species~richness~LO","Species~richness~HI",
                            "Total~realised~links",
                            "Connectance",
                            "Nestedness",
                            #"NODF",
                            #"Specialisation~asymmetry",
                            "H[2]'",  # Apply bquote() for H2 in fct_relevel
                            "Generality",
                            "Vulnerability"),
                        labels = 
                          c("Species~richness~LO","Species~richness~HI",
                            "Total~realised~links",
                            "Connectance",
                            "Nestedness",
                            #                          "NODF",
                            #"Specialisation~asymmetry",
                            bquote(italic("H"[2]*"'")),  # Apply bquote() for H2 in fct_relevel
                            "Generality",
                            "Vulnerability")))

datasets <- datasets %>% filter(!is.na(index))

summary_nulls_1 <- 
  datasets %>% 
  filter(null != 0) %>% 
  group_by(dataset_name, gr, index) %>% 
  summarise(lower_quantile = quantile(index_value, 0.025, na.rm = T),
            upper_quantile = quantile(index_value, 0.95, na.rm = T),
            sd = sd(index_value, na.rm = T),
            mean = mean(index_value, na.rm = T))

summary_nulls_2 <- 
  summary_nulls_1 %>% 
  left_join(datasets %>% filter(null != 0)) %>% 
  group_by(dataset_name, index, gr, null) %>% 
  summarise(ses_null_gr = (abs(index_value - mean)/abs(sd))) %>% 
  distinct() %>% 
  group_by(dataset_name, index, null) %>% 
  summarise(ses_null = sum(ses_null_gr)) %>% 
  summarise(upper_quantile_ses = quantile(ses_null, 0.95, na.rm = T),
            lower_quantile_ses = quantile(ses_null, 0.025, na.rm = T))

index_points <- datasets %>% 
  filter(null == 0) %>% 
  left_join(summary_nulls_1) %>% 
  group_by(dataset_name, index, gr, null) %>% 
  mutate(index_shape = case_when(index_value > upper_quantile |
                                   index_value < lower_quantile ~ "2",
                                 .default = "1"))

stars_true <-
  datasets %>% 
  filter(null == 0) %>% 
  left_join(summary_nulls_1) %>% 
  group_by(dataset_name, index, gr, null) %>% 
  summarise(ses_true_gr = (abs(index_value - mean)/abs(sd))) %>% 
  distinct() %>% 
  group_by(dataset_name, index, null) %>% 
  summarise(ses_true_gr = sum(ses_true_gr)) %>% 
  left_join(summary_nulls_2) %>% 
  mutate(sig = 
           case_when(ses_true_gr > upper_quantile_ses |
                       ses_true_gr < lower_quantile_ses ~ "*",
                     .default = "")) %>% 
  distinct(dataset_name, index, sig)




plot1 <- 
  datasets %>% 
  filter(null != 0) %>% 
  group_by(dataset_name, gr, index) %>% 
  summarise(lower_quantile = quantile(index_value, 0.025, na.rm = T),
            upper_quantile = quantile(index_value, 0.95, na.rm = T)) %>% 
  
  ggplot(aes(x = gr, fill = index, colour = index)) +
  
  # Confidence intervals
  geom_ribbon(aes(ymin = lower_quantile, ymax = upper_quantile), alpha = 0.25, colour = NA) +
  geom_line(aes(y = (upper_quantile-(upper_quantile - lower_quantile)/2)),
            # alpha = 0.75,
            linewidth = 1,
            linetype = "dotted") +
  
  # True points
  geom_point(data = index_points,
             aes(x = gr, y = index_value, shape = index_shape), colour = "black",
             size = 3) +
  scale_shape_manual(values = c("1" = 1, "2" = 24)) +
  geom_line(data = index_points,
            aes(x = gr, y = index_value)) +
  
  # Significance stars
  geom_text(data = stars_true, aes(x = Inf, y = Inf, label = sig),
            hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, size = 12,
            colour = "red") +
  
  # Faceting
  facet_grid2(index ~ dataset_name,
              labeller = labeller(
                index = label_parsed,  # Parse 'index' variable
                dataset_name = label_value),  # Leave 'dataset_name' as regular text
              scales = "free",
              independent = "y",
              switch = "y",
              axes = "all") +
  coord_cartesian(clip = 'off') +
  # Theming
  theme_classic(base_size = 15) +
  theme(
    plot.margin = unit(c(0, 1, 1, 0), "lines"),  # Widens the right margin
    legend.position = "bottom",
    strip.placement = "outside",
    strip.switch.pad.grid = unit(0, "cm"),
    strip.background = element_blank(), 
    strip.text.y.left = element_text(angle = 0),
    panel.spacing.x = unit(0.5, "cm"),
    # panel.spacing.y = unit(0.5, "cm")
  ) +
  labs(x = "Elevation (m a.s.l.)", y = "")
plot1<-plot1+scale_shape(name="Within null model",labels=c("Yes","No (     for whole gradient)"))+
  guides(fill="none")+guides(colour="none")
plot1


ggsave(filename = paste0("realfw2_EL", ".pdf"), plot = plot1, device = NULL, width = 12, height = 10.69, dpi = 300)


fwlist<-c("Plowmanetal","Krkonose","CameroonW","CameroonD")
fwlist2<-c("Mount Wilhelm","Krkonose","Mount Cameroon Wet","Mount Cameroon Dry")
fwlist2 <- c("Myrmecophytic\nant-plant dataset",
             "Temperate\nplant-pollinator dataset",
             "Tropical wet-season\nplant-pollinator dataset",
             "Tropical dry-season\nplant-pollinator dataset")

fwlen<-c(4,4,4,4)
dohist<-TRUE
pdf(paste0("noi_allfws6",ifelse(dohist,"_hist",""),"2_EL.pdf"),width=11,height=12)
#par(mfrow=c(8,5))
nf <- layout( matrix(c(1:40), nrow=8,ncol=5, byrow=TRUE),width=c(1,3,3,3,3),
              heights=c(1,rep(3,6),1))
wcols<-c("black",rainbow(8))
#names(wnetsdall)
i=1;fi=1
#wnetsdall$nestedness<-wnetsdall$nestedness/max(wnetsdall$nestedness,na.rm=T)
#wnetsdall$generality.HL<-wnetsdall$generality.HL/max(wnetsdall$generality.HL,na.rm=T)
#wnetsdall$vulnerability.LL<-wnetsdall$vulnerability.LL/max(wnetsdall$vulnerability.LL,na.rm=T)
#wnetsdall$NOconnections<-wnetsdall$NOconnections/max(wnetsdall$NOconnections,na.rm=T)
#8
par(mar=c(2,4,2,0.5))
plot(1:10,typ="n",bty='n',axes=F,xlab="",ylab="")
for (fi in 1:length(fwlist) ) {
  plot(1:10,typ="n",bty='n',axes=F,xlab="",ylab="")
  (onfw=fwlist[fi])
  mtext(fwlist2[grep(onfw,fwlist)],3,-2,cex=1.1)
}
#names(wnetsdall)[1:8]<-c("Total realized links","Connectance","nest","Nestedness",
#                    "Specialization asymmetry", "H2", "Generality","Vulnerability")
for (i in c(1:2,4,6:8) ){
  plot(1:10,typ="n",bty='n',axes=F,xlab="",ylab="")
  names(wnetsdall)[1:8]<-c("Total realized links","Connectance","nest","Nestedness",
                           "Specialization asymmetry", "H2", "Generality","Vulnerability")
  mtext(names(wnetsdall)[i],2,0,cex=1.2)
  for (fi in 1:length(fwlist) ) {
    wnetsdall<-read.csv(paste0("realfood/netsgrad1000noi_",fwlist[fi],"x1.csv"))
    
    mygr<-fwlen[fi]
    grale<-mygr
    (onfw=fwlist[fi])
    
    ane<-aggregate(as.formula(paste0(names(wnetsdall)[i],"~gr")),FUN=mean,data=wnetsdall[wnetsdall$null !=0,],na.rm=T )
    anesd<-aggregate(as.formula(paste0(names(wnetsdall)[i],"~gr")),FUN=sd,data=wnetsdall[wnetsdall$null !=0,],na.rm=T )
    anel<-aggregate(as.formula(paste0(names(wnetsdall)[i],"~gr")),FUN=quantile,data=wnetsdall[wnetsdall$null !=0,],na.rm=T,probs=c(0.025) )
    aneh<-aggregate(as.formula(paste0(names(wnetsdall)[i],"~gr")),FUN=quantile,data=wnetsdall[wnetsdall$null !=0,],na.rm=T,probs=c(0.975) )
    ane0<-wnetsdall[wnetsdall$null ==0,]
    #random ses
    rasum<-c()
    for (oneg in unique(wnetsdall[wnetsdall$null>0,]$null) ) {
      anex<-wnetsdall[wnetsdall$null == oneg,]
      rasum <- c(rasum, sum( abs( ane[,2] - anex[,i] )/abs(anesd[,2]) )   )
    }
    
    (rasumq<-quantile(rasum,probs=c(0.025,0.95), na.rm=T))
    # get sum of SES
    (obssum <- sum( abs( ane[,2] - ane0[,i])/abs(anesd[,2]) ))
    
    if (! dohist) {
      plot( y=seq(min(wnetsdall[,i],na.rm=T), max(wnetsdall[,i],na.rm=T),length=grale ),
            x=1:grale,col="white",ylab=ifelse(fi==12,names(wnetsdall)[i],""),xlab="",xaxt='n',
            main=ifelse(i==12,fwlist2[grep(onfw,fwlist)],"") )
      axis(1, sort(unique(wnetsdall$gr)))
      lines(x=ane$gr,y=ane[,2],col=wcols[i],lwd=1,lty=3)
      lines(x=anel$gr,y=anel[,2],col=wcols[i],lwd=1)
      lines(x=aneh$gr,y=aneh[,2],col=wcols[i],lwd=1)
      lines(x=ane0$gr,y=ane0[,i],col=wcols[i],lwd=3,lty=2)
      points(x=ane0$gr,y=ane0[,i],col=wcols[i],lwd=2,lty=2,
             cex= ifelse( (ane0[,i]<anel[,2]) | (ane0[,i]>aneh[,2]),3,1 ),
             pch=ifelse( (ane0[,i]<anel[,2]) | (ane0[,i]>aneh[,2]),1,2 ))
    } else {
      #if( any((ane0[,i]<anel[,2]) | (ane0[,i]>aneh[,2])) ) mtext("*",side=3,line=-2,adj=1,cex=3,col=wcols[i])
      #if ((i==1)& (fi==1))
      #  legend("topright",lty=c(2,1,3),lwd=c(3,1,1),legend=c("observed","95% conf ints","mean conf ints"))
      
      hist(rasum,main="",xlab="")
      abline(v=obssum,col="red",lwd=2)
      abline(v=rasumq,col="gray",lwd=2)
      mtext(paste0( "Observed sum=",round(obssum,3)),3,1,col="red",adj=0)
    }
    
    if( any( (obssum < rasumq[1]) | (obssum > rasumq[2])  ) ) mtext("*",side=3,line=-2,adj=1,cex=3,col="red")
  }
  #if ((i==1)& (fi==2))
  #  legend("bottomleft",legend=row.names(wnetsd),col=wcols,lty=1,bg =NA)
}
plot(1:10,typ="n",bty='n',axes=F,xlab="",ylab="")
for (fi in 1:length(fwlist) ) {
  plot(1:10,typ="n",bty='n',axes=F,xlab="",ylab="")
  (onfw=fwlist[fi])
  #mtext("elevation",3,0,cex=1.1)
}

dev.off()

