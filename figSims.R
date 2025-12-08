library(tidyverse)
library(broom)
library(ggrepel)
library(ggh4x)
library(glue)
rm(list=ls())
speclo<-c(2,20,50)
spechi<-c(20,50,90)
fun2run<-mean
onstep<-10
ffo<-"simout"

polyn<-4
grfilepre<-"mdeFeas2_12345_20_NOSPEC"

# if to show second or fourth polynom
# if to show random or increasing specialization
for (polyn in c(2,4)) 
  for(grfilepre in c("mdeFeas2_12345_20_NOSPEC","mdeFeas2_12345_20_DEC") ){


# res 9–4 in one 
#"res103"
(list.files(ffo, pattern = "*.csv$"))

#grfilepre<-"mdeFeas2_12345_20_NOSPEC"
#grfilepre<-"mdeFeas_12345_20_DEC"

(list.files(ffo, pattern = paste0("^", grfilepre, ".*", "500reps", "*.csv$")))
#(list.files(ffo, pattern = paste0("^", grfilepre, ".*", "100reps", "*.csv$")))

# HERE SELECT WHICH FILES TO DRAW MOSTLY ITS is 1:2 or 3:4
# filera<-1:2#3:4 #2:5 #7:10
filera<-1:3
(fili <- list.files(ffo, pattern = paste0("^", grfilepre, ".*", "500reps", "*.csv$"), full.names = T))
#(fili <- list.files(ffo, pattern = paste0("^", grfilepre, ".*", "100reps", "*.csv$"), full.names = T))

# (fili<-list.files(ffo,pattern = "*.csv$", full.names = T)[filera])
minlofoc<-speclo
maxlo<-spechi

# Read in files----

fili2 <- data.frame(range = 
                      paste("Range span 1–",
                            round(100 / 
                                    as.numeric(
                                      gsub("(.*)grad_(.*)rangered(.*)","\\2",fili),0))),
                    id = as.character(1:length(fili))) # id for joining

fili2 <- data.frame(range = as_factor(c("Range span 1–100", "Range span 1–50", "Range span 1–25")),
                    id = as.character(1:length(fili)))

data <- 
  fili %>% 
  lapply(read_csv) %>% 
  bind_rows(.id = "id") %>% 
  left_join(fili2)

data <-
  data %>% 
  rename("old_nestedness" = "nestedness")

#data$spec_asym <- (data$HIints-data$LOints)/(data$HIints+data$LOints)
#data <-
#  data %>% mutate("sa" = (HIints-LOints)/(HIints+LOints))  


data <-
  data %>% 
  rename("dist" = "gr",
         "Species richness LO" = "LOints",
         "Species richness HI" = "HIints",
         "Total realised links" = "allInts",
         #         "Specialisation asym." = "spec_asym",
         "nestedness" = "NODF")

table(data$minlo)
table(data$grstep)
table(data$id)
datax<-data[ (data$minlo ==20) & (data$id ==1) & (data$grstep==10), ]
names(datax)
#plot(`Specialisation asymmetry`~dist,data=datax)
#anova(lm(H2~poly(dist),data=datax))
#cor(datax$`Total realised links`,datax$)

lm1<-lm(H2~dist,data=datax)
lm2<-lm(H2~poly(dist,2),data=datax)

#plot(`Specialisation asym.`~ dist,data=datax)

#datax$sa<-(datax$`Species richness LO`-datax$`Species richness HI`)/
#  (datax$`Species richness LO`+datax$`Species richness HI`)
#plot(sa~ dist,data=datax)

#plot(`Specialisation asym.`~ dist,data=datax)
#plot(sa~ dist,data=datax)

#lm1<-lm(connectance~dist,data=datax)
#lm2<-lm(connectance~poly(dist,2),data=datax)

c(AIC(lm1), AIC(lm2)) - min(c(AIC(lm1), AIC(lm2)))


models <- data %>% 
  group_by(minlo, dist, range) %>% 
  filter(grstep == 10) %>% 
  filter(dist %in% seq(10,100,10)) %>% 
  # summarise(across("Total realised links":vulnerability, mean)) %>%
  mutate(Connections = `Total realised links`) %>% 
  pivot_longer(`Total realised links`:vulnerability, names_to = "index") %>% 
  filter(index %in% c("Species richness LO", "Species richness HI","Total realised links", 
                      "connectance",
                      "nestedness",
                      "NODF",
                      #                      "Specialisation asym.",
                      "H2",
                      "generality",
                      "vulnerability")) %>% 
  mutate(index = factor(index,
                        levels = 
                          c("Species richness LO","Species richness HI",
                            "Total realised links", 
                            "connectance",
                            "nestedness",
                            "NODF",
                            #                             "Specialisation asym.",
                            "H2",
                            "generality",
                            "vulnerability"),
                        labels = 
                          c("Species~richness~LO","Species~richness~HI",
                            "Total~realised~links",
                            "Connectance",
                            "Nestedness",
                            "NODF",
                            #                              "Specialisation~asym.",
                            bquote(italic("H"[2]*"'")),  # Apply bquote() for H2 in fct_relevel
                            "Generality",
                            "Vulnerability")
  )) %>%  
  filter(!is.na(value)) %>% 
  group_by(range, index, minlo)  %>%
  nest()
#rm(data)



# z scores
#i=1
#for (i in 1:nrow(models)){
#  models[i,]$data[[1]]<- models[i,]$data[[1]] %>% 
#    mutate( value= (value-mean(value,na.rm=T))/sd(value,na.rm=T))
#}

# remove outliers
#i=1
#for (i in 1:nrow(models)){
#  models[i,]$data[[1]]<- models[i,]$data[[1]] %>% 
#    mutate( IQR=IQR(value),value= ifelse( value < quantile(value, probs=c( .975), na.rm = FALSE)+1.5*IQR |
#              value > quantile(value, probs=c( .025), na.rm = FALSE)-1.5*IQR,value,NA) ) 
#}



models <- 
  models %>% 
  mutate(
    fit_lm1 = map(data, ~lm(value ~ dist, data = drop_na(., value, dist))),  # Linear model
    fit_lm2 = map(data, ~lm(value ~ poly(dist, polyn), data = drop_na(., value, dist)))  # Cubic model
  ) %>% 
  mutate(anova_result = list(anova(fit_lm1[[1]], fit_lm2[[1]]))    # ANOVA between the two models
  ) %>% 
  mutate(aic_result = list( which.min(c( AIC(fit_lm1[[1]]), AIC(fit_lm2[[1]]))) )    # ANOVA between the two models       
  ) %>%
  mutate(
    p_value = sapply(anova_result, function(x) x$`Pr(>F)`[2]),  # Extract p-value from ANOVA result
    f_stat_anova = sapply(anova_result, function(x) x$`F`[2]),  # Get the F-statistic for the comparison
    coef_lm1 = list(coef(fit_lm1[[1]])),  # Coefficients from the linear model
    coef_lm2 = list(coef(fit_lm2[[1]]))   # Coefficients from the quadratic model
  ) %>%
  mutate(model = case_when( ((p_value <= 0.05) & (aic_result==2)) ~ "lm2", .default = "lm1"),
         letter = 
           case_when(model == "lm2" & coef_lm2[[1]][3] <= 0 ~ "H",
                     model == "lm2" & coef_lm2[[1]][3] > 0 ~ "U",
                     model == "lm1" & coef_lm1[[1]][2] <= 0 ~ "Dec.",
                     .default = "Inc.")) 

model1 <- 
  models %>%
  mutate(
    augmentedlm1 = map2(fit_lm1, data, ~broom::augment(.x, newdata = .y)),
    augmentedlm2 = map2(fit_lm2, data, ~broom::augment(.x, newdata = .y))) %>% 
  mutate(correct_augment = 
           case_when(model == "lm1" ~ augmentedlm1, .default = augmentedlm2))

(plot1 <- 
    model1 %>% 
    unnest(correct_augment) %>% 
    group_by(minlo, range, index, dist, .fitted, letter) %>% 
    summarise(Connections = mean(Connections, na.rm = TRUE)) %>% 
    rename(fitted = .fitted) %>% 
    group_by(minlo, index, range) %>%
    mutate(model_label = 
             case_when(dist == max(dist) ~ letter, 
                       TRUE ~ NA_character_)) %>%  
    ggplot(aes(x = dist, y = fitted, group = as_factor(minlo))) +
    geom_line(aes(colour = as_factor(minlo), linetype = as_factor(minlo)), linewidth = 1.5,show.legend = T) +
    #ggtext(x=rep(100,3),label=letter)+
    # geom_point(aes(shape = as_factor(minlo), fill = Connections),
    #            colour = "white", size = 3) +
    geom_text_repel(
      aes(label = model_label, colour = as_factor(minlo)),
      # colour = "navy",
      # hjust = 0,
      size = 4,show.legend = F,
      xlim = c(106, 125),
      direction = "y"
    ) +coord_cartesian(clip = "off") +
    # scale_fill_gradient(low = "blue", high = "red") +  # Adjust colours as needed
    # scale_shape_manual(values = c(21, 22, 23)) +
    
    # scale_colour_manual(values = c("black", "grey", "orange"), name = "Text Labels") +  # For text
    # guides(
    #   colour = guide_legend(override.aes = list(linetype = 0, size = 5))  # Separate legend for text
    # ) +
    facet_grid2(index ~ range,
                labeller = labeller(index = label_parsed,
                                    range = function(x) str_wrap(x, width = 10)),
                scales = "free",
                independent = "y", switch = "y",
                axes = "all") +
    # coord_cartesian(xlim = c(5, 100), clip = 'off') + # Focuses the x-axis
    #    theme_classic(base_size = 12) +
    #    geom_hline(yintercept = 0.5 + 3, color = "red", linewidth = 1) 
    
    theme(
      plot.margin = unit(c(0.5, 2, 0.5, 0.5), "lines"),  # Widens the right margin
      text = element_text(size = 12),
      #      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 10),
      #      legend.title = element_text(size = 12),
      legend.position = "bottom",
      strip.placement = "outside",
      strip.switch.pad.grid = unit(0.25, "cm"),
      strip.background = element_blank(), 
      strip.text.y.left = element_text(angle = 0),
      panel.spacing.x = unit(0.5, "cm")
    ) 
)


(plot1 <- plot1 + labs(x ="Artificial gradient",y = "")+
    scale_colour_viridis_d(name="Specialisation",labels=c("High","Middle","Low"))+ 
    scale_linetype(name="Specialisation",labels=c("High","Middle","Low")) )

ggsave(filename = paste0(grfilepre, "4x_poly",polyn,"_500_EL.pdf"), plot = plot1, device = NULL, width = 8, height = 10.69, dpi = 300)
}