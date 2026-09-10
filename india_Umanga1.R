library(INLA)
library(maptools)
library(spdep)
library(rgeos);library(ggpubr); library(tidyverse)
library(INLAOutputs); library(INLA)
library(spdep); library(gridExtra)
library(sf); library(ggsn); 
library(GGally)

rm(list = ls())
cat("\08")

## Data
data<-read_csv("india5_fixed.csv")


# Map of India 
area <- st_read("./map/INDIA.shp")
area$ID <- str_to_title(area$ID)

ggplot(area) +
  geom_sf(size = .05,color="Grey") +
  geom_sf_label(aes(label = ID))+
  scale_fill_distiller(palette="PuRd", direction=1)+
  theme(text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10),
        legend.key.size = unit(10, "points"))
#__________________________________________________________________________________________#

# neighbour joining matrix 
india_nb <- poly2nb(as(area, "Spatial"), row.names = area$ID) #to create adjacency matrix 
nb2INLA("india_adj", india_nb)
india_adj <- "india_adj"

#making factor variable

data$cattledensity<-factor(data$cattledensity)
data$pden<-factor(data$pden)
data$gden<-factor(data$gden)
data$vet<-factor(data$vet)
data$inpro<-factor(data$inpro)
data$forest_c<-factor(data$forest_c)

data<-data%>%
  arrange(ID, year)%>%
  group_by(year) %>%
  mutate(expected_cases = totcattlepop * sum(outbreaknumber) / sum(totcattlepop),
         incidence_cases = outbreaknumber / totcattlepop * 1e5,
         rr_cases = outbreaknumber / expected_cases)%>%
  ungroup() %>%
  mutate(id =as.numeric(as.factor(ID)),
         id2 = as.numeric(as.factor(ID)),
         id3 = as.numeric(as.factor(ID)),
         id_year = year - min(year) + 1,
         id_year1 = year - min(year) + 1,
         id_year2 = year - min(year) + 1,
         id_year3 = year - min(year) + 1,
         id_id_year = 1:n())

## visual check of outbreak number##
data$ID<-as.factor(data$ID)
mp<-area%>%
  dplyr::select(ID)%>%
  left_join(data, by="ID")

ggplot(mp, aes(fill=outbreaknumber))+
  geom_sf(size=.5)+
  facet_wrap(~ year, nrow = 3)+
  scale_fill_viridis_c(option="C",name = "FMD")+
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(10,"points"))

sum(data$outbreaknumber == 0) / nrow(data)

##___________________Different models______________________________________________________________

##Priors 
prec <- list(prec = list(prior = "pc.prec", param = c(.3/.31, .01)))
prec_bym2 <- list(phi = list(prior = "pc", param = c(.5, 2/3)),
                  prec = list(prior = "pc.prec", param = c(.3/.31, .01)))

# Linear combination for spatio temporal interactions
lcs <- inla.make.lincombs(id_year = diag(length(unique(data$id_year))),
                          id_year2 = diag(length(unique(data$id_year2))))

RR01 <- inla(outbreaknumber ~ 1 +
               f(id,
                 model = "iid",
                 hyper = prec,
                 constr = TRUE) +
               f(id_year,
                 model = "iid",
                 hyper = prec,
                 constr = TRUE),
             family = "poisson",
             data = data,
             E = expected_cases,
             control.predictor = list(link = 1, compute = TRUE),
             control.compute = list(dic = TRUE), verbose = T)
summary(RR01)


qplot(data$outbreaknumber/data$expected_cases, RR01$summary.fitted.values$mean) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_smooth(method = "lm", color = "red", size = .5)

## Model checking## 

pred_p <- PredPValue(RR01)
pred_p$p_tails * 100 ## doing very bad


#Structured Model for Visceral[RR02] (Poisson)


RR02 <- inla(outbreaknumber ~ 1  +
               f(id,
                 model = "bym2",
                 hyper = prec_bym2,
                 graph = india_adj,
                 constr = TRUE,
                 scale.model = TRUE) +
               f(id_year,
                 model = "iid",
                 hyper = prec,
                 constr = TRUE) +
               f(id_year2,
                 model = "rw1",
                 hyper = prec,
                 constr = TRUE),
             family = "poisson",
             data = data,
             E = expected_cases,
             control.predictor = list(link = 1, compute = TRUE),
             control.compute = list(dic = TRUE), verbose = TRUE)

summary(RR02)

qplot(data$outbreaknumber/data$expected_cases, RR02$summary.fitted.values$mean) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_smooth(method = "lm", color = "red", size = .5)

## Model checking 

pred_p <- PredPValue(RR02)
pred_p$p_tails * 100 ## doing very bad

## Interaction 1: nu_i x phi_t
RR3 <- inla(outbreaknumber ~ 1 +
              f(id,
                model = "bym2",
                hyper = prec_bym2,
                graph = india_adj,
                scale.model = TRUE,
                constr = TRUE) +
              f(id_year,
                model = "iid",
                hyper = prec,
                constr = TRUE) +
              f(id_year2,
                model = "rw1",
                hyper = prec,
                constr = TRUE,
                scale.model = TRUE) +
              f(id_id_year,
                model = "iid",
                hyper = prec,
                constr = TRUE),
            family = "poisson",
            data = data,
            E = expected_cases,
            lincomb = lcs, 
            control.predictor = list(link = 1, compute = TRUE),
            control.compute = list(dic = TRUE, config=TRUE))
summary(RR3)

## Intercation 2: nu_i x gamma_t
RR4 <- inla(outbreaknumber ~ 1 +
              f(id,
                model = "bym2",
                hyper = prec_bym2,
                graph = india_adj,
                scale.model = TRUE,
                constr = TRUE) +
              f(id_year,
                model = "iid",
                hyper = prec,
                constr = TRUE) +
              f(id_year2,
                model = "rw1",
                hyper = prec,
                constr = TRUE,
                scale.model = TRUE) +
              f(id2,
                model = "iid",
                hyper = prec,
                constr = TRUE,
                group = id_year3,
                control.group = list(model = "rw1")),
            family = "poisson",
            data = data,
            E = expected_cases,
            lincomb = lcs,
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE))

## Interaction 3: phi_t x mu_i
RR5 <- inla(outbreaknumber ~ 1 +
              f(id,
                model = "bym2",
                hyper = prec_bym2,
                graph = india_adj,
                scale.model = TRUE,
                constr = TRUE) +
              f(id_year,
                model = "iid",
                hyper = prec,
                constr = TRUE) +
              f(id_year2,
                model = "rw1",
                hyper = prec,
                constr = TRUE,
                scale.model = TRUE) +
              f(id_year3,
                model = "iid",
                hyper = prec,
                constr = TRUE,
                group = id2,
                control.group = list(model = "besag", graph = india_adj)),
            family = "poisson",
            data = data,
            E = expected_cases,
            lincomb = lcs,
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE))

## Interaction 4: mu_i x gamma_t
RR6 <- inla(outbreaknumber ~ 1 +
              f(id,
                model = "bym2",
                hyper = prec_bym2,
                graph = india_adj,
                scale.model = TRUE,
                constr = TRUE) +
              f(id_year,
                model = "iid",
                hyper = prec,
                constr = TRUE) +
              f(id_year2,
                model = "rw1",
                hyper = prec,
                constr = TRUE,
                scale.model = TRUE) +
              f(id2,
                model = "besag",
                hyper = prec,
                graph = india_adj,
                group = id_year3,
                constr = TRUE,
                scale.model = TRUE,
                control.group = list(model = "rw1")),
            family = "poisson",
            data = data,
            E = expected_cases,
            lincomb = lcs,
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE))

## ________________________Model checking ___________________________________________________________

DIC(RR01, RR02, RR3, RR4, RR5, RR6)
pred_p <- PredPValue(RR01, RR02,RR3, RR5)
pred_p$p_tails * 100

## Explained variance
ExplainedVariance(RR3)

## model validation and further checks
data$muRR3 <- RR3$summary.fitted.values[,"mean"]
head(cbind(data$muRR3, data$outbreaknumber/data$expected_cases))## not bad

data$RR<-as.numeric(data$outbreaknumber/data$expected_cases)
datamicro2$RR<-datamicro2$RR%>%
  replace_na(0)
cor(data$muRR3,data$RR) # 0.93 getting near 1

#Plotting the calculated mean
spatio_temporal_fitted <- data %>%
  select(state, year) %>%
  mutate(RR = data$muRR3)

lepto2 <- mp %>%
  left_join(spatio_temporal_fitted, by = c("state", "year"))

ggplot(lepto2, aes(fill=RR))+
  geom_sf(size=.5)+
  facet_wrap(~ year, nrow = 5)+
  scale_fill_viridis_c(option="E",name = "mean")+
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(10,"points"))

#____________________________________________________
ERR1<- (data$RR-data$muRR3) / sqrt(data$muRR3)

# Residuals Vs parameters
plot(x = data$muRR3, 
     y = ERR1 ,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

plot(x = data$muRR3, # note that the PCA are explaning ok the zeros
     y = data$RR,
     xlab = "Fitted values",
     ylab = "Observed Outbreaks ",
     xlim = c(0, 12),
     ylim = c(0, 12))

plot(x = data$id, 
     y = ERR1,
     xlab = "Area",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)

#Priority Index calculation

pr_idx <- PriorityIndex(RR3, effect = "fitted", cutoff = 1, rescale_by = "year")
lepto <- data %>%
  mutate(pr_index = pr_idx)
lepto %>%
  select(ID, year, pr_index) %>%
  filter(year == 2008) %>%
  print(n = 32)

#Excess Risk calculation

exc_risk <- FittedExcess(RR3,cutoff = 1)

spatio_temporal_exc_risk <- data %>% 
  select(state, year) %>%
  mutate(ER = exc_risk) %>%
  spread(year, ER)

ggplot(lepto3, aes(fill = ER)) +
  geom_sf(size=.5)+
  facet_wrap(~ year, nrow = 3)+
  scale_fill_viridis_c(option="E",name = "RR")+
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(10,"points"))




























