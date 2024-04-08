# PotatoMETAbiom Field trial experiment
# NL - Structural Equation Modeling analysis (SEM) - Week5

#library####
library(nlme)
library(lme4)
library(piecewiseSEM)
library(readxl)
library(tidyverse)
#library(lavaan),FIML
library(ggplot2)
library(semPlot)# save model


#Belowground biomass####

##All management####

#use numeric for variety and treatment name. add one column "block" in data for random = list(~ 1 |Plot,~ 1 |Block). block is correspond with the number of treatment in each management. Management 1 only has one treatment, so block is 1；Management 2 has two treatments,block is 1 and 2；Management 3 has three treatments，block is 1,2,and 3

#load data
mydata1<- read_excel("all_matadata.xlsx",sheet = 1,col_names = T)

#composite predictor of microbial shannon
model1 <- lm(Belowground_biomass ~ Bac_shannon+Fun_shannon, mydata1)
coefs(model1, standardize = 'scale') #record as the bacterial,fungal,and protist significance to plant growth
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]
shannon <- beta_Bac * mydata1$Bac_shannon + beta_Fun * mydata1$Fun_shannon
mydata1$shannon <- shannon
summary(lm(Belowground_biomass ~ shannon,mydata1))
coefs(lm(Belowground_biomass ~ shannon,mydata1))

#composite predictor of microbial composition
model1 <- lm(Belowground_biomass ~ Bac_PCoA1+Fun_PCoA1, mydata1)
coefs(model1, standardize = 'scale')
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]

composition <- beta_Bac * mydata1$Bac_PCoA1 + beta_Fun * mydata1$Fun_PCoA1
mydata1$composition <- composition
summary(lm(Belowground_biomass ~ composition,mydata1))
coefs(lm(Belowground_biomass ~ composition,mydata1))


#multiple regression
microbe.list <- list(
  lme(Belowground_biomass~ Management + Variety + composition + shannon , random = list(~ 1 |Plot,~ 1 |Block), na.action = na.omit,
      data = mydata1),
  lme(composition~ Management+Variety, random = list(~ 1 |Plot,~ 1 |Block), na.action = na.omit,
      data = mydata1),
  lme(shannon~ Management+Variety, random = list(~ 1 |Plot,~ 1 |Block), na.action = na.omit,
      data = mydata1)
)


microbe.psem <- as.psem(microbe.list)
(new.summary <- summary(microbe.psem, .progressBar = F))
plot(microbe.psem,return = FALSE,alpha = 0.05,show = "std")



##Control management####

#load data
mydata2 <- mydata1 %>% filter(Management == "1")

#composite predictor of microbial shannon
model1 <- lm(Belowground_biomass ~ Bac_shannon+Fun_shannon, mydata2)
coefs(model1, standardize = 'scale')
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]

shannon <- beta_Bac * mydata2$Bac_shannon + beta_Fun * mydata2$Fun_shannon
mydata2$shannon <- shannon
summary(lm(Belowground_biomass ~ shannon,mydata2))
coefs(lm(Belowground_biomass ~ shannon,mydata2))

#composite predictor of microbial composition
model1 <- lm(Belowground_biomass ~ Bac_PCoA1+Fun_PCoA1, mydata2)
coefs(model1, standardize = 'scale')
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]

composition <- beta_Bac * mydata2$Bac_PCoA1 + beta_Fun * mydata2$Fun_PCoA1 
mydata2$composition <- composition
summary(lm(Belowground_biomass ~ composition,mydata2))
coefs(lm(Belowground_biomass ~ composition,mydata2))


#multiple regression
microbe.list <- list(
  lme(Belowground_biomass~Variety + composition + shannon ,random = ~ 1 | Plot, na.action = na.omit,
      data = mydata2),
  lme(composition~ Variety, random = ~ 1 | Plot, na.action = na.omit,
      data = mydata2),
  lme(shannon~ Variety, random = ~ 1 | Plot, na.action = na.omit,
      data = mydata2)
)


microbe.psem <- as.psem(microbe.list)
(new.summary <- summary(microbe.psem, .progressBar = F))
plot(microbe.psem,return = FALSE,alpha = 0.05,show = "std")


##Biologiacl management####

#load data
mydata3 <- mydata1 %>% filter(Management == "2")

#composite predictor of microbial shannon
model1 <- lm(Belowground_biomass ~ Bac_shannon+Fun_shannon, mydata3)
coefs(model1, standardize = 'scale')
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]

shannon <- beta_Bac * mydata3$Bac_shannon + beta_Fun * mydata3$Fun_shannon 
mydata3$shannon <- shannon
summary(lm(Belowground_biomass ~ shannon,mydata3))
coefs(lm(Belowground_biomass ~ shannon,mydata3))

#composite predictor of microbial composition
model1 <- lm(Belowground_biomass ~ Bac_PCoA1+Fun_PCoA1, mydata3)
coefs(model1, standardize = 'scale')
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]

composition <- beta_Bac * mydata3$Bac_PCoA1 + beta_Fun * mydata3$Fun_PCoA1 
mydata3$composition <- composition
summary(lm(Belowground_biomass ~ composition,mydata3))
coefs(lm(Belowground_biomass ~ composition,mydata3))


#multiple regression
microbe.list <- list(
  lme(Belowground_biomass~ Treatment + Variety + composition + shannon , random = ~ 1 | Plot, na.action = na.omit,
      data = mydata3),
  lme(composition~ Treatment + Variety, random = ~ 1 | Plot, na.action = na.omit,
      data = mydata3),
  lme(shannon~ Treatment + Variety, random = ~ 1 | Plot, na.action = na.omit,
      data = mydata3)
)


microbe.psem <- as.psem(microbe.list)
(new.summary <- summary(microbe.psem, .progressBar = F))
plot(microbe.psem,return = FALSE,alpha = 0.05,show = "std")


##Chemicial management####

#load data
mydata4 <- mydata1 %>% filter(Management == "3")

#composite predictor of microbial shannon
model1 <- lm(Belowground_biomass ~ Bac_shannon+Fun_shannon, mydata4)
coefs(model1, standardize = 'scale')
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]

shannon <- beta_Bac * mydata4$Bac_shannon + beta_Fun * mydata4$Fun_shannon
summary(lm(Belowground_biomass ~ shannon,mydata4))
coefs(lm(Belowground_biomass ~ shannon,mydata4))

#composite predictor of microbial composition
model1 <- lm(Belowground_biomass ~ Bac_PCoA1+Fun_PCoA1, mydata4)
coefs(model1, standardize = 'scale')
beta_Bac <-  summary(model1)$coefficients[2, 1]
beta_Fun <- summary(model1)$coefficients[3, 1]

composition <- beta_Bac * mydata4$Bac_PCoA1 + beta_Fun * mydata4$Fun_PCoA1 
mydata4$composition <- composition
summary(lm(Belowground_biomass ~ composition,mydata4))
coefs(lm(Belowground_biomass ~ composition,mydata4))


#multiple regression
microbe.list <- list(
  lme(Belowground_biomass~ Treatment + Variety + composition + shannon , random = ~ 1 | Plot, na.action = na.omit,
      data = mydata4),
  lme(composition~ Treatment + Variety, random = ~ 1 | Plot, na.action = na.omit,
      data = mydata4),
  lme(shannon~ Treatment + Variety, random = ~ 1 | Plot, na.action = na.omit,
      data = mydata4)
)


microbe.psem <- as.psem(microbe.list)
(new.summary <- summary(microbe.psem, .progressBar = F))
plot(microbe.psem,return = FALSE,alpha = 0.05,show = "std")



