################################################################################
################ structural effects of cluster policies ########################

#### set up
rm(list = ls())
getwd()

#### packages
library(tidyverse)
library(data.table)

#### load data
## treatment intensity
tr.int <- read.csv("T:/These_GATE/Paper_2/Treatment/treatment_intensity.csv",
                   colClasses = "character")
glimpse(tr.int)
tr.int <- tr.int %>% select(dep, period, treatment.int.period) %>%
        mutate(dep = as.character(dep), period = as.character(period)) %>%
        distinct()
# data with dep code
deps <- read.delim("T:/These_GATE/Paper_2/Other data/depts2013.txt",
                   header = T, stringsAsFactors = F)
glimpse(deps)
colnames(deps) <- tolower(colnames(deps))
tr.int <- left_join(tr.int, deps[, c("dep", "nccenr")], by = "dep")
rm(deps)

## outcomes
outcomes <-  read.csv("T:/These_GATE/Paper_2/Outcomes/res2.csv",
                      stringsAsFactors = F)
glimpse(outcomes)
outcomes$nccenr <- substr(outcomes$net_name, 1, nchar(outcomes$net_name)-2)
outcomes$period <- substr(outcomes$net_name, nchar(outcomes$net_name), nchar(outcomes$net_name))

## final df
df1 <- left_join(outcomes, tr.int[, -1], by = c("nccenr", "period"))
glimpse(df1)
df1 <- df1[, c(1, 20:22, 2:19)]
colnames(df1)[2] <- "region"
colnames(df1)[4] <- "treatment_int"
df1$region <- as.factor(df1$region)
df1$period <- as.factor(df1$period)
df1$treatment_int <- as.numeric(df1$treatment_int)


#### estimations
glimpse(df1)
summary(df1)
cor(df1[, 5:22])

library(plm)
library(AER)
library(tseries)
library(lmtest)

### 1) Network embeddedness : density
## fixed / random effects regression: Hausman test
phtest(plm(net_density ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(net_density ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(net_density ~ treatment_int, data = df1, index = c("region", "period"),
              model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(net_density ~ treatment_int + period, data = df1, index = c("region", "period"),
            model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects. Here, no need to use time-fixed effects
## testing for unit roots/stationarity
adf.test(df1$net_density) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(net_density ~ treatment_int, data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov. = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)


### 1) Network embeddedness : fragmentation index
## fixed / random effects regression: Hausman test
phtest(plm(fragmentation_index ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(fragmentation_index ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(fragmentation_index ~ treatment_int, data = df1, index = c("region", "period"),
              model = "random")
summary(model1)
## testing for random effects: Breusch-Pagan Lagrange multiplier (LM)
plmtest(plm(fragmentation_index ~ treatment_int, data = df1, index = c("region", "period"),
            model = "pooling"), type = c("bp")) # If p-value < 0.05 then random effects model is not
                                                # appropriate (compare to a simple OLS regression)
## testing for unit roots/stationarity
adf.test(df1$fragmentation_index) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(fragmentation_index ~ treatment_int, data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for random effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov. = vcovHC(model1, type = "HC1")) # robust standard errors (heteroskedasticity consistent coefficients)

### 1) Network embeddedness : share of the network’s main component
## fixed / random effects regression: Hausman test
phtest(plm(share_net_main_comp ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(share_net_main_comp ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(share_net_main_comp ~ treatment_int, data = df1, index = c("region", "period"),
              model = "within")
summary(model1)


### 2) Network efficiency : clustering coefficient (ratio)
## fixed / random effects regression: Hausman test
phtest(plm(CC_ratio ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(CC_ratio ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(share_net_main_comp ~ treatment_int, data = df1, index = c("region", "period"),
              model = "random")
summary(model1)


### 2) Network efficiency : average path length (ratio)
## fixed / random effects regression: Hausman test
phtest(plm(PL_ratio ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(PL_ratio ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(share_net_main_comp ~ treatment_int, data = df1, index = c("region", "period"),
              model = "random")
summary(model1)

### 3) Network resilience : network hierarchy
## fixed / random effects regression: Hausman test
phtest(plm(net_hierarchy ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(net_hierarchy ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(net_hierarchy ~ treatment_int, data = df1, index = c("region", "period"),
              model = "random")
summary(model1)

### 3) Network resilience : network assortativity
## fixed / random effects regression: Hausman test
phtest(plm(net_assortativity ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(net_assortativity ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(net_assortativity ~ treatment_int, data = df1, index = c("region", "period"),
              model = "random")
summary(model1)

### 4) Network geographical anchoring : share of local nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_local_nodes ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(share_local_nodes ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(share_local_nodes ~ treatment_int, data = df1, index = c("region", "period"),
              model = "within")
summary(model1)

### 4) Network geographical anchoring : share of regional nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_regional_nodes ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(share_regional_nodes ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(share_regional_nodes ~ treatment_int, data = df1, index = c("region", "period"),
              model = "random")
summary(model1)

### 4) Network geographical anchoring : share of national nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_national_nodes ~ treatment_int, data = df1, index = c("region", "period"),
           model = "within"),
       plm(share_national_nodes ~ treatment_int, data = df1, index = c("region", "period"),
           model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(share_national_nodes ~ treatment_int, data = df1, index = c("region", "period"),
              model = "within")
summary(model1)









################################################################### STOP
################################################################### TEST

# estimating the fixed effects regression with lm()
model2 <- lm(fragmentation_index ~ treatment_int + region - 1, data = df1)
summary(model2)


################################################################### DiD estimation
# creating a dummy variable to indicate the time when the treatment started
df1$time <- ifelse(as.character(df1$period) %in% c("1", "2"), 0, 1)
# creating a dummy variable to identify the group exposed to the treatment. 
head(df1$treatment_int, 12) # treatment is continuous
# creating an interaction between time and treated. We will call this interaction ‘did’.
df1$did <- df1$time * df1$treatment_int

# estimating the DID estimator (with lm)
model2 <- lm(fragmentation_index ~ treatment_int * time, data = df1)
summary(model2)

summary(lm(fragmentation_index ~ treatment_int + time + did, data = df1))



head(df1$did, 12)
head(df1$time + df1$treatment_int, 12)
cor(df1$did, df1$treatment_int)

# estimating the DID estimator (with plm)
model3 <- plm(fragmentation_index ~ did + time, data = df1,
              index = c("region", "period"), model = "within")
summary(model3)
coeftest(model3, vcov. = vcovHC, type = "HC1")

























