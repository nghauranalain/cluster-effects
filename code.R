################################################################################
################ structural effects of cluster policies ########################

#### set up
rm(list = ls())
getwd()
options(digits = 5)

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


## outcomes and controls
# outcomes
outcomes <-  read.csv("T:/These_GATE/Paper_2/Outcomes/res2.csv",
                      stringsAsFactors = F)
glimpse(outcomes)
outcomes$nccenr <- substr(outcomes$net_name, 1, nchar(outcomes$net_name)-2)
outcomes$period <- substr(outcomes$net_name, nchar(outcomes$net_name), nchar(outcomes$net_name))
# controls
controls <- read.csv("T:/These_GATE/Paper_2/Controls/Data Eurolio/data_eurolio1.csv",
                     header = T, colClasses = c(rep("character", 2),
                                                rep("numeric", 6)))
glimpse(controls)
colnames(controls)[1] <- "dep"
glimpse(deps)
controls <- controls %>% left_join(deps[, c("dep", "nccenr")], by = "dep")

# add group variable
grps <- read.csv("T:/These_GATE/Paper_2/Controls/Data Eurolio/Data/donnees_departements.csv",
                 header = T, colClasses = "character")
glimpse(grps)
grps <- select(grps, dep, groupe)
grps$dep <- ifelse(nchar(grps$dep) == 1, paste0("0", grps$dep), grps$dep)
controls <- controls %>% left_join(grps, by = "dep") %>% select(-dep)
colnames(controls)[9] <- "group"
rm(deps, grps)


## final df
df1 <- left_join(outcomes, tr.int[, -1], by = c("nccenr", "period"))
df1 <- left_join(df1, controls, by = c("nccenr", "period"))
glimpse(df1)


df1 <- df1[, c(1, 20:22, 23:29, 2:19)]
colnames(df1)[2] <- "region"
colnames(df1)[4] <- "treatment_int"
df1$region <- as.factor(df1$region)
df1$group <- as.factor(df1$group)
df1$period <- as.factor(df1$period)
df1$treatment_int <- as.numeric(df1$treatment_int)
df1 <- select(df1, -tot_sub) # not useful

glimpse(df1)
summary(df1)

table(df1$sub_region == 0)
table(df1$sub_nat == 0)
table(df1$sub_cee == 0)

## log values
# add 1 euro (0.001) to sub_region, sub_nat, sub_cee and 
df1$gdp <- log(df1$gdp)
df1$dird <- log(df1$dird)
df1$sub_region <- log(df1$sub_region + 0.001)
df1$sub_nat <- log(df1$sub_nat + 0.001)
df1$sub_cee <- log(df1$sub_cee + 0.001)

cor(df1[, 5:10])
cor(df1[, 11:28])



library(plm)
library(AER)
library(tseries)
library(lmtest)

################################################################################
################################################################## Estimations 1 
################################################################################

#cntrls <- c("gdp", "dird", "sub_region", "sub_nat", "sub_cee")
#as.formula(paste("net_density ~ ", paste(cntrls, collapse = "+")))


### 1) Network embeddedness : density
## fixed / random effects regression: Hausman test
phtest(plm(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects. Here, no need to use time-fixed effects
## testing for unit roots/stationarity
adf.test(df1$net_density) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano")) 
# assessing multicollinearity 
car::vif(plm(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
             data = df1, index = c("region", "period"), model = "pooling"))




### 1) Network embeddedness : fragmentation index
## fixed / random effects regression: Hausman test
phtest(plm(fragmentation_index ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(fragmentation_index ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(fragmentation_index ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(fragmentation_index ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model1 <- plm(fragmentation_index ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
## testing for unit roots/stationarity
adf.test(df1$fragmentation_index) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(fragmentation_index ~ treatment_int, data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano"))




### 1) Network embeddedness : share of the network’s main component
## fixed / random effects regression: Hausman test
phtest(plm(share_net_main_comp ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_net_main_comp ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(share_net_main_comp ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(share_net_main_comp ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects. Here, no need to use time-fixed effects
## testing for unit roots/stationarity
adf.test(df1$share_net_main_comp) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(share_net_main_comp ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano")) 




### 2) Network efficiency : clustering coefficient (ratio)
## fixed / random effects regression: Hausman test
phtest(plm(CC_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(CC_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(CC_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "random")
summary(model1)
## testing for random effects: Breusch-Pagan Lagrange multiplier (LM)
plmtest(plm(CC_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
            data = df1, index = c("region", "period"),
            model = "pooling"), type = c("bp")) # If p-value < 0.05 then random effects model is not
# appropriate (compare to a simple OLS regression)
## testing for unit roots/stationarity
adf.test(df1$CC_ratio) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(CC_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano"))




### 2) Network efficiency : average path length (ratio)
## fixed / random effects regression: Hausman test
phtest(plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model1 <- plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "random")
summary(model1)
## testing for random effects: Breusch-Pagan Lagrange multiplier (LM)
plmtest(plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
            data = df1, index = c("region", "period"),
            model = "pooling"), type = c("bp")) # If p-value < 0.05 then random effects model is not
# appropriate (compare to a simple OLS regression)

model1 <- plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"),
              model = "pooling") # pooled OLS

summary(plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
            data = df1, index = c("region", "period"),
            model = "within"))

summary(model1)
## testing for unit roots/stationarity
adf.test(df1$PL_ratio) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano"))




### 3) Network resilience : network hierarchy
## fixed / random effects regression: Hausman test
phtest(plm(net_hierarchy ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(net_hierarchy ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(net_hierarchy ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(net_hierarchy ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects. Here, no need to use time-fixed effects
## testing for unit roots/stationarity
adf.test(df1$net_hierarchy) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(net_hierarchy ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano")) 




### 3) Network resilience : network assortativity
## fixed / random effects regression: Hausman test
phtest(plm(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model1 <- plm(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
## testing for unit roots/stationarity
adf.test(df1$net_assortativity) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano")) 




### 4) Network geographical anchoring : share of local nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model1 <- plm(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
## testing for unit roots/stationarity
adf.test(df1$share_local_nodes) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano")) 




### 4) Network geographical anchoring : share of regional nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_regional_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_regional_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the random effects regression with plm()
model1 <- plm(share_regional_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(share_regional_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects. Here, no need to use time-fixed effects
## testing for unit roots/stationarity
adf.test(df1$share_regional_nodes) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(share_regional_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano")) 




### 4) Network geographical anchoring : share of national nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model1 <- plm(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model1)
## testing for time-fixed effects
pFtest(plm(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model1) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model1 <- plm(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
## testing for unit roots/stationarity
adf.test(df1$share_national_nodes) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model1) # original coefficients
coeftest(model1, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model1, vcovHC(model1, method = "arellano")) 





################################################################################
################################################################## Estimations 2 
################################################################################
glimpse(df1)
cor(df1[, 5:9])



### 1) Network embeddedness : density
## fixed / random effects regression: Hausman test
phtest(plm(net_density ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(net_density ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model2 <- plm(net_density ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "random")
summary(model2)
## testing for random effects: Breusch-Pagan Lagrange multiplier (LM)
plmtest(plm(net_density ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
            data = df1, index = c("region", "period"),
            model = "pooling"), type = c("bp")) # If p-value < 0.05 then random effects model is not
# appropriate (compare to a simple OLS regression)
model2 <- plm(net_density ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"),
              model = "pooling") # pooled OLS
summary(model2)
## testing for unit roots/stationarity
adf.test(df1$net_density) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(net_density ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano"))





### 1) Network embeddedness : fragmentation index
## fixed / random effects regression: Hausman test
phtest(plm(fragmentation_index ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(fragmentation_index ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the random effects regression with plm()
model2 <- plm(fragmentation_index ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model2)
## testing for time-fixed effects
pFtest(plm(fragmentation_index ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model2) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model2 <- plm(fragmentation_index ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
## testing for unit roots/stationarity
adf.test(df1$fragmentation_index) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(fragmentation_index ~ treatment_int * group, data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for random effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC(model2, type = "HC1")) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 




### 1) Network embeddedness : share of the network’s main component
## fixed / random effects regression: Hausman test
phtest(plm(share_net_main_comp ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_net_main_comp ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value > 0.05, then use random effects
## estimating the fixed effects regression with plm()
model2 <- plm(share_net_main_comp ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "random")
summary(model2)
## testing for random effects: Breusch-Pagan Lagrange multiplier (LM)
plmtest(plm(share_net_main_comp ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
            data = df1, index = c("region", "period"),
            model = "pooling"), type = c("bp")) # If p-value < 0.05 then random effects model is not
# appropriate (compare to a simple OLS regression)
model2 <- plm(share_net_main_comp ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"),
              model = "pooling") # pooled OLS
summary(model2)
## testing for unit roots/stationarity
adf.test(df1$share_net_main_comp) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(share_net_main_comp ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 




### 2) Network efficiency : clustering coefficient (ratio)
## fixed / random effects regression: Hausman test
phtest(plm(CC_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(CC_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model2 <- plm(CC_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model2)
## testing for time-fixed effects
pFtest(plm(CC_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model2) # If p-value < 0.05 then use time-fixed effects
## testing for unit roots/stationarity
adf.test(df1$CC_ratio) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(CC_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 



### 2) Network efficiency : average path length (ratio)
## fixed / random effects regression: Hausman test
phtest(plm(PL_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(PL_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model2 <- plm(PL_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model2)
## testing for time-fixed effects
pFtest(plm(PL_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model2) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model2 <- plm(PL_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
## testing for unit roots/stationarity
adf.test(df1$PL_ratio) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(PL_ratio ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 




### 3) Network resilience : network hierarchy
## fixed / random effects regression: Hausman test
phtest(plm(net_hierarchy ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(net_hierarchy ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model2 <- plm(net_hierarchy ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model2)
## testing for time-fixed effects
pFtest(plm(net_hierarchy ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model2) # If p-value < 0.05 then use time-fixed effects. Here, no need to use time-fixed effects
## testing for unit roots/stationarity
adf.test(df1$net_hierarchy) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(net_hierarchy ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 




### 3) Network resilience : network assortativity
## fixed / random effects regression: Hausman test
phtest(plm(net_assortativity ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(net_assortativity ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model2 <- plm(net_assortativity ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model2)
## testing for time-fixed effects
pFtest(plm(net_assortativity ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model2) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model2 <- plm(net_assortativity ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
## testing for unit roots/stationarity
adf.test(df1$net_assortativity) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(net_assortativity ~ treatment_int * group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 




### 4) Network geographical anchoring : share of local nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model2 <- plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model2)
## testing for random effects: Breusch-Pagan Lagrange multiplier (LM)
plmtest(plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
            data = df1, index = c("region", "period"),
            model = "pooling"), type = c("bp")) # If p-value < 0.05 then random effects model is not
# appropriate (compare to a simple OLS regression)
model2 <- plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"),
              model = "pooling") # pooled OLS
summary(model2)
## testing for unit roots/stationarity
adf.test(df1$share_local_nodes) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 



### 4) Network geographical anchoring : share of regional nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value > 0.05, then use random effects
## estimating the random effects regression with plm()
model2 <- plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within")
summary(model2)
## testing for random effects: Breusch-Pagan Lagrange multiplier (LM)
plmtest(plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
            data = df1, index = c("region", "period"),
            model = "pooling"), type = c("bp")) # If p-value < 0.05 then random effects model is not
# appropriate (compare to a simple OLS regression)
model2 <- plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"),
              model = "pooling") # pooled OLS
summary(model2)
## testing for heteroskedasticity 
bptest(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov. = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 
# assessing multicollinearity 
car::vif(plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
             data = df1, index = c("region", "period"),
             model = "pooling"))




### 4) Network geographical anchoring : share of national nodes
## fixed / random effects regression: Hausman test
phtest(plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "within"),
       plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
           data = df1, index = c("region", "period"), model = "random")) # p-value < 0.05, then use fixed effects
## estimating the fixed effects regression with plm()
model2 <- plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "random")
summary(model2)
## testing for time-fixed effects
pFtest(plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee + period,
           data = df1, index = c("region", "period"), model = "within"),
       model2) # If p-value < 0.05 then use time-fixed effects
## time fixed effects
model2 <- plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
## testing for unit roots/stationarity
adf.test(df1$share_national_nodes) # If p-value < 0.05 then no unit roots present
## testing for heteroskedasticity 
bptest(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
       data = df1, studentize = FALSE) # if p-value < 0.05, presence of heteroskedasticity 
# correction (for fixed effects)
coeftest(model2) # original coefficients
coeftest(model2, vcov = vcovHC) # robust standard errors (heteroskedasticity consistent coefficients)
coeftest(model2, vcovHC(model2, method = "arellano")) 



cor(df1[, 5:9])
glimpse(df1)



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

























