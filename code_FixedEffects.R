################################################################################
################ structural effects of cluster policies ########################

#### set up
rm(list = ls())
getwd()
options(digits = 5)

#install.packages("")

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
# change group name
grps <- grps %>% mutate(groupe = case_when(
        groupe == "3" ~ "G1",
        groupe == "2" ~ "G2",
        groupe == "4" ~ "G3",
        groupe == "1" ~ "G4")
        )
controls <- controls %>% left_join(grps, by = "dep")
colnames(controls)[10] <- "group"
# save grps for plot
write.csv(grps, "./carto/grps.csv")
rm(deps, grps)


## final df
df1 <- left_join(outcomes, tr.int[, -1], by = c("nccenr", "period"))
df1 <- left_join(df1, controls, by = c("nccenr", "period"))
glimpse(df1)


df1 <- df1[, c(1, 20, 23, 21, 22, 24:30, 2:19)]
colnames(df1)[2] <- "region"
colnames(df1)[5] <- "treatment_int"
df1$dep <- as.factor(df1$dep)
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

cor(df1[, 6:10])
cor(df1[, 12:29])



library(plm)
library(AER)
library(tseries)
library(lmtest)

################################################################################
######################################################### Descriptive statistics 
################################################################################
glimpse(df1)

## data for descriptive statistics  
df1.desc <- df1 %>% select(region, period, group, # independent variables
                           nb_nodes, nb_patents, # network composition
                           gdp, dird, sub_region, sub_nat, sub_cee, treatment_int, 
                           net_density, fragmentation_index, share_net_main_comp, # dependent variables
                           CC_ratio, PL_ratio, net_hierarchy, net_assortativity,
                           share_local_nodes, share_regional_nodes, share_national_nodes)

df1.desc$CC_ratio <- log(df1.desc$CC_ratio)
glimpse(df1.desc)

## global per period
mat1 <- psych::describeBy(df1.desc$nb_patents, df1.desc$period, digits = 5, mat = T)[, -c(1, 3)]
colnames(mat1)[1] <- "period"
row.names(mat1) <- NULL
mat1

## group per period
mat2 <- psych::describeBy(df1.desc$nb_patents, list(df1.desc$period, df1.desc$group),
                          digits = 5, mat = T)[, -c(1, 4)]
colnames(mat2)[1:2] <- c("period", "group")
row.names(mat2) <- NULL
mat2


########################################################################## TEST
model0 <- plm(log(nb_patents) ~ gdp + dird + sub_region + sub_nat + sub_cee +
                      net_density + fragmentation_index + share_net_main_comp +
                      CC_actual + PL_actual + share_local_nodes + share_regional_nodes,
              data = df1, index = c("region", "period"), model = "within",
              effect = "twoways")
summary(model0)
coeftest(model0, vcov = vcovHC, type = "HC1")



################################################################################
############################################ Estimations 1 : fixed effects model 
################################################################################

### 1) Network embeddedness : density
## estimating the fixed effects regression with plm()
model1 <- plm(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : fragmentation index
model1 <- plm(fragmentation_index ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : share of the network’s main component
model1 <- plm(share_net_main_comp ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

################################################################################

### 2) Network efficiency : clustering coefficient (ratio)
model1 <- plm(log(CC_ratio) ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 2) Network efficiency : average path length (ratio)
model1 <- plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

################################################################################

### 3) Network resilience : network hierarchy
model1 <- plm(net_hierarchy ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 3) Network resilience : network assortativity
model1 <- plm(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

################################################################################

### 4) Network geographical anchoring : share of local nodes
model1 <- plm(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of regional nodes
model1 <- plm(share_regional_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of national nodes
model1 <- plm(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")



################################################################################
######################################### Estimations 2: treatment heterogeneity 
################################################################################

### 1) Network embeddedness : density
## estimating the fixed effects regression with plm()
model2 <- plm(net_density ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : fragmentation index
model2 <- plm(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : share of the network’s main component
model2 <- plm(share_net_main_comp ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

################################################################################

### 2) Network efficiency : clustering coefficient (ratio)
model2 <- plm(log(CC_ratio) ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 2) Network efficiency : average path length (ratio)
model2 <- plm(PL_ratio ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

################################################################################

### 3) Network resilience : network hierarchy
model2 <- plm(net_hierarchy ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 3) Network resilience : network assortativity
model2 <- plm(net_assortativity ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

################################################################################

### 4) Network geographical anchoring : share of local nodes
model2 <- plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of regional nodes
model2 <- plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of national nodes
model2 <- plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "individual")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

################################################################################
################### Estimations 3 : Time fixed effects + treatment heterogeneity 
################################################################################

### 1) Network embeddedness : density
## estimating the fixed effects regression with plm()
model3 <- plm(net_density ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : fragmentation index
model3 <- plm(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : share of the network’s main component
model3 <- plm(share_net_main_comp ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

################################################################################

### 2) Network efficiency : clustering coefficient (ratio)
model3 <- plm(log(CC_ratio) ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

### 2) Network efficiency : average path length (ratio)
model3 <- plm(PL_ratio ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

################################################################################

### 3) Network resilience : network hierarchy
model3 <- plm(net_hierarchy ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

### 3) Network resilience : network assortativity
model3 <- plm(net_assortativity ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

################################################################################

### 4) Network geographical anchoring : share of local nodes
model3 <- plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of regional nodes
model3 <- plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of national nodes
model3 <- plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model3)
coeftest(model3, vcov = vcovHC, type = "HC1")



################################################################################
###################################### Estimations 4: Spatial Durbin model (SDM) 
################################################################################
library(splm)
library(spdep)

glimpse(df1)
df2 <- select(df1, -net_name, -region)
glimpse(df2)

################################################### STOP : export df2 for MATLAB
str(df2)
df2.matlab <- df2
df2.matlab <- df2.matlab[, c(2, 1, 3:27)]
str(df2.matlab)
df2.matlab <- df2.matlab %>% arrange(period, dep)

# KSI data
ksi <- read.csv("T:/These_GATE/Paper_2/KSI/KSI.csv")
glimpse(ksi)
ksi$dep <- ifelse(nchar(ksi$dep) == 1, paste0("0", ksi$dep), ksi$dep)
ksi$dep  <- as.factor(ksi$dep)
# join df2.matlab and KSI data
glimpse(df2.matlab)
glimpse(ksi)
df2.matlab <- left_join(df2.matlab, ksi, by = "dep")


# export en xls
library(xlsx)
write.xlsx(df2.matlab, "T:/These_GATE/Paper_2/Estimation1/df_with_group_ksi.xls",
           row.names = F)

################################################################################

sp.mat <- read.table("../spatial_matrix.txt")
#sp.mat <- read.table("../net_matrix2.txt")

row.names(sp.mat)
colnames(sp.mat) <- substr(colnames(sp.mat), 2, 3)
sp.mat <- as.matrix(sp.mat)
rowSums(sp.mat)
sp.matl <- mat2listw(sp.mat, style = "W")


# SDM (for geo spatial matrix)
glimpse(df2)
df3 <- pdata.frame(df2, index = c("dep", "period"))
df3$treatment_int_SL <- slag(df3$treatment_int, sp.matl)
df3$gdp_SL <- slag(df3$gdp, sp.matl)
df3$dird_SL <- slag(df3$dird, sp.matl)
df3$sub_region_SL <- slag(df3$sub_region, sp.matl)
df3$sub_nat_SL <- slag(df3$sub_nat, sp.matl)
df3$sub_cee_SL <- slag(df3$sub_cee, sp.matl)
glimpse(df3)
colnames(df3)


# SDM
summary(
        model4 <- spml(net_density ~ treatment_int : group + gdp + dird +
                                   sub_region + sub_nat + sub_cee + treatment_int_SL : group +
                                   gdp_SL + dird_SL + sub_region_SL + sub_nat_SL + 
                                   sub_cee_SL, data = df3, index = c("dep", "period"),
                           model = "within", effect = "twoways", lag = TRUE,
                           listw = sp.matl, spatial.error = "none", LeeYu = TRUE, Hess = FALSE)
        )
summary(model4)$rsqr
effects.splm(model4)


## calculate impact measures
time <- length(unique(df3$period))
set.seed(1234)
imps <- impacts(model4, listw = sp.matl, time = time)
impacts(model4)
summary(imps, zstats = TRUE, short = TRUE)

# test
model4$spat.coef
iIrW  <- invIrW(sp.matl, model4$spat.coef)
model4$coefficients
model4$coefficients[2]
model4$coefficients[7]

S_gdp <- iIrW %*% ((model4$coefficients[2] * diag(94)) - # diag(nrow(df3))
                           (model4$coefficients[7] * listw2mat(sp.matl)))


### Bivand's answer 2
Time <- length(unique(df3$period))
N <- length(unique(df3$dep))

library(Matrix)
s.lws <- kronecker(Diagonal(Time) , listw2dgCMatrix(sp.matl))
IrW <- Diagonal(N * Time) - model4$spat.coef * s.lws
IrWi <- solve(IrW)
S_gdp <- IrWi * (Diagonal(N * Time) * model4$coefficients[2] +
                         s.lws * model4$coefficients[7])

# diret impacts
mean(diag(S_gdp))
sd(diag(S_gdp))
t.test(diag(S_gdp), mu = 0)$p.value
LaplacesDemon::p.interval(diag(S_gdp), prob = 0.95)

# total impacts
mean(rowSums(S_gdp)) # sum(S_gdp)/376
sd(rowSums(S_gdp))
t.test(rowSums(S_gdp), mu = 0)$p.value
LaplacesDemon::p.interval(rowSums(S_gdp), prob = 0.95)

# indirect impacts
mean(rowSums(S_gdp)) - mean(diag(S_gdp)) 
        







### Bivand's answer 1

#S_gdp <- iIrW %*% (model4$coefficients[2] * diag(94)) # gives wrong impact got from impacts()

# Then you can get the direct and total impacts in the usual way: 
dir_gdp <-  mean(diag(S_gdp)) # or sum(diag(S_gdp))/94
tot_gdp <- mean(rowSums(S_gdp)) # or sum(S_gdp)/94
indir_gdp <- sum(S_gdp[lower.tri(S_gdp) | upper.tri(S_gdp)]) / 94 # or tot_gdp - dir_gdp 


sd(diag(S_gdp))
LaplacesDemon::p.interval(diag(S_gdp), prob = 0.95)

sd(rowSums(S_gdp))
LaplacesDemon::p.interval(c(S_gdp), prob = 0.95)


sd(c(S_gdp[lower.tri(S_gdp) | upper.tri(S_gdp)]))
LaplacesDemon::p.interval(c(S_gdp[lower.tri(S_gdp) | upper.tri(S_gdp)]), prob = 0.95)


t.test(diag(S_gdp), mu = 0)$p.value
t.test(rowSums(S_gdp), mu = 0)$p.value



impacts
methods(impacts)
getS3method("impacts", "sarlm") 

library(spatialreg)
spatialreg::intImpacts
spatialreg::lagDistrImpacts
spatialreg::lagImpacts
spatialreg::lagImpacts_e

MASS::mvrnorm

matrix(c(10,3,3,2),2,2)

#

#### Seline's suggestion
# SLX
summary(
        model3.SLX <- plm(net_density ~ treatment_int : group + gdp + dird + 
                                  sub_region + sub_nat + sub_cee + treatment_int_SL : group +
                                  gdp_SL + dird_SL + sub_region_SL + sub_nat_SL + sub_cee_SL,
                          quiet = FALSE, data = df3, index = c("dep", "period"),
                          model = "within", effect = "twoways")
        )


# SAR
summary(
        model3.SAR <- spml(net_density ~ treatment_int : group + gdp + dird +
                                   sub_region + sub_nat + sub_cee,
                           data = df3, index = c("dep", "period"),
                           model = "within", effect = "twoways", lag = TRUE,
                           listw = sp.matl, spatial.error = "none",
                           LeeYu = TRUE, Hess = FALSE)
)
model3.SAR$type

summary(model3.SAR)$rsqr
model3.SAR$logLik


### Models comparison: SDM vs SAR and SDM vs SLX
## SDM vs SAR 
AICsplm(model4, criterion = "AIC")
AICsplm(model3.SAR, criterion = "AIC")
AICsplm(model4, criterion = "BIC")
AICsplm(model3.SAR, criterion = "BIC")
# SDM wins SAR

## SDM vs SLX 
AICsplm(model4, criterion = "AIC")
AICplm(model3.SLX, criterion = "AIC")
AICsplm(model4, criterion = "BIC")
AICplm(model3.SLX, criterion = "BIC")
#



## SDM (GM estimation)
summary(
        model4 <- spgm(net_density ~ treatment_int : group + gdp + dird + sub_region +
                               sub_nat + sub_cee +
                               treatment_int_SL : group + gdp_SL + dird_SL + sub_region_SL +
                               sub_nat_SL + sub_cee_SL,
             data = df3, index = c("dep", "period"), model = "within",
             listw = sp.matl, moments = "fullweights", lag = TRUE,
             spatial.error = FALSE)
)
summary(model4)$rsqr

# https://cran.r-project.org/web/packages/export/export.pdf


