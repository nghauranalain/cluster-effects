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
controls <- controls %>% left_join(grps, by = "dep")
colnames(controls)[10] <- "group"
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


########################################################################## TEST
model0 <- plm(log(nb_patents) ~ gdp + dird + sub_region + sub_nat + sub_cee +
                      net_density + fragmentation_index + share_net_main_comp +
                      CC_actual + PL_actual + share_local_nodes + share_regional_nodes,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model0)
coeftest(model0, vcov = vcovHC, type = "HC1")

################################################################################
################################################################## Estimations 1 
################################################################################

### 1) Network embeddedness : density
## estimating the fixed effects regression with plm()
model1 <- plm(net_density ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : fragmentation index
model1 <- plm(fragmentation_index ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : share of the network’s main component
model1 <- plm(share_net_main_comp ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

################################################################################

### 2) Network efficiency : clustering coefficient (ratio)
model1 <- plm(CC_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 2) Network efficiency : average path length (ratio)
model1 <- plm(PL_ratio ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

################################################################################

### 3) Network resilience : network hierarchy
model1 <- plm(net_hierarchy ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 3) Network resilience : network assortativity
model1 <- plm(net_assortativity ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

################################################################################

### 4) Network geographical anchoring : share of local nodes
model1 <- plm(share_local_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of regional nodes
model1 <- plm(share_regional_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of national nodes
model1 <- plm(share_national_nodes ~ treatment_int + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model1)
coeftest(model1, vcov = vcovHC, type = "HC1")



################################################################################
################################################################## Estimations 2 
################################################################################

### 1) Network embeddedness : density
## estimating the fixed effects regression with plm()
model2 <- plm(net_density ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 1) Network embeddedness : fragmentation index
model2 <- plm(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")
coeftest(model2, vcov = vcovHC(model2, type = "HC1", cluster = "group"))

### 1) Network embeddedness : share of the network’s main component
model2 <- plm(share_net_main_comp ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

################################################################################

### 2) Network efficiency : clustering coefficient (ratio)
model2 <- plm(CC_ratio ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 2) Network efficiency : average path length (ratio)
model2 <- plm(PL_ratio ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

################################################################################

### 3) Network resilience : network hierarchy
model2 <- plm(net_hierarchy ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 3) Network resilience : network assortativity
model2 <- plm(net_assortativity ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

################################################################################

### 4) Network geographical anchoring : share of local nodes
model2 <- plm(share_local_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of regional nodes
model2 <- plm(share_regional_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")

### 4) Network geographical anchoring : share of national nodes
model2 <- plm(share_national_nodes ~ treatment_int : group + gdp + dird + sub_region + sub_nat + sub_cee,
              data = df1, index = c("region", "period"), model = "within", effect = "twoways")
summary(model2)
coeftest(model2, vcov = vcovHC, type = "HC1")



################################################################################
################################################### Test spatial autocorrelation 
################################################################################
library(splm)
library(spdep)

#data("Produc", package = "Ecdat")
#data("usaww")


glimpse(df1)
df2 <- select(df1, -net_name, -region)
glimpse(df2)

sp.mat <- read.table("../spatial_matrix.txt")
#sp.mat <- read.table("../net_matrix2.txt")

row.names(sp.mat)
colnames(sp.mat) <- substr(colnames(sp.mat), 2, 3)
sp.mat <- as.matrix(sp.mat)
rowSums(sp.mat)
sp.matl <- mat2listw(sp.mat, style = "W")


## tests
# RLMlag robust version
slmtest(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region +
                sub_nat + sub_cee,
        data = df2, listw = sp.matl, model = "within", effect = "twoways",
        test = "rlml")
# RLMerr robust version
slmtest(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region +
                sub_nat + sub_cee,
        data = df2, listw = sp.matl, model = "within", effect = "twoways",
        test = "rlme")
## choose the most significative test

## models
# SAR (for network spatial matrix)
summary(
        spml(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region +
                     sub_nat + sub_cee,
             data = df2, index = c("dep", "period"), model = "within", effect = "twoways",
             lag = TRUE, listw = sp.matl, error = "none")
)

# SEM (for geo spatial matrix)
summary(
        spml(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region +
                     sub_nat + sub_cee,
             data = df2, index = c("dep", "period"), model = "within", effect = "twoways",
             lag = FALSE, listw = sp.matl, error = "b")
)



AICsplm = function(object, k=2, criterion=c("AIC", "BIC")){
        sp = summary(object)
        l = sp$logLik
        np = length(coef(sp))
        N = nrow(sp$model)
        if (sp$effects=="sptpfe") {
                n = length(sp$res.eff[[1]]$res.sfe)
                T = length(sp$res.eff[[1]]$res.tfe)
                np = np+n+T
        }
        if (sp$effects=="spfe") {
                n = length(sp$res.eff[[1]]$res.sfe)
                np = np+n+1
        }
        if (sp$effects=="tpfe") {
                T = length(sp$res.eff[[1]]$res.tfe)
                np = np+T+1
        }
        if(criterion=="AIC"){
                aic = -2*l+k*np
                names(aic) = "AIC"
                return(aic)
        }
        if(criterion=="BIC"){
                bic = -2*l+log(N)*np
                names(bic) = "BIC"
                return(bic)
        }
}



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

# (SDEM)
summary(
        model3.SDEM <- spml(net_density ~ treatment_int : group + gdp + dird +
                                    sub_region + sub_nat + sub_cee + treatment_int_SL +
                                    gdp_SL + dird_SL + sub_region_SL + sub_nat_SL + 
                                    sub_cee_SL,
                       data = df3, index = c("dep", "period"), model = "within",
                       effect = "twoways", lag = FALSE, listw = sp.matl, error = "b",
                       LeeYu = TRUE, Hess = FALSE)
        )
summary(model3.SDEM)$rsqr



# SDM
summary(
        model3.SDM <- spml(net_density ~ treatment_int : group + gdp + dird +
                                   sub_region + sub_nat + sub_cee + treatment_int_SL +
                                   gdp_SL + dird_SL + sub_region_SL + sub_nat_SL + 
                                   sub_cee_SL, data = df3, index = c("dep", "period"),
                           model = "within", effect = "twoways", lag = TRUE,
                           listw = sp.matl, error = "none", quiet = FALSE,
                           LeeYu = TRUE, Hess = FALSE)
        )
(4*94)-(16)-94-4+1
s2.model3.SDM <- sum(model3.SDM$resid^2)/263
s2.model3.SDM
s.model3.SDM = sqrt(s2.model3.SDM)
s.model3.SDM
R2.model3.SDM = 1 - var(model3.SDM$resid) / var(df3$net_density)
R2.model3.SDM
#or
summary(model3.SDM)$rsqr
summary(model3.SDM)$effects

AICsplm(model3.SDM, criterion = "BIC")





# SLX
summary(
        model3.SLX <- plm(net_density ~ treatment_int : group + gdp + dird + 
                                  sub_region + sub_nat + sub_cee + treatment_int_SL +
                                  gdp_SL + dird_SL + sub_region_SL + sub_nat_SL + sub_cee_SL,
                          data = df3, index = c("dep", "period"), model = "within",
                          effect = "twoways")
        )

# SAR
summary(
        model3.SAR <- spml(net_density ~ treatment_int : group + gdp + dird +
                                   sub_region + sub_nat + sub_cee,
                           data = df3, index = c("dep", "period"),
                           model = "within", effect = "twoways", lag = TRUE,
                           listw = sp.matl, error = "none", LeeYu = TRUE, Hess = FALSE)
)

coeftest(model3)

## calculate impact measures
time <- length(unique(df3$period))
set.seed(1234)
imps <- impacts(model3, listw = sp.matl, time = time)
summary(imps, zstats = TRUE, short = TRUE)

## fixed effects
effects(model3)

## moments

summary(
        model4 <- spgm(fragmentation_index ~ treatment_int : group + gdp + dird + sub_region +
                               sub_nat + sub_cee +
                               treatment_int_SL + gdp_SL + dird_SL + sub_region_SL +
                               sub_nat_SL + sub_cee_SL,
             data = df3, index = c("dep", "period"), model = "within",
             listw = sp.matl, moments = "fullweights", lag = FALSE,
             spatial.error = TRUE)
)

effects(model4)

# https://cran.r-project.org/web/packages/export/export.pdf











