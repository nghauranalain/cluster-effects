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


# function from Spatial_Econometrics_with_R (Prof. Dr. Reinhold Kosfeld)
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


# function from https://stackoverflow.com/questions/46186527/how-to-calculate-bic-and-aic-for-a-gmm-model-in-r-using-plm
AICplm <- function(object, criterion) {
        # object is "plm", "panelmodel" 
        # Lets panel data has index :index = c("dep", "period")
        
        sp = summary(object)
        
        if(class(object)[1]=="plm"){
                u.hat <- residuals(sp) # extract residuals
                df <- cbind(as.vector(u.hat), attr(u.hat, "index"))
                names(df) <- c("resid", "dep", "period")
                c = length(levels(df$dep)) # extract country dimension 
                t = length(levels(df$period)) # extract time dimension 
                np = length(sp$coefficients[,1]) # number of parameters
                n.N = nrow(sp$model) # number of data
                s.sq  <- log( (sum(u.hat^2)/(n.N))) # log sum of squares
                
                # effect = c("individual", "time", "twoways", "nested"),
                # model = c("within", "random", "ht", "between", "pooling", "fd")
                
                # I am making example only with some of the versions:
                
                if (sp$args$model == "within" & sp$args$effect == "individual"){
                        n = c
                        np = np+n+1 # update number of parameters
                }
                
                if (sp$args$model == "within" & sp$args$effect == "time"){
                        T = t
                        np = np+T+1 # update number of parameters
                }
                
                if (sp$args$model == "within" & sp$args$effect == "twoways"){
                        n = c
                        T = t
                        np = np+n+T # update number of parameters
                }
                aic <- round(2*np  +  n.N * (log(2*pi) + s.sq  + 1 ),1)
                bic <- round(log(n.N)*np  +  n.N * (log(2*pi) + s.sq  + 1 ), 1)
                
                if(criterion=="AIC"){
                        names(aic) = "AIC"
                        return(aic)
                }
                if(criterion=="BIC"){
                        names(bic) = "BIC"
                        return(bic)
                }
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

# SDEM
summary(
        model3.SDEM <- spml(net_hierarchy ~ treatment_int : group + gdp + dird +
                                    sub_region + sub_nat + sub_cee + treatment_int_SL : group +
                                    gdp_SL + dird_SL + sub_region_SL + sub_nat_SL + 
                                    sub_cee_SL,
                       data = df3, index = c("dep", "period"), model = "within",
                       effect = "twoways", lag = FALSE, listw = sp.matl, spatial.error = "b",
                       LeeYu = TRUE, Hess = FALSE)
        )
summary(model3.SDEM)$rsqr


# SDM
summary(
        model3.SDM <- spml(net_density ~ treatment_int : group + gdp + dird +
                                   sub_region + sub_nat + sub_cee + treatment_int_SL : group +
                                   gdp_SL + dird_SL + sub_region_SL + sub_nat_SL + 
                                   sub_cee_SL, data = df3, index = c("dep", "period"),
                           model = "within", effect = "twoways", lag = TRUE,
                           listw = sp.matl, spatial.error = "none", LeeYu = TRUE, Hess = FALSE)
        )
summary(model3.SDM)$rsqr



## calculate impact measures
time <- length(unique(df3$period))
set.seed(1234)
imps <- impacts(model3.SDM, listw = sp.matl, time = time)
summary(imps, zstats = TRUE, short = TRUE)

# test
model3.SDM$spat.coef
iIrW  <- invIrW(sp.matl, model3.SDM$spat.coef)
model3.SDM$coefficients
model3.SDM$coefficients[2]
model3.SDM$coefficients[7]

S_gdp <- iIrW %*% ((model3.SDM$coefficients[2] * diag(94)) - # diag(nrow(df3))
                           (model3.SDM$coefficients[7] * listw2mat(sp.matl)))

S_gdp <- iIrW %*% (model3.SDM$coefficients[2] * diag(94)) # gives wrong impact got from impacts()

# Then you can get the direct and total impacts in the usual way: 
dir_gdp <- sum(diag(S_gdp))/94 # mean
tot_gdp <- sum(c(S_gdp))/94
indir_gdp <- tot_gdp - dir_gdp

mean(diag(S_gdp))
summary(diag(S_gdp))
sd(diag(S_gdp))
LaplacesDemon::p.interval(diag(S_gdp), prob = 0.95)
LaplacesDemon::p.interval(c(S_gdp), prob = 0.95)

impacts
methods(impacts)
getS3method("impacts", "sarlm") 


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
AICsplm(model3.SDM, criterion = "AIC")
AICsplm(model3.SAR, criterion = "AIC")
AICsplm(model3.SDM, criterion = "BIC")
AICsplm(model3.SAR, criterion = "BIC")
# SDM wins SAR

## SDM vs SLX 
AICsplm(model3.SDM, criterion = "AIC")
AICplm(model3.SLX, criterion = "AIC")
AICsplm(model3.SDM, criterion = "BIC")
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


