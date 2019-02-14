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
library(plm)
library(AER)

## Fixed Effects Regression
# estimate the fixed effects regression with plm()
# all data set
model1 <- plm(share_national_nodes ~ treatment_int, data = df1, index = c("region", "period"),
              model = "within")
summary(model1)
# print summary using robust standard errors
coeftest(model1, vcov. = vcovHC, type = "HC1")
# data set for period = 3 or 4
model2 <- plm(share_national_nodes ~ treatment_int, data = subset(df1, period %in% c("3", "4")),
              index = c("region", "period"), model = "within")
summary(model2)
# print summary using robust standard errors
coeftest(model2, vcov. = vcovHC, type = "HC1")



## DiD estimation
# Create a dummy variable to indicate the time when the treatment started
df1$time <- ifelse(as.character(df1$period) %in% c("1", "2"), 0, 1)
# Create a dummy variable to identify the group exposed to the treatment. 
df1$treatment_int
# Create an interaction between time and treated. We will call this interaction ‘did’.
df1$did <- df1$time * df1$treatment_int

# Estimating the DID estimator (with lm)
model1 <- lm(share_national_nodes ~ did + time, data = df1)
summary(model1)


# Estimating the DID estimator (with plm)
model1 <- plm(share_national_nodes ~ did + time, data = df1,
              index = c("region", "period"), model = "within")
summary(model1)
coeftest(model1, vcov. = vcovHC, type = "HC1")

























