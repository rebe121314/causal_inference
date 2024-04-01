
#install.packages('sensemakr')
# loads package
library(sensemakr)

# loads data
data("darfur")
?darfur

#install.packages("Matching")
library(Matching)

##Propensity 

# Assuming 'darfur' data is loaded and the 'Matching' library is loaded
glm1 <- glm(directlyharmed ~ age + female + hhsize_darfur + herder_dar+ farmer_dar + pastvoted + female*hhsize_darfur + herder_dar*pastvoted + farmer_dar*pastvoted, 
                family = binomial, data = darfur)



# Predict propensity scores using the fitted model
Xp  <- glm1$fitted
Y  <- darfur$peacefactor
Tr  <- darfur$directlyharmed

rr <- Match(Y = NULL, Tr = Tr, X = Xp, estimand = "ATT", M = 1)

mb <- MatchBalance(directlyharmed ~  age + female + hhsize_darfur + herder_dar+ farmer_dar + pastvoted + female*hhsize_darfur + herder_dar*pastvoted + farmer_dar*pastvoted,
                    data=darfur, match.out=rr, nboots=500)

rr_t <- Match(Y = Y, Tr = Tr, X = Xp, estimand = "ATT", M = 1)
summary(rr_t)

#Sensitivity

mi <- rr_t$index.treated
mci<- rr_t$index.control

md <- darfur[c(mi, mci),]
md


#linear model
lm1 <- lm(peacefactor ~ directlyharmed + age + female + hhsize_darfur + herder_dar+ 
            farmer_dar + pastvoted + female*hhsize_darfur + herder_dar*pastvoted 
          + farmer_dar*pastvoted, data = md)
summary(lm1)

lm1.sensitivity <- sensemakr(model = lm1, 
                              treatment = "directlyharmed",
                              benchmark_covariates = "female",
                              kd = 1:5,
                              alpha = 0.05, 
                              reduce = TRUE)
lm1.sensitivity
plot(lm1.sensitivity)
plot(lm1.sensitivity, type = "extreme")
plot(lm1.sensitivity, sensitivity.of = "t-value")
#ovb_minimal_reporting(lm1.sensitivity, format = "latex")



### GENETIC MATCHING - Original
data("darfur")
set.seed(333)


X <- cbind(darfur$age, darfur$female, darfur$hhsize_darfur, darfur$herder_dar,  darfur$farmer_dar, darfur$pastvoted)
BalanceMat <- cbind(darfur$age ,darfur$female, darfur$hhsize_darfur,  darfur$herder_dar,  darfur$farmer_dar, darfur$pastvoted, darfur$female*darfur$hhsize_darfur, darfur$herder_dar*darfur$pastvoted,darfur$farmer_dar*darfur$pastvoted)
#install.packages('rgenoud')

genout <- GenMatch(Tr=Tr, X=X, BalanceMatrix=BalanceMat, estimand="ATT", M=1,  pop.size = 20, nboots = 100)
genout

mout <- Match(Y=NULL, Tr=Tr, X=X, estimand="ATT", Weight.matrix=genout)
summary(mout)
MatchBalance(directlyharmed ~  age + female + hhsize_darfur + herder_dar+ farmer_dar + pastvoted + female*hhsize_darfur + herder_dar*pastvoted + farmer_dar*pastvoted,
             data=darfur, match.out=mout, nboots=500)

mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=genout)
summary(mout)

mi2 <- mout$index.treated
mci2<- mout$index.control

md2 <- darfur[c(mi2, mci2),]



lm2 <- lm(peacefactor ~ directlyharmed + age + female + hhsize_darfur + herder_dar+ farmer_dar + pastvoted + female*hhsize_darfur + herder_dar*pastvoted + farmer_dar*pastvoted, data = md2)

lm2.sensitivity <- sensemakr(model = lm2, 
                             treatment = "directlyharmed",
                             benchmark_covariates = "female",
                             kd = 1:5,
                             reduce = TRUE)
summary(lm2)
lm2.sensitivity
plot(lm2.sensitivity)
ovb_minimal_reporting(lm2.sensitivity, format = "latex")
plot(lm2.sensitivity, type = "extreme")
plot(lm2.sensitivity, sensitivity.of = "t-value")


##GENMATCHING 2 - Extended 
data("darfur")
set.seed(333)
X <- cbind(darfur$age, darfur$female, darfur$hhsize_darfur, darfur$herder_dar,  darfur$farmer_dar, darfur$pastvoted)

BalanceMat2 <- cbind(darfur$age, darfur$age ,darfur$female, darfur$hhsize_darfur,  darfur$herder_dar,  darfur$farmer_dar, darfur$pastvoted, darfur$female*darfur$hhsize_darfur, darfur$herder_dar*darfur$pastvoted, darfur$farmer_dar*darfur$pastvoted, darfur$female*darfur$age ,I(darfur$age^2),darfur$hhsize_darfur*darfur$hhsize_darfur)

genout2 <- GenMatch(Tr=Tr, X=X, BalanceMatrix=BalanceMat2, estimand="ATT", M=1, pop.size = 20, nboots = 100)
mout2 <- Match(Y=NULL, Tr=Tr, X=X, estimand="ATT", Weight.matrix=genout2)
summary(mout2)

MatchBalance(directlyharmed ~  age + female + hhsize_darfur + herder_dar+ farmer_dar + pastvoted + female*hhsize_darfur + herder_dar*pastvoted + farmer_dar*pastvoted,
             data=darfur, match.out=mout2, nboots=500)

mout2 <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=genout2)
summary(mout2)

mi3 <- mout2$index.treated
mci3<- mout2$index.control

md3 <- darfur[c(mi3, mci3),]


lm3 <- lm(peacefactor ~ directlyharmed + age +female*hhsize_darfur + herder_dar*pastvoted + farmer_dar*pastvoted + I(hhsize_darfur^2) + female*age + I(age^2),data = md3)

lm3.sensitivity <- sensemakr(model = lm3, 
                             treatment = "directlyharmed",
                             benchmark_covariates = "female",
                             kd = 1:5,
                             reduce = TRUE)
summary(lm3)
lm3.sensitivity
plot(lm3.sensitivity)
plot(lm3.sensitivity, type = "extreme")
plot(lm3.sensitivity, sensitivity.of = "t-value")
ovb_minimal_reporting(lm3.sensitivity, format = "latex")



###PART II
## Syntehetic cosntrol
#install.packages('Synth')
library(Synth)
data("basque")
?basque
library(tidyverse)

data(basque)

#where's cantabria?
unique(basque$regionname)
cantabria_id <- with(basque, basque$regionno[regionname == "Cantabria"])
unique(cantabria_id)




# dataprep: prepare data for synth
dataprep.out <-
  dataprep(
    foo = basque
    ,predictors= c("school.illit",
                   "school.prim",
                   "school.med",
                   "school.high",
                   "school.post.high"
                   ,"invest"
    )
    ,predictors.op = c("mean")
    ,dependent = c("gdpcap")
    ,unit.variable = c("regionno")
    ,time.variable = c("year")
    ,special.predictors = list(
      list("gdpcap",1960:1969,c("mean")),
      list("sec.agriculture",seq(1961,1969,2),c("mean")),
      list("sec.energy",seq(1961,1969,2),c("mean")),
      list("sec.industry",seq(1961,1969,2),c("mean")),
      list("sec.construction",seq(1961,1969,2),c("mean")),
      list("sec.services.venta",seq(1961,1969,2),c("mean")),
      list("sec.services.nonventa",seq(1961,1969,2),c("mean")),
      list("popdens",1969,c("mean")))
    ,treatment.identifier = 7
    ,controls.identifier = c(2:6, 8:16, 18)
    ,time.predictors.prior = c(1964:1969)
    ,time.optimize.ssr = c(1960:1969)
    ,unit.names.variable = c("regionname")
    ,time.plot = c(1955:1997)
  )
# 1. combine highest and second highest
# schooling category and eliminate highest category
dataprep.out$X1["school.high",] <-
  dataprep.out$X1["school.high",] +
  dataprep.out$X1["school.post.high",]
dataprep.out$X1 <-
  as.matrix(dataprep.out$X1[
    -which(rownames(dataprep.out$X1)=="school.post.high"),])
dataprep.out$X0["school.high",] <-
  dataprep.out$X0["school.high",] +
  dataprep.out$X0["school.post.high",]
dataprep.out$X0 <-
  dataprep.out$X0[
    -which(rownames(dataprep.out$X0)=="school.post.high"),]
# 2. make total and compute shares for the schooling catgeories
lowest <- which(rownames(dataprep.out$X0)=="school.illit")
highest <- which(rownames(dataprep.out$X0)=="school.high")
dataprep.out$X1[lowest:highest,] <-
  (100 * dataprep.out$X1[lowest:highest,]) /
  sum(dataprep.out$X1[lowest:highest,])
dataprep.out$X0[lowest:highest,] <-
  100 * scale(dataprep.out$X0[lowest:highest,],
              center=FALSE,
              scale=colSums(dataprep.out$X0[lowest:highest,])
  )
# run synth
synth.out <- synth(data.prep.obj = dataprep.out)
# Get result tables
synth.tables <- synth.tab(
  dataprep.res = dataprep.out,
  synth.res = synth.out
  )
# results tables:
print(synth.tables)
# plot results:
# path
path.plot(
  synth.res = synth.out,
  dataprep.res = dataprep.out,
  Ylab = "real per-capita GDP (1986 USD, thousand)",
  Xlab = "year",
  Ylim = c(0, 13),
  Legend = c("Cantabria", "synthetic Cantabria")
)
abline(v=1970, col="red", lwd=2) 

# Plot results: Gaps plot
gaps.plot(
  synth.res = synth.out,
  dataprep.res = dataprep.out,
  Ylab = "gap in real per-capita GDP (1986 USD, thousand)",
  Xlab = "year",
  Ylim = c(-1.5, 1.5)
)
abline(v=1970, col="red", lwd=2) 


          