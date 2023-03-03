# https://smolski.github.io/livroavancado/index.html
require(ggplot2)
require(sandwich)
require(msm)




glm.RP <- function(GLM.RESULT, digits = 3) {
  
  if (GLM.RESULT$family$family == "binomial") {
    LABEL <- "OR"
  } else if (GLM.RESULT$family$family == "poisson") {
    LABEL <- "RP"
  } else {
    stop("Not logistic or Poisson model")
  }
  
  COEF      <- stats::coef(GLM.RESULT)
  CONFINT   <- stats::confint(GLM.RESULT)
  TABLE     <- cbind(coef=COEF, CONFINT)
  TABLE.EXP <- round(exp(TABLE), digits)
  
  colnames(TABLE.EXP)[1] <- LABEL
  
  TABLE.EXP


### Cobertura mediana


m_cob <- glm(casos ~ offset(log(pop2021)) + cob_mediana, family="poisson", micro.sf)

cov.m1 <- vcovHC(m_cob, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(glm.RP(m_cob),
               LL = round(exp(coef(m_cob) - 1.96 * std.err), 3),
               UL = round(exp(coef(m_cob) + 1.96 * std.err), 3),
               "p-value" = round(2 * pnorm(abs(coef(m_cob)/std.err), lower.tail=FALSE), 2))

r.est

dispersiontest(m_cob)



### IVS


m_ivs <- glm(casos ~ offset(log(pop2021)) + ivs, family="poisson", micro.sf)

cov.m1 <- vcovHC(m_ivs, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(glm.RP(m_ivs),
               LL = round(exp(coef(m_ivs) - 1.96 * std.err), 3),
               UL = round(exp(coef(m_ivs) + 1.96 * std.err), 3),
               "p-value" = round(2 * pnorm(abs(coef(m_ivs)/std.err), lower.tail=FALSE), 2))

r.est


### IDHM


m_idhm <- glm(casos ~ offset(log(pop2021)) + idhm, family="poisson", micro.sf)

cov.m1 <- vcovHC(m_idhm, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(glm.RP(m_idhm),
               LL = round(exp(coef(m_idhm) - 1.96 * std.err), 3),
               UL = round(exp(coef(m_idhm) + 1.96 * std.err), 3),
               "p-value" = round(2 * pnorm(abs(coef(m_idhm)/std.err), lower.tail=FALSE), 2))

r.est


### Mortalidade < 5 anos

m_mort5 <- glm(casos ~ offset(log(pop2021)) + mort5, family="poisson", micro.sf)

cov.m1 <- vcovHC(m_mort5, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(glm.RP(m_mort5),
               LL = round(exp(coef(m_mort5) - 1.96 * std.err), 3),
               UL = round(exp(coef(m_mort5) + 1.96 * std.err), 3),
               "p-value" = round(2 * pnorm(abs(coef(m_mort5)/std.err), lower.tail=FALSE), 2))

r.est

### Proporção de analfabetos > 18 anos

m_analf18 <- glm(casos ~ offset(log(pop2021)) + analf18, family="poisson", micro.sf)

cov.m1 <- vcovHC(m_analf18, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(glm.RP(m_analf18),
               LL = round(exp(coef(m_analf18) - 1.96 * std.err), 3),
               UL = round(exp(coef(m_analf18) + 1.96 * std.err), 3),
               "p-value" = round(2 * pnorm(abs(coef(m_analf18)/std.err), lower.tail=FALSE), 2))

r.est

### Gini

m_gini <- glm(casos ~ offset(log(pop2021)) + gini, family="poisson", micro.sf)

cov.m1 <- vcovHC(m_gini, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(glm.RP(m_gini),
               LL = round(exp(coef(m_gini) - 1.96 * std.err), 3),
               UL = round(exp(coef(m_gini) + 1.96 * std.err), 3),
               "p-value" = round(2 * pnorm(abs(coef(m_gini)/std.err), lower.tail=FALSE), 2))

r.est


### Modelo cheio 
  
m_cheio <- glm(casos ~ offset(log(pop2021)) + cob_mediana + ivs + idhm + mort5 + analf18 + gini, family="poisson", micro.sf) 
  
cov.m1 <- vcovHC(m_cheio, type="HC0") 
std.err <- sqrt(diag(cov.m1)) 
r.est <- cbind(glm.RP(m_cheio), 
                        LL = round(exp(coef(m_cheio) - 1.96 * std.err), 3), 
                        UL = round(exp(coef(m_cheio) + 1.96 * std.err), 3), 
                        "p-value" = round(2 * pnorm(abs(coef(m_cheio)/std.err), lower.tail=FALSE), 2)) 
  
r.est

# https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
# Install the package jtools if not already installed
library("jtools")
# you may be asked to install 'broom' and 'ggstance' packages as well
library("broom")
library("ggstance")

plot_coefs(m_cob, m_ivs, m_idhm, m_mort5, exp = TRUE, escale=T,
           coefs = c("Cob Vacinal %" = "cob_mediana", "IVS" = "ivs", "IDHM" = "idhm", "Mortalidade 5 anos" = "mort5"),
           scale = TRUE, robust = TRUE)

# https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
# https://andrewproctor.github.io/rcourse/module5.html


# Modelagem

|>data.matrix()


DM <-gw.dist(dp.locat=coordinates(micro.sp))

bw.gwr <- bw.ggwr(casos ~ 1,
                 data = micro.sp,
                 family = "poisson",
                 approach = "AICc",
                 kernel = "bisquare",
                 longlat=F,
                 adaptive = T,
                 dMat = DM)
bw.gwr

micro.coords <- st_coordinates(micro.utm) |>
  as.data.frame() |>
  select(X, Y) 


modelo <- gw.zi(casos ~ offset(log(pop2021)) + cob_mediana + ivs + idhm + mort5 + analf18 + gini,"poisson", micro.sp,400, kernel = 'adaptive', coordinates(micro.sp))

install.packages("geobr")
library(geobr)

regiao <- read_region(year = 2018)


