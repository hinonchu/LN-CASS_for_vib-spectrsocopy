library('lhs')
# simulation data ----
set.seed(55)
library(tidyverse)

int.10 <- 0.11
int.50 <- -0.6
int.100 <- -0.3
int.150 <- 0.4
int.1000 <- -0.2

beta.10 <- c(0,0,-2,-1.6,-2,0,0,1,1.2,1)

beta.50 <- c(rep(0,5),rep(2,5),rep(0,5),(rep(-1,5) + rnorm(5,mean = 0,sd = 0.1)),(rep(1,5) + rnorm(5,mean = 0,sd = 0.1)),(rep(-1,5) + rnorm(5,mean = 0,sd = 0.1)),rep(randomLHS(1,1),5),rep(randomLHS(1,1),5),rep(randomLHS(1,1),5),rep(randomLHS(1,1),5))

beta.100 <- c(rep(0,30),rep(2,10),rep(0,20),rep(2,10),(rep(-1,10) + rnorm(10,mean = 0,sd = 0.1)),rep(0,15),-0.5,-0.5,0,-0.5,-0.5)

beta.150 <- c(rep(0,25),rep(1,5),rep(0,25),(rep(2,5) + rnorm(5,mean = 0,sd = 0.1)), rep(0,25),-1,-1,-1,-2,-1,0.5,.5,.5,2,.5,rep(0,25),rep(3,15),rep(0,15))

beta.1000 <- c(rep(beta.100,2), rep(beta.150,2), rep(beta.50,10))

cov.10 <- randomLHS(100,10)

cov.50 <- randomLHS(100,50)

cov.100 <- randomLHS(100,100)

cov.150<- randomLHS(100,150)

cov.1000 <- randomLHS(100,1000)

y.10 <- int.10 + cov.10 %*% beta.10 + rnorm(100,mean = 0,sd = 0.5)

y.50 <- int.50 + cov.50 %*% beta.50 + rnorm(100,mean = 0,sd = 0.5)

y.100 <- int.100 + cov.100 %*% beta.100 +rnorm(100,mean = 0, sd = 0.5)

y.150 <- int.150 + cov.150 %*% beta.150 + rnorm(100,mean = 0, sd = 0.5)

y.1000 <- int.1000 + cov.1000 %*% beta.1000 +rnorm(100,mean = 0, sd = 0.5)


# p = 10, n = 100
library('SGL')
grplassofit <- cvSGL(data = list(x = cov.10,y = y.10),
                     index = c(rep(1,2),rep(2,2),rep(3,2),rep(4,2),rep(5,2)),
                     maxit = 10000)

library(rstanarm)
hsfit <- stan_glm(y.10 ~ . ,
                  data = data.frame(cov.10,y.10),
                  prior = hs(df = 2,global_scale = 0.2),
                  family = gaussian(),
                  prior_intercept = normal(0,10),
                  adapt_delta = 0.9999)

median.pars.hs <- apply(as.matrix(hsfit),2,median)[2:11]

lmfit <- lm(y.10 ~ . , data = data.frame(cov.10,y.10))

library(glmnet)
glmfit <- cv.glmnet(x = cov.10,y = y.10)
library(rstan)

N_obs <- 100
ncov <- length(cov.10[1,])
ngroup_n <- ncov/5

cassfit <- stan("simcasshered.stan",
                data = list(N = 100,
                            ncov = 10,
                            y = y.10[,1],
                            X = cov.10,
                            sigma_indic = 10,
                            mu_indic = 0,
                            ngroup = ngroup_n),
                control = list(adapt_delta = 0.999, max_treedepth = 15),
                chains = 4, cores = 4, iter = 2000)

median.pars.cass <- summary(cassfit,pars = "beta")$summary[,6]


newdf <- data.frame(actual = rep(beta.10,5),
                    fitted = c(grplassofit$fit$beta[,8],
                               median.pars.hs,
                               coef(lmfit)[-1],
                               coef(glmfit)[-1],
                               median.pars.cass))

library('reshape2')
newdf <- data.frame(melt(newdf),class = rep(rep(c("SGL","Horseshoe","OLS","LASSO","LN-CASS"),each = 10),2))

ggplot(data = newdf, 
       aes(x = rep(1:10,10),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "")


library(pROC)

AUC_CASS <- multiclass.roc(beta.10, median.pars.cass)$auc

MAE_CASS <- mean(((beta.10-median.pars.cass)^2)^0.5)

AUC_HS <- multiclass.roc(beta.10,median.pars.hs)$auc

MAE_HS <- mean(((beta.10-median.pars.hs)^2)^0.5)

AUC_SGL <- multiclass.roc(beta.10, grplassofit$fit$beta[,8])$auc

MAE_SGL <- mean(((beta.10-grplassofit$fit$beta[,8])^2)^0.5)

AUC_OLS <- multiclass.roc(beta.10,coef(lmfit)[-1])$auc

MAE_OLS <- mean(((beta.10-coef(lmfit)[-1])^2)^0.5)

AUC_LASSO <- multiclass.roc(beta.10,coef(glmfit)[-1])$auc

MAE_LASSO <- mean(((beta.10-coef(glmfit)[-1])^2)^0.5)

# run some fitting algorithms

#### p = 50, n = 100 ----
#uses SGL library
# p = 50, n = 500
library('SGL')
grplassofit_p50 <- cvSGL(data = list(x = cov.50,y = y.50),
                     index = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10)),
                     maxit = 50000)

library(rstanarm)
hsfit_p50 <- stan_glm(y.50 ~ . ,
                  data = data.frame(cov.50,y.50),
                  prior = hs(df = 2,global_scale = 0.2),
                  family = gaussian(),
                  prior_intercept = normal(0,50),
                  adapt_delta = 0.9999)

median.pars.hs_p50 <- apply(as.matrix(hsfit_p50),2,median)[2:(length(beta.50)+1)]

lmfit_p50 <- lm(y.50 ~ . , data = data.frame(cov.50,y.50))

library(glmnet)
glmfit_p50 <- cv.glmnet(x = cov.50,y = y.50)
library(rstan)

N_obs <- 100
ncov <- length(cov.50[1,])
ngroup_n <- ncov/5

cassfit_p50 <- stan("simcasshered.stan",
                data = list(N = 100,
                            ncov = 50,
                            y = y.50[,1],
                            X = cov.50,
                            sigma_indic = 50,
                            mu_indic = 0,
                            ngroup = ngroup_n),
                control = list(adapt_delta = 0.999, max_treedepth = 15),
                chains = 4, cores = 4, iter = 2000)

median.pars.cass_p50 <- summary(cassfit_p50,pars = "beta")$summary[,6]


newdf_p50 <- data.frame(actual = rep(beta.50,5),
                    fitted = c(grplassofit_p50$fit$beta[,8],
                               median.pars.hs_p50,
                               coef(lmfit_p50)[-1],
                               coef(glmfit_p50)[-1],
                               median.pars.cass_p50))

library('reshape2')
newdf_p50 <- data.frame(melt(newdf_p50),class = rep(rep(c("SGL","Horseshoe","OLS","LASSO","LN-CASS"),each = 50),2))

ggplot(data = newdf_p50, 
       aes(x = rep(1:50,10),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "")


library(pROC)

AUC_CASS_p50<- multiclass.roc(beta.50, median.pars.cass_p50)$auc

MAE_CASS_p50 <- mean(((beta.50-median.pars.cass_p50)^2)^0.5)

AUC_HS_p50 <- multiclass.roc(beta.50,median.pars.hs_p50)$auc

MAE_HS_p50 <- mean(((beta.50-median.pars.hs_p50)^2)^0.5)

AUC_SGL_p50 <- multiclass.roc(beta.50, grplassofit_p50$fit$beta[,8])$auc

MAE_SGL_p50 <- mean(((beta.50-grplassofit_p50$fit$beta[,8])^2)^0.5)

AUC_OLS_p50 <- multiclass.roc(beta.50,coef(lmfit_p50)[-1])$auc

MAE_OLS_p50 <- mean(((beta.50-coef(lmfit_p50)[-1])^2)^0.5)

AUC_LASSO_p50 <- multiclass.roc(beta.50,coef(glmfit_p50)[-1])$auc

MAE_LASSO_p50 <- mean(((beta.50-coef(glmfit_p50)[-1])^2)^0.5)


#### p = 100, n = 100 ----

# grplassofit1 <- cvSGL(data = list(x = cov.100,y = y.100),
#                       index = rep(1:14,each = 5),
#                       maxit = 10000)
# 
# hsfit1 <- stan_glm(y.100 ~ . ,
#                    data = data.frame(cov.100,y.100),
#                    prior = hs(df = 2),
#                    family = gaussian(),
#                    prior_intercept = normal(0,10),
#                    adapt_delta = 0.9999)
# 
# median.pars.hs1 <- apply(as.matrix(hsfit1),2,median)[2:71]
# 
# 
# lmfit1 <- lm(y.100 ~ . , data = data.frame(cov.100,y.100))
# 
# glmfit1 <- cv.glmnet(x = cov.100,y = y.100)
# 
# cassfit1 <- stan("simcasshered.stan",
#                  data = list(N = 100,
#                              ncov = 100,
#                              y = y.100,
#                              X = cov.100,
#                              sigma_indic = 10,
#                              mu_indic = 0,
#                              ngroup = 14),
#                  control = list(adapt_delta = 0.9999),
#                  chains = 4, cores = 4, iter = 5000)
# 
# median.pars.cass1 <- summary(cassfit1,pars = "beta")$summary[,6]
# 
# 
# newdf <- data.frame(actual = rep(beta.100,5),
#                     fitted = c(grplassofit1$fit$beta[,8],
#                                median.pars.hs1,
#                                coef(lmfit1)[-1],
#                                coef(glmfit1)[-1],
#                                median.pars.cass1))
# newdf <- data.frame(melt(newdf),class = rep(rep(c("SGL","Horseshoe","OLS","LASSO","LN-CASS"),each = 100),2))
# 
# ggplot(data = newdf, 
#        aes(x = rep(1:100,10),y = value,col = variable)) + 
#   geom_point(size = 1) + 
#   geom_line() + 
#   facet_wrap(~ class,nrow = 1,scales = "free_y") + 
#   scale_color_manual(labels=c("Actual", "Fitted"),
#                      values = c("blue","red")) +
#   labs(x = "Coefficient index", 
#        y = "Value of coefficient",
#        color = "")

library('SGL')
grplassofit_p100 <- cvSGL(data = list(x = cov.100,y = y.100),
                          index = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10)),
                          maxit = 100000)

library(rstanarm)
hsfit_p100 <- stan_glm(y.100 ~ . ,
                       data = data.frame(cov.100,y.100),
                       prior = hs(df = 2,global_scale = 0.2),
                       family = gaussian(),
                       prior_intercept = normal(0,100),
                       adapt_delta = 0.9999)

median.pars.hs_p100 <- apply(as.matrix(hsfit_p100),2,median)[2:(length(beta.100)+1)]

lmfit_p100 <- lm(y.100 ~ . , data = data.frame(cov.100,y.100))

library(glmnet)
glmfit_p100 <- cv.glmnet(x = cov.100,y = y.100)
library(rstan)

N_obs <- 100
ncov <- length(cov.100[1,])
ngroup_n <- ncov/5

cassfit_p100 <- stan("simcasshered.stan",
                     data = list(N = 100,
                                 ncov = 100,
                                 y = y.100[,1],
                                 X = cov.100,
                                 sigma_indic = 10,
                                 mu_indic = 0,
                                 ngroup = ngroup_n),
                     control = list(adapt_delta = 0.999, max_treedepth = 15),
                     chains = 4, cores = 4, iter = 2000)

median.pars.cass_p100 <- summary(cassfit_p100,pars = "beta")$summary[,6]


newdf_p100 <- data.frame(actual = rep(beta.100,5),
                         fitted = c(grplassofit_p100$fit$beta[,8],
                                    median.pars.hs_p100,
                                    coef(lmfit_p100)[-1],
                                    coef(glmfit_p100)[-1],
                                    median.pars.cass_p100))

library('reshape2')
newdf_p100 <- data.frame(melt(newdf_p100),class = rep(rep(c("SGL","Horseshoe","OLS","LASSO","LN-CASS"),each = 100),2))

ggplot(data = newdf_p100, 
       aes(x = rep(1:100,10),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "") +
  theme(text = element_text(size = 20))


library(pROC)

AUC_CASS_p100<- multiclass.roc(beta.100, median.pars.cass_p100)$auc

MAE_CASS_p100 <- mean(((beta.100-median.pars.cass_p100)^2)^0.5)

AUC_HS_p100 <- multiclass.roc(beta.100,median.pars.hs_p100)$auc

MAE_HS_p100 <- mean(((beta.100-median.pars.hs_p100)^2)^0.5)

AUC_OLS_p100 <- multiclass.roc(beta.100,coef(lmfit_p100)[-1])$auc

MAE_OLS_p100 <- mean(((na.exclude(beta.100-coef(lmfit_p100)[-1]))^2)^0.5)

AUC_SGL_p100 <- multiclass.roc(beta.100, grplassofit_p100$fit$beta[,8])$auc

MAE_SGL_p100 <- mean(((beta.100-grplassofit_p100$fit$beta[,8])^2)^0.5)

AUC_LASSO_p100 <- multiclass.roc(beta.100,coef(glmfit_p100)[-1])$auc

MAE_LASSO_p100 <- mean(((beta.100-coef(glmfit_p100)[-1])^2)^0.5)

### n = 100 p = 150
library('SGL')
grplassofit_p150 <- cvSGL(data = list(x = cov.150,y = y.150),
                          index = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10),rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,10)),
                          maxit = 100000)


library(rstanarm)
hsfit_p150 <- stan_glm(y.150 ~ . ,
                       data = data.frame(cov.150,y.150),
                       prior = hs(df = 2,global_scale = 0.2),
                       family = gaussian(),
                       prior_intercept = normal(0,150),
                       adapt_delta = 0.9999)

median.pars.hs_p150 <- apply(as.matrix(hsfit_p150),2,median)[2:(length(beta.150)+1)]

lmfit_p150 <- lm(y.150 ~ . , data = data.frame(cov.150,y.150))

library(glmnet)
glmfit_p150 <- cv.glmnet(x = cov.150,y = y.150)
library(rstan)

N_obs <- 150
ncov <- length(cov.150[1,])
ngroup_n <- ncov/5

cassfit_p150 <- stan("simcasshered.stan",
                     data = list(N = 100,
                                 ncov = 150,
                                 y = y.150[,1],
                                 X = cov.150,
                                 sigma_indic = 10,
                                 mu_indic = 0,
                                 ngroup = ngroup_n),
                     control = list(adapt_delta = 0.999, max_treedepth = 15),
                     chains = 4, cores = 4, iter = 2000)

median.pars.cass_p150 <- summary(cassfit_p150,pars = "beta")$summary[,6]


newdf_p150 <- data.frame(actual = rep(beta.150,4),
                         fitted = c(grplassofit_p150$fit$beta[,8],
                                    median.pars.hs_p150,
                                    coef(glmfit_p150)[-1],
                                    median.pars.cass_p150))

library('reshape2')
newdf_p150 <- data.frame(melt(newdf_p150),class = rep(rep(c("SGL","Horseshoe","LASSO","LN-CASS"),each = 150),2))

ggplot(data = newdf_p150, 
       aes(x = rep(1:150,8),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "") +
  theme(text = element_text(size = 20))


library(pROC)

AUC_CASS_p150<- multiclass.roc(beta.150, median.pars.cass_p150)$auc

MAE_CASS_p150 <- mean(((beta.150-median.pars.cass_p150)^2)^0.5)

AUC_HS_p150 <- multiclass.roc(beta.150,median.pars.hs_p150)$auc

MAE_HS_p150 <- mean(((beta.150-median.pars.hs_p150)^2)^0.5)

AUC_SGL_p150 <- multiclass.roc(beta.150, grplassofit_p150$fit$beta[,8])$auc

MAE_SGL_p150 <- mean(((beta.150-grplassofit_p50$fit$beta[,8])^2)^0.5)

AUC_LASSO_p150 <- multiclass.roc(beta.150,coef(glmfit_p150)[-1])$auc

MAE_LASSO_p150 <- mean(((beta.150-coef(glmfit_p150)[-1])^2)^0.5)

# # n = 100 p = 1000
# 
# library('SGL')
# grplassofit_p1000 <- cvSGL(data = list(x = cov.100,y = y.100),
#                            index = c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100)),
#                            maxit = 100000)
# 
# library(rstanarm)
# hsfit_p1000 <- stan_glm(y.1000 ~ . ,
#                        data = data.frame(cov.1000,y.1000),
#                        prior = hs(df = 2,global_scale = 0.2),
#                        family = gaussian(),
#                        prior_intercept = normal(0,1000),
#                        adapt_delta = 0.9999)
# 
# median.pars.hs_p1000 <- apply(as.matrix(hsfit_p1000),2,median)[2:(length(beta.1000)+1)]
# 
# lmfit_p1000 <- lm(y.1000 ~ . , data = data.frame(cov.1000,y.1000))
# 
# library(glmnet)
# glmfit_p1000 <- cv.glmnet(x = cov.1000,y = y.1000)
# library(rstan)
# 
# N_obs <- 100
# ncov <- length(cov.1000[1,])
# ngroup_n <- ncov/5
# 
# cassfit_p1000 <- stan("simcasshered.stan",
#                      data = list(N = 100,
#                                  ncov = 1000,
#                                  y = y.1000[,1],
#                                  X = cov.1000,
#                                  sigma_indic = 10,
#                                  mu_indic = 0,
#                                  ngroup = ngroup_n),
#                      control = list(adapt_delta = 0.9999, max_treedepth = 15),
#                      chains = 4, cores = 8, iter = 2000)
# 
# median.pars.cass_p1000 <- summary(cassfit_p1000,pars = "beta")$summary[,6]
# 
# 
# newdf_p1000 <- data.frame(actual = rep(beta.1000,4),
#                          fitted = c(median.pars.hs_p1000,
#                                     coef(lmfit_p1000)[-1],
#                                     coef(glmfit_p1000)[-1],
#                                     median.pars.cass_p1000))
# 
# library('reshape2')
# newdf_p1000 <- data.frame(melt(newdf_p1000),class = rep(rep(c("Horseshoe","OLS","LASSO","LN-CASS"),each = 1000),2))
# 
# ggplot(data = newdf_p1000, 
#        aes(x = rep(1:1000,8),y = value,col = variable)) + 
#   geom_point(size = 1) + 
#   geom_line() + 
#   facet_wrap(~ class,nrow = 1,scales = "free_y") + 
#   scale_color_manual(labels=c("Actual", "Fitted"),
#                      values = c("blue","red")) +
#   labs(x = "Coefficient index", 
#        y = "Value of coefficient",
#        color = "")
# 
# 
# library(pROC)
# 
# AUC_CASS_p1000 <- multiclass.roc(beta.1000, median.pars.cass_p1000)$auc
# 
# MAE_CASS_p1000 <- mean(((beta.1000-median.pars.cass_p1000)^2)^0.5)
# 
# AUC_HS_p1000 <- multiclass.roc(beta.1000,median.pars.hs_p1000)$auc
# 
# MAE_HS_p1000 <- mean(((beta.1000-median.pars.hs_p1000)^2)^0.5)
# 
# AUC_OLS_p1000 <- multiclass.roc(beta.1000,coef(lmfit_p1000)[-1])$auc
# 
# MAE_OLS_p1000 <- mean(((na.exclude(beta.1000-coef(lmfit_p1000)[-1]))^2)^0.5)
# 
# AUC_SGL_p1000 <- multiclass.roc(beta.1000, grplassofit_p1000$fit$beta[,8])$auc
# 
# MAE_SGL_p1000 <- mean(((beta.1000-grplassofit_p50$fit$beta[,8])^2)^0.5)
# 
# AUC_LASSO_p1000 <- multiclass.roc(beta.1000,coef(glmfit_p1000)[-1])$auc
# 
# MAE_LASSO_p1000 <- mean(((beta.1000-coef(glmfit_p1000)[-1])^2)^0.5)


## n = 100, p = 250 ##

int.250 <- -0.4

beta.250 <- c(rep(0,30),rep(2,10),rep(0,20),rep(2,10),(rep(-1,10) + rnorm(10,mean = 0,sd = 0.1)),rep(0,15),-0.5,-0.5,0,-0.5,-0.5,rep(0,25),rep(1,5),rep(0,25),(rep(2,5) + rnorm(5,mean = 0,sd = 0.1)), rep(0,25),-1,-1,-1,-2,-1,0.5,.5,.5,2,.5,rep(0,25),rep(3,15),rep(0,15))

cov.250 <- randomLHS(100,250)

y.250 <- int.250 + cov.250 %*% beta.250 + rnorm(100,mean = 0,sd = 0.5)

#SGL#
grplassofit_p250 <- cvSGL(data = list(x = cov.250,y = y.250),
                          index = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10),rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,10),rep(16,10),rep(17,10),rep(18,10),rep(19,10),rep(20,10),rep(21,10),rep(22,10),rep(23,10),rep(24,10),rep(25,10)),
                          maxit = 100000)


library(rstanarm)
hsfit_p250 <- stan_glm(y.250 ~ . ,
                       data = data.frame(cov.250,y.250),
                       prior = hs(df = 2,global_scale = 0.2),
                       family = gaussian(),
                       prior_intercept = normal(0,250),
                       adapt_delta = 0.9999,
                       chains = 5)

median.pars.hs_p250 <- apply(as.matrix(hsfit_p250),2,median)[2:(length(beta.250)+1)]



library(glmnet)
glmfit_p250 <- cv.glmnet(x = cov.250,y = y.250)
library(rstan)

N_obs <- 100
ncov <- length(cov.250[1,])
ngroup_n <- ncov/5

cassfit_p250 <- stan("simcasshered.stan",
                     data = list(N = 100,
                                 ncov = 250,
                                 y = y.250[,1],
                                 X = cov.250,
                                 sigma_indic = 10,
                                 mu_indic = 0,
                                 ngroup = ngroup_n),
                     control = list(adapt_delta = 0.99, max_treedepth = 16),
                     chains = 4, cores = 4, iter = 2000)

median.pars.cass_p250 <- summary(cassfit_p250,pars = "beta")$summary[,6]


newdf_p250 <- data.frame(actual = rep(beta.250,4),
                         fitted = c(grplassofit_p250$fit$beta[,8],
                                    median.pars.hs_p250,
                                    coef(glmfit_p250)[-1],
                                    median.pars.cass_p250))

library('reshape2')
newdf_p250 <- data.frame(melt(newdf_p250),class = rep(rep(c("SGL","Horseshoe","LASSO","LN-CASS"),each = 250),2))

ggplot(data = newdf_p250, 
       aes(x = rep(1:250,8),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "") +
  theme(text = element_text(size = 20))


library(pROC)

AUC_CASS_p250<- multiclass.roc(beta.250, median.pars.cass_p250)$auc

MAE_CASS_p250 <- mean(((beta.250-median.pars.cass_p250)^2)^0.5)

AUC_HS_p250 <- multiclass.roc(beta.250,median.pars.hs_p250)$auc

MAE_HS_p250 <- mean(((beta.250-median.pars.hs_p250)^2)^0.5)

AUC_SGL_p250 <- multiclass.roc(beta.250, grplassofit_p250$fit$beta[,8])$auc

MAE_SGL_p250 <- mean(((beta.250-grplassofit_p250$fit$beta[,8])^2)^0.5)

AUC_LASSO_p250 <- multiclass.roc(beta.250,coef(glmfit_p250)[-1])$auc

MAE_LASSO_p250 <- mean(((beta.250-coef(glmfit_p250)[-1])^2)^0.5)

## n = 100, p = 500 ##

int.500 <- 0.9

beta.500 <- c(rep(0,30),rep(2,10),rep(0,20),rep(2,10),(rep(-1,10) + rnorm(10,mean = 0,sd = 0.1)),rep(0,15),-0.5,-0.5,0,-0.5,-0.5,rep(0,25),rep(1,5),rep(0,25),(rep(2,5) + rnorm(5,mean = 0,sd = 0.1)), rep(0,25),-1,-1,-1,-2,-1,0.5,.5,.5,2,.5,rep(0,25),rep(3,15),rep(0,15),rep(0,30),rep(2,10),rep(0,20),rep(2,10),(rep(-1,10) + rnorm(10,mean = 0,sd = 0.1)),rep(0,15),-0.5,-0.5,0,-0.5,-0.5,rep(0,25),rep(1,5),rep(0,25),(rep(2,5) + rnorm(5,mean = 0,sd = 0.1)), rep(0,25),-1,-1,-1,-2,-1,0.5,.5,.5,2,.5,rep(0,25),rep(3,15),rep(0,15))

cov.500 <- randomLHS(100,500)

y.500 <- int.500 + cov.500 %*% beta.500 + rnorm(100,mean = 0,sd = 0.5)

#SGL#
grplassofit_p500 <- cvSGL(data = list(x = cov.500,y = y.500),
                          index = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10),rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,10),rep(16,10),rep(17,10),rep(18,10),rep(19,10),rep(20,10),rep(21,10),rep(22,10),rep(23,10),rep(24,10),rep(25,10),rep(26,10),rep(27,10),rep(28,10),rep(29,10),rep(30,10),rep(31,10),rep(32,10),rep(33,10),rep(34,10),rep(35,10),rep(36,10),rep(37,10),rep(38,10),rep(39,10),rep(40,10),rep(41,10),rep(42,10),rep(43,10),rep(44,10),rep(45,10),rep(46,10),rep(47,10),rep(48,10),rep(49,10),rep(50,10)),
                          maxit = 100000)


library(rstanarm)
hsfit_p500 <- stan_glm(y.500 ~ . ,
                       data = data.frame(cov.500,y.500),
                       prior = hs(df = 2,global_scale = 0.2),
                       family = gaussian(),
                       prior_intercept = normal(0,500),
                       adapt_delta = 0.9999)

median.pars.hs_p500 <- apply(as.matrix(hsfit_p500),2,median)[2:(length(beta.500)+1)]



library(glmnet)
glmfit_p500 <- cv.glmnet(x = cov.500,y = y.500)
library(rstan)

N_obs <- 100
ncov <- length(cov.500[1,])
ngroup_n <- ncov/5

cassfit_p500 <- stan("simcasshered.stan",
                     data = list(N = 100,
                                 ncov = 500,
                                 y = y.500[,1],
                                 X = cov.500,
                                 sigma_indic = 10,
                                 mu_indic = 0,
                                 ngroup = ngroup_n),
                     control = list(adapt_delta = 0.99, max_treedepth = 16),
                     chains = 4, cores = 10, iter = 2000)

median.pars.cass_p500 <- summary(cassfit_p500,pars = "beta")$summary[,6]


newdf_p500 <- data.frame(actual = rep(beta.500,4),
                         fitted = c(grplassofit_p500$fit$beta[,8],
                                    median.pars.hs_p500,
                                    coef(glmfit_p500)[-1],
                                    median.pars.cass_p500))

library('reshape2')
newdf_p500 <- data.frame(melt(newdf_p500),class = rep(rep(c("SGL","Horseshoe","LASSO","LN-CASS"),each = 500),2))

ggplot(data = newdf_p500, 
       aes(x = rep(1:500,8),y = value,col = variable)) + 
  geom_point(size = 1) + 
  geom_line() + 
  facet_wrap(~ class,nrow = 1,scales = "free_y") + 
  scale_color_manual(labels=c("Actual", "Fitted"),
                     values = c("blue","red")) +
  labs(x = "Coefficient index", 
       y = "Value of coefficient",
       color = "") +
  theme(text = element_text(size = 20))


library(pROC)

AUC_CASS_p500<- multiclass.roc(beta.500, median.pars.cass_p500)$auc

MAE_CASS_p500 <- mean(((beta.500-median.pars.cass_p500)^2)^0.5)

AUC_HS_p500 <- multiclass.roc(beta.500,median.pars.hs_p500)$auc

MAE_HS_p500 <- mean(((beta.500-median.pars.hs_p500)^2)^0.5)

AUC_SGL_p500 <- multiclass.roc(beta.500, grplassofit_p500$fit$beta[,8])$auc

MAE_SGL_p500 <- mean(((beta.500-grplassofit_p500$fit$beta[,8])^2)^0.5)

AUC_LASSO_p500 <- multiclass.roc(beta.500,coef(glmfit_p500)[-1])$auc

MAE_LASSO_p500 <- mean(((beta.500-coef(glmfit_p500)[-1])^2)^0.5)