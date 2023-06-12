#load all the required libraries

library(dplyr)
library(data.table)
library(pROC)
library(rstan)
library(plotly)
library(purrr)
library(float)
library(foreach)
library(doParallel)

#define variables for Raman data in this case.Please change the file for IR data to do this with that type of data.

raman_dat_all <- as.matrix(read.csv('test_sugars_raman_n20.csv'))
raman_data_transposed <- t(raman_dat_all)
raman_dat_type <- as.matrix(read.csv('type_dat_n20.csv',header = FALSE))

rs <- read.table('rs.txt')

#rs for Raman Shift. Below adapts code from Thompson et al.

alondat <- t(raman_dat_all)
alondat <- (scale(alondat))
alondat <- alondat-min(alondat)
alony <- as.numeric(raman_dat_type)
alonnames <- as.matrix(rs)

alondatlognorm <- apply(alondat,2,log)
alondatlognorm <- apply(alondatlognorm,2,function(x) (x - mean(x))/sd(x))
alondatlognorm[is.nan(alondatlognorm)] <- 0

pvals <- rep(NA,35)
zscores <- rep(NA,35)


cores=detectCores()
c1 <- makeCluster(cores[1]-1)
registerDoParallel(c1)


for (i in 1:35) {
  X.i <- alondatlognorm[,i]
  
}  
for (i in 1:35){
  my.summary <- summary(glm(y ~ X.i + 0,
                            data = data.frame(X.i,y = alony),
                            family = "gaussian"))$coefficients
  
  zscores[i] <- my.summary[3]
  pvals[i] <- my.summary[4]
}

#stopCluster(c1)

alonX <- alondatlognorm[,order(pvals)[1:35]]
alonnames.pared <- alonnames[order(pvals)[1:35]]

registerDoParallel(c1)

inv_logit <- function(x) 1/(1 + exp(-x))
mean.pred <- rep(NA,length(alony))
pred <- rep(NA,length(alony))

aloncass <- stan("CASS_RAMAN2.stan",
                 data = list(N = 20,
                             ncov = 35,
                             y = alony,
                             X = alonX,
                             sigma_indic = 10,
                             mu_indic = 0,
                             tau = 5),
                 chains = 4,
                 cores = 8,
                 iter = 1000,
                 control = list(adapt_delta = 0.99, max_treedepth = 11))

## LOOCV

indices.pos <- which(alony == 1)
indices.neg <- which(alony == 0)

to.remove <- rep(NA,length(alony))

for (i in 1:length(alony)){
  if (alony[i] == 1){
    to.remove[i] <- sample(indices.neg,1)
  } else{
    to.remove[i] <- sample(indices.pos,1)
  }
}


loolist <- list()

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

for (i in 1:length(alony)){
  
  newX <- alonX[-c(i),]
  newy <- alony[i]
}




for (i in 1:length(alony)){
  loolist[[i]] <- stan("CASS_RAMAN2.stan",
                       data = list(N = length(alony),
                                   ncov = 35,
                                   y = alony,
                                   X = alonX,
                                   sigma_indic = 10,
                                   mu_indic = 0,
                                   tau = 5),
                       cores = 8,
                       chains = 1,
                       iter = 500,
                       control = list(adapt_delta = 0.99, max_treedepth = 12))

}

## loo cv for random forest

library(randomForest)

mean.pred.rf <- rep(NA,20)

for (i in 1:20){
  
  newX <- alonX[-c(i,to.remove[i]),]
  newy <- alony[-c(i,to.remove[i])]
  
  mean.pred.rf[i] <- predict(randomForest(newX,as.factor(newy),ntree = 1000),newdata = alonX[i,],type = "prob")[2]
  
}

## loo cv for lasso

library(glmnet)

mean.pred.lasso <- rep(NA,20)

for (i in 1:20){
  
  newX <- alonX[-c(i,to.remove[i]),]
  newy <- alony[-c(i,to.remove[i])]
  
  mean.pred.lasso[i] <- predict(cv.glmnet(newX,newy,family = "binomial"),newx = t(as.matrix(alonX[i,])),type = "response")
  
}




### Horseshoe

library(rstanarm)

for (i in 1:100){
  newx_hs <- alonX[,-c(i,to.remove[i])]
  newy_hs <- alony[-c(i,to.remove[i])]
  
hsfit <- stan_glm(alony ~ . ,
                  data = data.frame(newx_hs[,i],alony),
                  prior = hs(df = 2,global_scale = 0.2),
                  family = gaussian(),
                  prior_intercept = normal(0,60),
                  adapt_delta = 0.99)

median.pars.hs[i] <- apply(as.matrix(hsfit),2,median)[2:101]
}
lmfit <- lm(alony ~ . , data = data.frame(alondat,alony))

## compute AUCS
library(pROC)

alon.auc.cass <- roc(alony,mean.pred)
alon.auc.lasso <- auc(roc(alony,mean.pred.lasso))
alon.auc.rf <- auc(roc(alony,mean.pred.rf))



## posterior predictions


for (i in 1:20){
  
  beta.matrix <- extract(loolist[[i]],pars = "beta")$beta
  
  X.loo <- alonX[i,]
  
  for (j in 1:35){
    pred[j] <- rbinom(1,1,inv_logit(X.loo%*%beta.matrix[j,]))
  }
  
  mean.pred[i] <- mean(pred)
}


alon.auc.cass <- roc(alony,mean.pred)

pal <- c('red','blue')

f <- list(
  family = "Courier New, monospace",
  size = 16,
  color = "#7f7f7f"
)

xax <- list(
  title = "Sample",
  titlefont = f
)

yax <- list(
  title = "Mean Predicted Probability",
  titlefont = f
)

fig <- plot_ly(x=(1:length(alony)),y=mean.pred,color=alony,colors=pal,showlegend=FALSE)
fig <- fig %>% layout(xaxis = xax, yaxis = yax, title="Probability of Being Glucose (LN-CASS)")
fig <- fig %>% hide_colorbar()
fig

fig2 <- plot_ly(x=1:length(alony),y=mean.pred.lasso,color=alony,colors=pal,showlegend=TRUE)
fig2 <- fig2 %>% layout(xaxis = xax, yaxis = yax, title="Probability of Being Glucose (LASSO)")
fig2 <- fig2 %>% hide_colorbar()
fig2

fig3 <- plot_ly(x=1:length(alony),y=mean.pred.rf,color=alony,colors=pal,showlegend=TRUE)
fig3 <- fig3 %>% layout(xaxis = xax, yaxis = yax, title="Probability of Being Glucose (RF)")
fig3 <- fig3 %>% hide_colorbar()
fig3

fig4 <- plot_ly(x=1:length(alony),y=median.pars.hs,color=alony,colors=pal,showlegend=TRUE)
fig4 <- fig4 %>% layout(xaxis = xax, yaxis = yax, title="Probability of Being Glucose (HS)")
fig4 <- fig4 %>% hide_colorbar()
fig4

zscores_dat <- data.frame(zscores)

zscore_fig <- plot_ly(x=alonnames[,1],y=zscores_dat)
zscore_fig

mae_cass <- mean(((alony-mean.pred)^2)^0.5)
mae_lasso <- mean(((alony-mean.pred.lasso)^2)^0.5)
mae_rf<- mean(((alony-mean.pred.rf)^2)^0.5)
