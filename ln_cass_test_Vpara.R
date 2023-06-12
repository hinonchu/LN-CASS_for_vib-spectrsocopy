library(cepp)
library(dplyr)
library(data.table)
library(pROC)
library(rstan)
library(plotly)
library(purrr)
library(float)
library(foreach)
library(doParallel)

# FRI.txt contains the Raman intensity data, while rs.txt contains the associated wavenumber for the Raman shift.
# person_type contains the classification for each data sample.

raman_dat_all <- as.matrix(read.table('FRI.txt'))
raman_data_transposed <- t(raman_dat_all)
raman_dat_type <- as.matrix(read.table('person_type.txt'))

rs <- read.table('rs.txt')


alondat <- raman_data_transposed
alondat <- (alondat/(max(alondat)))
alony <- as.numeric(raman_dat_type)
alonnames <- as.matrix(rs)

alondatlognorm <- apply(alondat,2,log)
alondatlognorm <- apply(alondatlognorm,2,function(x) (x - mean(x))/sd(x))

pvals <- rep(NA,length(alondat[1,]))
zscores <- rep(NA,length(pvals))


cores=detectCores()
c1 <- makeCluster(cores[1]-1)
registerDoParallel(c1)


for (i in 1:length(pvals)) {
  X.i <- alondatlognorm[,i]

  }  
for (i in 1:length(pvals)){
  my.summary <- summary(glm(y ~ X.i + 0,
                            data = data.frame(X.i,y = alony),
                            family = "gaussian"))$coefficients
  
  zscores[i] <- my.summary[3]
  pvals[i] <- my.summary[4]
}

stopCluster(c1)

alonX <- alondatlognorm[,order(pvals)[1:100]]
alonnames.pared <- alonnames[order(pvals)[1:100]]

registerDoParallel(c1)

inv_logit <- function(x) 1/(1 + exp(-x))
mean.pred <- rep(NA,length(alony))
pred <- rep(NA,length(alony))

aloncass <- stan("CASS_RAMAN.stan",
                 data = list(N = length(alony),
                             ncov = 100,
                             y = alony,
                             x = alonX,
                             sigma_indic = 10,
                             mu_indic = 0,
                             tau = 5),
                 chains = 4,
                 cores = 8,
                 iter = 1000,
                 control = list(adapt_delta = 0.99, max_treedepth = 15))

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

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

for (i in 1:length(alony)){
  
  newX <- alonX[-c(i),]
  newy <- alony[i]
}



  
for (i in 1:length(alony)){
  loolist[[i]] <- stan("CASS_RAMAN.stan",
                       data = list(N = length(alony),
                                   ncov = 100,
                                   y = alony,
                                   x = alonX,
                                   sigma_indic = 10,
                                   mu_indic = 0,
                                   tau = 5),
                       cores = 10,
                       chains = 1,
                       iter = 500,
                       control = list(adapt_delta = 0.99, max_treedepth = 15))
  #if (i == 1) print("Started!", Sys.time()+Sys.Date())
  #if (i == 58) print("Done!", Sys.time()+Sys.Date())
}

## posterior predictions


for (i in 1:length(alony)){
  
  beta.matrix <- extract(loolist[[i]],pars = "beta")$beta
  
  X.loo <- alonX[i,]
  
  for (j in 1:100){
    pred[j] <- rbinom(1,1,inv_logit(X.loo%*%beta.matrix[j,]))
  }
  
  mean.pred[i] <- mean(pred)
}


alon.auc.cass <- auc(roc(alony,mean.pred))

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

fig <- plot_ly(x=(1:length(alony)),y=mean.pred,color=alony,colors=pal,showlegend=TRUE)
fig <- fig %>% layout(xaxis = xax, yaxis = yax, title="Mean Prediction Probability of Sex Being Male")
#fig <- fig %>% hide_colorbar()
fig


zscores_dat <- data.frame(zscores)

zscore_fig <- plot_ly(x=alonnames[,1],y=zscores_dat)
zscore_fig

