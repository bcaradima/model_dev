# Fit and cross-validate the individual models,
# deploy the joint models using deploy.jsdm()

# > SELECT influence factors ####
# Predictors for paper 1 (BDM species)
K <- c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD")
predictors <- prepare.inputs(K, sample.bdms, center=TRUE)
write.csv(predictors, "outputs/predictors.csv", row.names = F)

# Predictors for variable selection (paper 2, all programs)
predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.bdmf, center = TRUE)
write.csv(predictors.vs, 'outputs/predictors_vs_bdmf.csv', row.names = F)

predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.bdms, center = TRUE)
write.csv(predictors.vs, 'outputs/predictors_vs_bdms.csv', row.names = F)

predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.invf, center = TRUE)
write.csv(predictors.vs, 'outputs/predictors_vs_invf.csv', row.names = F)

predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.invf.plat, center = TRUE)
write.csv(predictors.vs, 'outputs/predictors_vs_invfp.csv', row.names = F)



# > Individual model calibration ####
isdm <- run.isdm(sample.bdms, predictors, trace = FALSE)

# How many models lowered the deviance?
table(sapply(isdm$models, function(m){
  m$deviance < m$null.deviance
}))

# How many taxa used optim?
table(sapply(isdm$models, function(m){
  m$optim
}))

# Extract beta/dev/prob ####
isdm.beta <- isdm$parameters
isdm.dev <- isdm$deviance
isdm.prob <- isdm$probability

### > k-fold cross-validation ####
# Manually calibrate (i.e., train) iSDMs with run.isdm(), then run each k-fold calibration
# with cv.isdm() to obtain predictions
# Train iSDMs ####
bdm.1 <- run.isdm(train1, predictors, trace = FALSE)
bdm.2 <- run.isdm(train2, predictors, trace = FALSE)
bdm.3 <- run.isdm(train3, predictors, trace = FALSE)

table(sapply(bdm.1$models, function(m){
  m$deviance < m$null.deviance
}))

table(sapply(bdm.2$models, function(m){
  m$deviance < m$null.deviance
}))

table(sapply(bdm.3$models, function(m){
  m$deviance < m$null.deviance
}))

bdm.1$probability$Fold <- 1
bdm.2$probability$Fold <- 2
bdm.3$probability$Fold <- 3

isdm.train <- bind_rows(bdm.1$probability, bdm.2$probability, bdm.3$probability)
isdm.train$Type <- "Training"

### Test iSDMs ####
# Note that in testing iSDMs, predictions are made on NA observations
bdm.1.test <- cv.isdm(bdm.1, test1, predictors); bdm.1.test$Fold <- 1
bdm.2.test <- cv.isdm(bdm.2, test2, predictors); bdm.2.test$Fold <- 2
bdm.3.test <- cv.isdm(bdm.3, test3, predictors); bdm.3.test$Fold <- 3
isdm.test <- bind_rows(bdm.1.test, bdm.2.test, bdm.3.test)
isdm.test$Type <- "Testing"
rm(bdm.1.test, bdm.2.test, bdm.3.test)


### Combine train/test results ####
isdm.test <- isdm.test[, colnames(isdm.train), with=FALSE]
isdm.cv.prob <- bind_rows(isdm.train, isdm.test)
rm(isdm.train, isdm.test)

### Deviance from list of list of GLMs ####
# Assumes bdm.1, bdm.2, bdm.3 are trained run.isdm() objects that exist in environment
isdm.cv.dev <- data.table()
for (k in 1:3){
  isdm <- get(paste("bdm.",k,sep=""))
  
  # TRAINING
  null.deviance.train <- isdm$deviance$null.dev
  residual.deviance.train <- isdm$deviance$res.dev
  d.train <- isdm$deviance$D2
  n.samples <- isdm$deviance$n.samples
  n.present <- isdm$taxa
  std.deviance.train <- residual.deviance.train/n.samples
  
  # Tidy data for output
  train.deviance <- data.table(Taxon = names(isdm$taxa),
                               # Statistics
                               null.deviance = null.deviance.train,
                               residual.deviance = residual.deviance.train,
                               std.deviance = std.deviance.train,
                               D2 = d.train,
                               # Sample information
                               n.samples = n.samples,
                               n = n.present,
                               # Validation information
                               Type = "Training",
                               Fold = k,
                               Model = "iSDM",
                               stringsAsFactors = F)
  
  train.deviance <- train.deviance[, c("Taxon", "Type", "Fold", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n")]
  
  # bind_rows() tidy training data to output
  isdm.cv.dev <- bind_rows(isdm.cv.dev, train.deviance)
  rm(train.deviance, null.deviance.train, residual.deviance.train, std.deviance.train, d.train, n.samples, n.present)
  
  # TESTING
  # Prepare input/output data
  taxa <- names(isdm$taxa) # get names of taxa found in training
  test.data <- get(paste("test", k, sep = "")) # get test data corresponding to k-fold
  test.data <- setDT(test.data)
  test.data <- select_(test.data, .dots=c("SiteId", "SampId", taxa)) # only keep taxa found in training
  inf.fact <- isdm$inf.fact
  
  # Calculate intercepts from test data
  n.present <- occur.freq(test.data)
  n.total <- apply(select(test.data, -SiteId, -SampId), 2, function(j){
    sum(!is.na(j), na.rm = TRUE)
  })
  freq <- n.present/n.total
  freq <- ifelse(freq==0, 1e-04, freq)
  alpha.taxa <- -log(1/freq-1)
  
  # get predictions to calculate testing deviance for each taxa
  for (j in 1:length(taxa)){
    taxon <- taxa[j]
    m <- isdm$models[[j]]
    
    # prepare parameters for null model
    null.parameters <- c(alpha.taxa[taxon], rep(0, length(inf.fact)))
    
    # get the observations
    obs <- select_(test.data, .dots = c("SiteId", "SampId", taxon))
    colnames(obs) <- c("SiteId", "SampId", "Obs")
    obs <- na.omit(obs)
    
    # join the inputs to the observations
    test.model.data <- left_join(obs, predictors, by = c("SiteId", "SampId"))
    test.model.data <- na.omit(test.model.data)
    
    y.test <- test.model.data[, "Obs"]
    env.cond.test <- as.matrix(test.model.data[, inf.fact])
    
    # null deviance
    z <- null.parameters[1] + env.cond.test%*%null.parameters[-1]
    null.prob.test <- 1/(1+exp(-z))
    null.deviance.taxon.test <- sum(-2*log(ifelse(y.test==1, null.prob.test, 1-null.prob.test)), na.rm = TRUE)
    # null.deviance.taxon.test <- sum(-2*(y.test*log(null.prob.test)+(1-y.test)*log(1-null.prob.test)),na.rm=TRUE)
    
    # residual deviance
    # calculate predicted probabilities
    if (m$optim){
      z <- m$par["Intercept"] + env.cond.test%*%m$par[-1]
      z <- as.vector(z)
      residual.prob.test <- 1/(1+exp(-z))
    }else{
      residual.prob.test <- predict.glm(m, newdata = test.model.data, type = "response") # get model predictions for test sites
    }
    
    residual.deviance.taxon.test <- sum(-2*log(ifelse(y.test == 1, residual.prob.test, 1-residual.prob.test)), na.rm = TRUE)
    
    dt <- data.table(Taxon = taxon,
                     # Statistics
                     null.deviance = null.deviance.taxon.test,
                     residual.deviance = residual.deviance.taxon.test,
                     std.deviance = residual.deviance.taxon.test/length(y.test),
                     D2 = (null.deviance.taxon.test - residual.deviance.taxon.test)/null.deviance.taxon.test,
                     # Sample information
                     n.samples = length(y.test),
                     n = n.present[taxon],
                     # Validation information
                     Type = "Testing",
                     Fold = k,
                     Model = "iSDM",
                     stringsAsFactors = F)
    
    dt <- dt[, c("Taxon", "Type", "Fold", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n")]
    isdm.cv.dev <- bind_rows(isdm.cv.dev, dt)
  }
}

rm(bdm.1, bdm.2, bdm.3)

# > Deploy jSDM workspaces ####
K <- c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD")
# BDM species deployments for calibration to entire dataset (cv=0) and 3-fold cross-validation (cv 1-3)
deploy.jsdm(K, sample.bdms, "bdms", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "bdms_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "bdms_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "bdms_train3", center=T, cv=3)
