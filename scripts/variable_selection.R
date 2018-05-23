# Script outline:
# - obtain correlation structure of pre-selected variables
# - get models with independent variables
# - Loop: for each possible number of parameters...
# - doParallel: Distribute taxa to the CPU cores...
# - For each taxa, run all models
# - For each model, run k-fold cross-validation
# - Obtain relative deviance statistics
# - Store/save data in a tidy output for extraction

# Packages ####
library(data.table) # improved data.frame class
library(tibble) # improved data frames (needed for debugging)
library(dplyr)
library(tidyr)
library(gtools) # needed for ?combn
library(doParallel)

# Global parameters ####
# Set number of CPU cores to utilize (recommended: n.cores - 1)
cores <- 5
# Set number of parameters to explore (recommended: 4-10)
npar <- 4:13

# Read datasets ####
i <- fread('inputs/predictors_centered_all.csv', header = T, sep = ',', stringsAsFactors = F); setkey(i, SampId)
full.dataset <- fread('inputs/bdm.sample.csv', header = T, sep = ',', stringsAsFactors = F); setkey(full.dataset, SampId)
train1 <- fread('inputs/bdm.train1.csv', header = T, sep = ',', stringsAsFactors = F); setkey(train1, SampId)
train2 <- fread('inputs/bdm.train2.csv', header = T, sep = ',', stringsAsFactors = F); setkey(train2, SampId)
train3 <- fread('inputs/bdm.train3.csv', header = T, sep = ',', stringsAsFactors = F); setkey(train3, SampId)

test1 <- fread('inputs/bdm.test1.csv', header = T, sep = ',', stringsAsFactors = F); setkey(test1, SampId)
test2 <- fread('inputs/bdm.test2.csv', header = T, sep = ',', stringsAsFactors = F); setkey(test2, SampId)
test3 <- fread('inputs/bdm.test3.csv', header = T, sep = ',', stringsAsFactors = F); setkey(test3, SampId)

n <- apply(full.dataset[, !c("SiteId", "SampId"), with=F], 2, sum, na.rm = TRUE)
n <- sort(n, decreasing = TRUE)
inf.fact <- colnames(i[,!c("SiteId", "SampId"), with=F])

### Correlation structure ####
correlation.structure <- cor(na.omit(i[,!c("SiteId", "SampId"), with=F]))
correlation.structure <- as.data.frame(correlation.structure)
correlation.structure <- rownames_to_column(correlation.structure, var = "Var1")
correlation.structure <- gather(correlation.structure, Var2, Value, -Var1)
correlation.structure <- as.data.table(correlation.structure)
output <- list()

### LOOP: Parameters ####
system.time(for (np in 1:length(npar)) { # number of parameters
  n.parameters <- npar[np]
  
  ### PRE-SELECT MODELS ####
  # Generate all combinations of the elements of [k] influence factors taken [p] at a time
  all.possible.models <- t(combn(inf.fact, n.parameters, simplify = TRUE))
  
  # Filter all possible models by variable(s) to candidate models
  ind <- apply(all.possible.models, 1, function(model){
    a <- all(c("Temp", "IAR", "FV") %in% model)
    # b <- all(c("Temp", "IAR", "FV") %in% model)
    # c <- all(c("Temp", "IAR", "FV") %in% model)
    # d <- a|b|c
    # d
  })
  candidate.models <- all.possible.models[ind, ]
  
  # Check for collinearity
  # Mark pairwise correlations among variables per model
  models.to.run.boolean <- apply(candidate.models, 1, function(model){
    pairs <- t(combn(model, 2, simplify = TRUE))
    pairs <- data.table(pairs); colnames(pairs) <- c("Var1", "Var2")
    pairs <- left_join(pairs, correlation.structure, by = c("Var1", "Var2"))
    setDT(pairs)
    pairs <- pairs[!(Var1=="Temp" & Var2=="Temp2"),]
    
    correlation <- ifelse(pairs$Value > 0.6 | pairs$Value < -0.6, FALSE, TRUE)
    model.run <- ifelse(all(correlation)==TRUE, TRUE, FALSE)
  })
  
  # Subset candidate models to obtain models for estimation
  models.to.run <- candidate.models[models.to.run.boolean, ]
  rm(all.possible.models, ind, candidate.models, models.to.run.boolean)

  ### CROSS-VALIDATION ####
  output[[np]] <- list() # k-fold list per parameter
  
  folds <- 1:3
  for (k in folds){
    fold <- folds[k]
    # Get the k-fold calibration/prediction data
    # Only model taxa that occur in calibration data
    train.data <- as.data.table(get(paste("train", k, sep = "")))
    n.train <- apply(train.data[,!c("SiteId", "SampId"), with=F], 2, sum, na.rm = T)
    n.train <- sort(n.train, decreasing = TRUE)
    n.train <- n.train[n.train > 0]
    taxa <- names(n.train)
    
    train.data <- train.data[, c("SampId", taxa), with = F]
    test.data <- as.data.table(get(paste("test", k, sep = "")))
    
    ### Parallellized loop through taxa ####
    # Taxa are parallelized to CPU cores; models per taxa are executed on that core
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    dev.comm <- foreach(j = 1:length(taxa), .combine=function(...)rbindlist(list(...)), .inorder = FALSE, .multicombine = TRUE, .packages=c("data.table", "dplyr")) %dopar% {
    # for (j in 1:length(taxa)) {
      # Get taxon
      taxon <- taxa[j]
      
      # Select calibration/prediction observations for taxon
      train.taxon.obs <- train.data[, c("SampId", taxon), with = F]
      test.taxon.obs <- test.data[, c("SampId", taxon), with = F]

      # For each model, calibrate and predict
      dev.taxon <- apply(models.to.run, 1, function(x){
        # Select variables as needed
        x.model <- i[, c("SampId", x), with=F]

        ### Train/test data ####
        # Prepare training data
        train.model.data <- left_join(train.taxon.obs, x.model, by = "SampId")
        setDT(train.model.data)
        train.model.data <- na.omit(train.model.data) # Drop NAs
        # Get training response by sample; use index position
        y.train <- train.model.data[[2]] # return vector of integers
        
        # Prepare testing data
        test.model.data <- left_join(test.taxon.obs, x.model, by = "SampId")
        setDT(test.model.data)
        test.model.data <- na.omit(test.model.data) # Drop NAs
        # Get training response by sample; use index position
        y.test <- test.model.data[[2]]
        
        # Check if taxon still occurs after removing NAs
        if (sum(y.train, na.rm = TRUE) > 0){
          # Define initial parameters
          n.present <- sum(y.train, na.rm = TRUE) # number of occurrences (i.e., number of 1s)
          n.total <- sum(!is.na(y.train)) # number of samples (i.e., number of 0s and 1s)
          freq <- n.present/n.total
          alpha <- -log(1/freq-1)
          start.values <- c(alpha, rep(0, length(x)))
          
          ### glm() ####
          # Run the model with base "glm"
          f <- paste(taxon, "~", paste(x, collapse = "+"))
          f <- as.formula(f)
          m <- glm(f, family = binomial(link = "logit"), data = train.model.data, na.action = na.omit, start = start.values)
          m$data <- NULL
          m$optim <- FALSE
          
          ### optim() ####
          # If residual deviance is still greater than null deviance, use more robust optimizer
          if (m$deviance > m$null.deviance){
            rm(m)
            # Prepare response/input data
            env.cond <- as.matrix(train.model.data[, x, with=F])
            # Calculate null deviance based on initial parameters
            z.start <- start.values[1] + env.cond%*%start.values[-1]
            p.start <- 1/(1+exp(-z.start))
            null.deviance <- sum(-2*log(ifelse(y.train == 1, p.start, 1-p.start)), na.rm = TRUE)

            # Define function to optimize
            fn <- function(par, env.cond, response){
              z <- par[1] + env.cond %*% par[-1]
              p <- 1/(1+exp(-z))
              dev <- sum(-2*(response*log(p)+(1-response)*log(1-p)),na.rm=TRUE)
              return(dev)
            }
            
            # Maximize likelihood of observations given inputs and initial parameters
            m <- optim(start.values, fn, env.cond = env.cond, response = y.train)
            
            # Prepare glm-like object
            names(m$par) <- c("Intercept", x)
            m$y <- y.train
            m$z.est <- m$par[1] + env.cond%*%m$par[-1]
            m$fitted.values <- 1/(1+exp(-m$z.est))
            m$deviance <- sum(-2*log(ifelse(y.train == 1, m$fitted.values, 1-m$fitted.values)), na.rm = TRUE)
            m$null.deviance <- null.deviance
            m$optim <- TRUE
          }
          # Independent predicted probabilities
          if (m$optim){
            env.cond <- as.matrix(test.model.data[, x, with=F])
            z <- m$par["Intercept"] + env.cond%*%m$par[-1]
            prob <- 1/(1+exp(-z))
          }else{
            prob <- predict.glm(m, newdata = test.model.data, type = "response") # get model predictions for test sites
          }
          ### Deviance statistics ####
          # Training deviance
          relative.deviance.train <- m$deviance/length(y.train)
          
          # Testing deviance
          residual.deviance.test <- sum(-2*log(ifelse(y.test == 1, prob, 1-prob)), na.rm = TRUE)
          relative.deviance.test <- residual.deviance.test/length(y.test)
          
          # Deviance statistics by taxon
          dev.taxon <- c(relative.deviance.train, relative.deviance.test, length(y.train), length(y.test))
        } else{
          # Taxon has no occurrences in the data
          dev.taxon <- c(NA, NA, NA, NA)
        }
        dev.taxon
      })
      ### Prepare output ####
      dev.taxon <- data.table(Taxon = taxon, t(dev.taxon), stringsAsFactors = F)
      colnames(dev.taxon) <- c("Taxon", "rdev.train", "rdev.test", "length.y.train", "length.y.test")
      dev.taxon
    }
    # Exit parallel loop with output: dev.comm
    output[[np]][[paste("fold", k, sep = '')]] <- dev.comm
    stopCluster(cl)
  }
  output[[np]][["models.to.run"]] <- models.to.run
  output[[np]][["n.parameters"]] <- n.parameters
  cat("Models completed for", n.parameters, "parameters \n")
})

save.image("variable_selection.RData")