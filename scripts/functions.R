### FUNCTIONS ####
# Description: this script defines all custom functions called in the other scripts. The functions are ordered as intuitively as possible to reflect the model development process. For the most part, these are considered to be stable functions.

# Utility functions ####
# Package installer and loader function: pass a vector of package names
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Get the occurrence frequency per taxon from a community matrix 
# Checks for site/sample IDs before calculation
occur.freq <- function(cm){
  if(!all(c("SiteId", "SampId") %in% colnames(cm))){
    n <- apply(cm, 2, sum, na.rm=TRUE)
  } else{
    n <- apply(select(cm, -SiteId, -SampId), 2, sum, na.rm=T)
  }
  n <- sort(n, decreasing = T)
  return(n)
}

# Print total NAs per column
check.na <- function(df){
  apply(df, 2, function(column){
    sum.na <- sum(is.na(column))
  })
}

# Calculate root mean squared error, given vector of residuals (y minus y hat)
rmse <- function(residuals){
  square.r <- (residuals)^2 # square the residuals
  mean.square.r <- sum(square.r)/length(residuals) # MSE: mean of the sum of squared residuals
  root.mean.square.error <- sqrt(mean.square.r) # RMSE: root mean squared error
  return(root.mean.square.error)
}

# Calculate the deviance for a taxon, given vectors of predicted probabilities and observations
# Note: jSDM must calculate deviance by matrix
deviance <- function(obs, pred){
  sum(-2*(obs*log(pred)+(1-obs)*log(1-pred)),na.rm=TRUE)
}

# Logistic link function outputs probabilities, given linear predictor
link.function <- function(z) {
  1 / (1 + exp(-z))
}

# Calculate D-squared statistic from glm() or optim()
D2 <- function(m){
  d.squared <- (m$null.deviance-m$deviance)/m$null.deviance
  return(d.squared)
}

# Returns expression labels given vector of influence factors
labeller <- function(inf.fact){
  labels <- c("A10m" = expression(paste("A10m (%)")),
              "Temp" = expression(paste("Temp (", degree, "C)")),
              "Temp2" = expression(paste("Temp"^2, " (", degree, "C"^2,")")),
              "FV" = expression(paste("FV (m/s)")),
              
              "F10m" = expression(paste("F10m (%)")),
              "F100m" = expression(paste("F100m (%)")),
              "FRI" = expression(paste("FRI (%)")),
              "bFRI" = expression(paste("bFRI (%)")),
              
              "IAR" = expression(paste("IAR (w"["c"]%*%"f"["c"],")")),
              "WV" = expression(paste("WV")),
              "Urban" = expression(paste("Urban (%)")),
              "LUD" = expression(paste("LUD (CE/km"^2,")"))
              
  )
  labels[inf.fact]
}

labeller.names <- function(inf.fact){
  labels <- c("A10m" = expression(paste("A10m")),
              "Temp" = expression(paste("Temp")),
              "Temp2" = expression(paste("Temp"^2)),
              "FV" = expression(paste("FV")),
              
              "F10m" = expression(paste("F10m")),
              "FRI" = expression(paste("FRI")),
              "bFRI" = expression(paste("bFRI")),
              
              "IAR" = expression(paste("IAR")),
              "Urban" = expression(paste("Urban")),
              "LUD" = expression(paste("LUD"))
              
  )
  labels[inf.fact]
}
# Returns expression label for linear model
labeller.lm <- function(x, y, df){
  f <- as.formula(paste(y, "~", x, sep=""))
  m <- lm(f, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

# Format species name (e.g., "Baetis_alpinus") as abbreviated genus and full species name (i.e., "B. alpinus")
labeller.species <- function(j){
  # taxon <- unlist(strsplit(j, "_"))
  taxon <- strsplit(j, "_")
  genus <- sapply(taxon, function(i){
    substring(i[[1]], 1, 1)
  })
  species <- sapply(taxon, function(i){
    i[2]
  })
  paste(genus,". ", species, sep="")
}

# Returns taxonomic level of a given taxon from BDMs
# Note: does not work for pooled taxonomic datasets...
taxon.rank <- function(taxon){
  cat('Identify: ', taxon, '\n')
  levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  hierarchy <- matrix(nrow=length(levels), ncol=length(levels))
  hierarchy <- lower.tri(hierarchy, diag = TRUE)
  row.names(hierarchy) <- levels
  colnames(hierarchy) <- levels
  rm(levels)
  
  taxon.rank.logical <- apply(inputs$taxonomy, 2, function(rank){
    taxon %in% rank
  })
  rank <- names(taxon.rank.logical[taxon.rank.logical]) # What is the highest taxonomic level of ID the taxon?
  
  # Get all taxa that belong to taxon at that rank 
  taxon.data <- dplyr::filter_(inputs$taxonomy, paste(rank, "=='",taxon,"'", sep=''))
  
  # Edit the taxonomy data and select one row,
  # unless rank is species (species always return one row)
  if(rank != "Species"){
    # Hierarchy matrix determines which taxonomic information to keep (TRUE) or discard (FALSE)
    # based on given taxonomic rank
    rank.hierarchy <- hierarchy[rank, ]
    rank.delete <- names(rank.hierarchy[!rank.hierarchy])
    
    taxon.data[, rank.delete] <- NA
    taxon.data <- taxon.data[1, ]
  }
  
  taxon.data$Rank <- rank
  taxon.data$Taxon <- taxon
  setDT(taxon.data)
  return(taxon.data) # Return list of data.tables for fast rbindlist()
}

# For a given sample, calls taxon.rank() on each taxa in the sample and returns complete taxonomy dataset
get.taxonomy <- function(sample.data){
  taxa <- colnames(sample.data)
  taxa <- taxa[!(taxa %in% c("SiteId", "SampId"))]
  sample.taxonomy <- lapply(taxa, function(j){
    taxon.rank(j)
  })
  sample.taxonomy <- rbindlist(sample.taxonomy)
  sample.taxonomy <- sample.taxonomy[, c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxon", "Rank")]
  return(sample.taxonomy)
}

# Copied from StackOverflow: remove NAs from data for nodes of a data.tree (for plotting taxonomic hierarchies)
# https://stackoverflow.com/questions/43803949/create-and-print-a-product-hierarchy-tree-without-na-from-data-frame-in-r-with
paste5 <- function(..., sep = " ", collapse = NULL, na.rm = F) {
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

# Run models ####
# prepare.inputs performs selection and transformation to prepare the influence factors 
# in a model, where [K] is a character vector of influence factors, [y] is the invertebrate data,
# [center] is a Boolean indicating whether the mean should be subtracted at each site
prepare.inputs <- function(K, y, center){
  # Get sites/samples of observation data
  site.samples <- select(y, SiteId, SampId)
  
  # Keep predictors at program sites
  dt <- p[, c("SiteId", "SampId", setdiff(K, "Temp2"))] # remove Temp2 if included in K
  
  dt <- left_join(site.samples, dt, by=c("SiteId", "SampId"))
  
  # If center argument TRUE, subtract the mean
  if(center){
    dt[, !(colnames(dt) %in% c("SiteId", "SampId"))] <- apply(dt[, !(colnames(dt) %in% c("SiteId", "SampId"))], 2, function(k){
      k - mean(k, na.rm = TRUE)
    })
  }
  
  # If temperature is an input, add a quadratic transformation
  if ("Temp" %in% K & !("Temp2" %in% colnames(dt))){
    dt$Temp2 <- dt$Temp^2
  }
  # # If temperature is an input, add a quadratic transformation
  # if ("FV" %in% K & !("FV2" %in% colnames(dt))){
  #   dt$FV2 <- dt$FV^2
  # }
  # Print number of NAs
  for (k in seq_along(K)){
    variable <- K[k]
    cat(variable, ": ", sum(is.na(dt[,variable])), " NAs \n")
  }
  return(dt)
}

# deploy.jsdm() creates the directories and data needed to deploy a joint model, where [K, y, and center]
# are arguments for prepare.inputs(), [trial] is the data.frame of observation data for sampling, and [cv]
# is an integer for controlling the sampling: 0 calibrates the model using the entire dataset, 1-3 writes
# 3-fold cross-validation training and testing data to the /input directory. Deploy.jsdm() also copies the 
# jsdm_inference.R script from model_dev/scripts to the top directory
deploy.jsdm <- function(K, y, trial, center, cv){
  # Create directory structure
  dir.top <- paste('jsdm_',trial,sep='')
  dir.input <- paste('jsdm_',trial,'/inputs',sep='')
  dir.create(dir.top)
  dir.create(dir.input)
  cat("Create directories -> \n")
  
  if(cv==0){
    # Write invertebrate data
    write.csv(y, paste(dir.input,'/community.csv',sep=''), row.names = F)
    cat("community samples: ", nrow(y), "\n")
  }else{
    set.seed(2017)
    y.sample <- y[sample(nrow(y)),]
    k.folds <- cut(seq(1, nrow(y.sample)), breaks=3, labels=FALSE)
    fold1 <- y.sample[which(k.folds == 1),]
    fold2 <- y.sample[which(k.folds == 2),]
    fold3 <- y.sample[which(k.folds == 3),]
    
    kcv <- matrix(c(1,2 , 2,3 , 1,3), nrow=2, ncol=3)
    
    # Training data: bind folds by row
    train.data <- bind_rows(get(paste("fold", kcv[1,cv], sep="")), get(paste("fold", kcv[2,cv], sep="")))
    write.csv(train.data, paste(dir.input,'/train', cv,'.csv', sep=''), row.names=F)
    cat("train", cv," samples: ", nrow(train.data), "\n")
    # Testing data
    test.data <- get(paste("fold", setdiff(1:3, c(kcv[1,cv], kcv[2,cv])), sep=""))
    write.csv(test.data, paste(dir.input,'/test', cv,'.csv', sep=''), row.names=F)
    cat("test", cv," samples: ", nrow(test.data), "\n")
  }
  
  # Write predictor data
  cat("NAs per predictor:\n")
  dt <- prepare.inputs(K, y, center)
  write.csv(dt, paste(dir.input,'/predictors.csv',sep=''), row.names = F)
  cat("predictor samples: ", nrow(dt)," -> ")
  
  # Copy script
  file.copy('scripts/jsdm_inference.R', dir.top)
  cat("Done\n")
}


# Calculate GLM likelihood given model inputs, observations, and parameters
likelihood <- function(par, env.cond, y){
  z <- par[1] + env.cond %*% par[-1]
  p <- 1/(1+exp(-z))
  dev <- sum(-2*(y*log(p)+(1-y)*log(1-p)),na.rm=TRUE)
  return(dev)
}

# Fits an individual generalized linear model (GLM) per taxon given site-species and predictor data,
# outputs a list storing the models and tidy data of (parameters, deviance, probabilities)
run.isdm <- function(cm, predictors, trace) {
  cm <- as.data.frame(cm)
  predictors <- as.data.frame(predictors)
  inf.fact <- colnames(predictors)[!colnames(predictors) %in% c("SiteId", "SampId")]
  
  # Left join observation sites/samples with inputs
  cm.sites <- select(cm, SiteId, SampId)
  
  if (length(cm.sites$SampId[duplicated(cm.sites$SampId)]) > 0){
    cat("WARNING:", length(cm.sites$SampId[duplicated(cm.sites$SampId)]), "duplicate sample IDs found!")
  }
  
  model.inputs <- left_join(cm.sites, predictors, by = c("SiteId", "SampId"))
  rm(cm.sites)
  # exclude rows from observations/inputs that do not have complete influence factors:
  # Logical vector of whether any influence factors are NA per sample
  ind <- !apply(is.na(model.inputs[,inf.fact]),1,FUN=any)
  ind <- ifelse(is.na(ind),FALSE,ind)
  
  cm <- cm[ind, ]
  model.inputs <- model.inputs[ind, ]
  
  # Count and filter taxa by cutoff
  n <- occur.freq(cm)
  n <- n[n > 0 & n < nrow(cm)]
  taxa <- names(n)
  
  # Prepare data for model input
  mdata <- cbind(cm, model.inputs)
  # Prepare the influence factors for the formulas
  variables <- paste(inf.fact, collapse = "+")
  
  fm <- lapply(taxa, function(j) {
    print.formula <- paste(j, "~", variables)
    f <- as.formula(print.formula)
    
    # Define initial parameters
    n.present <- sum(mdata[, j], na.rm=TRUE)
    n.total <- sum(!is.na(mdata[, j]))
    
    freq <- n.present/n.total
    alpha <- -log(1/freq-1)
    
    start.values <- c(alpha, rep(0, length(inf.fact)))
    
    m <- glm(f, family = binomial(link = "logit"), data = mdata, na.action = na.omit, control = glm.control(trace = trace), start = start.values)
    
    m$data <- NULL
    m$optim <- FALSE
    
    # If residual deviance is greater than null deviance, use a more robust optimizer
    if (m$deviance > m$null.deviance){
      rm(m)
      # Subset input/response data to omit NA rows in response
      ind <- !is.na(mdata[, j])
      y <- mdata[ind, j]
      env.cond <- as.matrix(mdata[ind, inf.fact])
      rm(ind)
      
      # Calculate null deviance based on initial parameters
      z.start <- start.values[1] + env.cond%*%start.values[-1]
      p.start <- 1/(1+exp(-z.start))
      null.deviance <- sum(-2*(y*log(p.start)+(1-y)*log(1-p.start)),na.rm=TRUE)
      
      # Maximize likelihood of observations given inputs and initial parameters
      m <- optim(start.values, likelihood, env.cond = env.cond, y = y)
      
      # Prepare glm-like object
      est <- m$par
      names(m$par) <- c("Intercept", inf.fact)
      m$y <- y
      m$z.est <- est[1] + env.cond%*%est[-1]
      m$fitted.values <- 1/(1+exp(-m$z.est))
      m$deviance <- sum(-2*(y*log(m$fitted.values)+(1-y)*log(1-m$fitted.values)),na.rm=TRUE)
      m$null.deviance <- null.deviance
      m$optim <- TRUE
    }
    cat("iSDM:", print.formula, "\n")
    return(m)
  })
  
  # Get parameters
  parameters <- sapply(fm, function(m) {
    if (m$optim){ # Process manually optimized GLMs
      b <- m$par[-1]
    }
    else{ # Process 'glm' objects
      coef <- coefficients(summary(m))[, "Estimate"]
      b <- coef[-1] # drop the intercept
    }
  })
  
  parameters <- as.data.frame(parameters)
  colnames(parameters) <- taxa
  parameters <- rownames_to_column(parameters, var = "Variable")
  
  # Spread the estimates to column by predictor
  parameters <- gather(parameters, Taxa, Estimate, -Variable)
  colnames(parameters) <- c("Variable", "Taxon", "Parameter")
  parameters$Model <- "iSDM"
  parameters <- parameters[, c("Taxon", "Variable", "Parameter", "Model")]
  setDT(parameters)
  
  # Get deviance statistics
  deviance <- lapply(taxa, function(j){
    m <- fm[[which(taxa == j)]]
    
    # Number of samples
    n.samples <- length(na.omit(mdata[, j]))
    
    null.dev <- m$null.deviance
    res.dev <- m$deviance
    std.res.dev <- m$deviance/n.samples
    D2 <- (null.dev-res.dev)/null.dev
    data.table(Taxon = j, null.dev = null.dev, res.dev = res.dev, std.res.dev = std.res.dev, n.samples = n.samples, n = n[j], D2 = D2, Model = "iSDM")
  })
  deviance <- rbindlist(deviance)
  
  # Get probabilities
  probability <- lapply(taxa, function(j){
    m <- fm[[which(taxa == j)]]
    # cat("Taxon:", j, "\n")
    prob <- m$fitted.values
    
    dt <- mdata[, c("SiteId", "SampId", j)]
    colnames(dt) <- c("SiteId", "SampId", "Obs")
    dt <- na.omit(dt)
    dt$Taxon <- j
    dt$Pred <- m$fitted.values
    dt$Model <- "iSDM"
    dt <- dt[, c("SiteId", "SampId", "Taxon", "Obs", "Pred", "Model")]
    return(dt)
  })
  probability <- rbindlist(probability)
  
  opt <- list("mdata" = mdata, "community" = cm, "predictors" = model.inputs, "inf.fact" = inf.fact, "taxa" = n, "models" = fm, "parameters" = parameters, "deviance" = deviance, "probability" = probability)
  return(opt)
}

# Select models ####
# Returns glm object of specified taxon from run.isdm() output 
select.isdm <- function(results, taxon) {
  x <- which(names(results$taxa) == taxon)
  m <- results$models[[x]]
  return(m)
}

# Locates and returns extracted Stan object from joint model results
select.jsdm <- function(results){
  ind <- sapply(results, function(i){
    # class(i)[1]=="list"
    i$jSDM==TRUE
  })
  ind <- unlist(ind)
  image <- results[[names(ind[ind==TRUE])]]
}


# Validation ####
# Note that independent predictions are made for samples with NA observations
cv.isdm <- function(results, data, predictors){
  # Create an output data.table
  pred.dt <- data.table()
  
  # Loop through the taxa
  taxa <- names(results$taxa)
  for (j in 1:length(taxa)){
    taxon <- taxa[j] # get taxon name
    m <- select.isdm(results, taxon)
    
    # get observations at test sites
    obs <- data[, c("SiteId", "SampId", taxon)]
    colnames(obs) <- c("SiteId", "SampId", "Obs")
    # obs <- na.omit(obs)
    
    # get the predictors for test sites
    test.data <- left_join(obs, predictors, by = c("SiteId", "SampId"))
    # test.data <- na.omit(test.data)
    
    # get the model predictions
    if(m$optim){
      env.cond <- as.matrix(test.data[, results$inf.fact])
      z <- m$par["Intercept"] + env.cond%*%m$par[-1]
      prob <- 1/(1+exp(-z))
    }else{
      prob <- predict.glm(m, newdata = test.data, type = "response") # get model predictions for test sites
    }
    # Combine predictions and observations for the test data
    dt <- data.frame(SiteId = test.data$SiteId, SampId = test.data$SampId, Pred = prob, Taxon = taxon, stringsAsFactors = F)
    dt <- left_join(dt, obs, by = c("SiteId", "SampId"))
    
    pred.dt <- bind_rows(pred.dt, dt)
  }
  pred.dt$Model <- "iSDM"
  return(pred.dt)
}

# Extracts tidy calibration (training) data (probability, parameters, deviance), then performs prediction (testing)
# against k-fold datasets. Warning: this function assumes a certain folder structure to read jSDM k-fold CV results
cv.jsdm <- function(directory, folders) {
  output = list("deviance" = data.table(), "probability" = data.table())
  
  # Loop through trials
  for (f in 1:length(folders)){
    folder <- folders[f]
    cat("Model run: ", folder, " | ")
    for (k in 1:3){
      model.image <- new.env()
      
      full.path <- paste(getwd(), directory, folder, sep="/")
      extensions <- identify.jsdm(full.path)
      load(paste(full.path, paste('_train', k, sep=""),"/", extensions$workspace, sep=""), envir = model.image)
      cat("Load image -> | ")
      
      # Extract objects from workspace image to local environment
      sites <- model.image$sites
      samples <- model.image$samples
      occur.taxa <- model.image$occur.taxa
      env.cond <- model.image$env.cond
      inf.fact <- model.image$inf.fact
      
      set.ident <- model.image$set.ident
      res <- model.image$res
      rm(model.image)
      
      ### Extract parameters
      res.extracted <- extract(res,permuted=TRUE,inc_warmup=FALSE)
      
      x <- as.matrix(env.cond[set.ident,inf.fact])
      colnames(x) <- inf.fact
      
      ind.maxpost <- which.max(res.extracted[["lp__"]])
      mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
      sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
      mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
      sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
      names(mu.beta.comm.maxpost) <- inf.fact
      names(sigma.beta.comm.maxpost) <- inf.fact
      alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,] # alpha at max posterior
      beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,] # "" ""
      rownames(beta.taxa.maxpost) <- inf.fact
      
      ### TRAIN prob/dev
      z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) + 
        x%*%beta.taxa.maxpost
      p.maxpost <- 1/(1+exp(-z))
      y.train <- occur.taxa
      n.train <- apply(y.train, 2, sum, na.rm = TRUE)
      n.obs.train <- apply(y.train, 2, function(j){sum(!is.na(j))})
      
      p.primitive.train <- apply(y.train, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
      p.primitive.train <- ifelse(p.primitive.train==0,1e-4,p.primitive.train)
      p.primitive.train <- matrix(rep(p.primitive.train,nrow(y.train)),nrow=nrow(y.train),byrow=TRUE)
      
      # Calculate deviance over the samples
      dev.resid.train <- sign(y.train-0.5)*sqrt(-2*(y.train*log(p.maxpost)+(1-y.train)*log(1-p.maxpost)))
      # Calculate the sum of squared residual deviance for one taxa
      # deviance.maxpost.train   <- sum(dev.resid.train^2, na.rm = T)
      # Equivalent to deviance calculated for iSDMs under cross-validation:
      # sum(-2*log(ifelse(isdm$mdata$Gammaridae == 1, m$fitted.values, 1-m$fitted.values)), na.rm = TRUE)
      
      # Calculate the null ("primitive") deviance
      deviance.primitive.train <- -2*sum(y.train*log(p.primitive.train)+(1-y.train)*log(1-p.primitive.train), na.rm = T)
      
      # Calculate the sum of squared residual deviance for all taxa
      deviance.maxpost.taxa.train <- apply(dev.resid.train^2,2,sum, na.rm=T) # NAs removed
      deviance.primitive.taxa.train <- -2*apply(y.train*log(p.primitive.train)+(1-y.train)*log(1-p.primitive.train),2,sum,na.rm=T)
      
      d.train <- 1-deviance.maxpost.taxa.train/deviance.primitive.taxa.train ###
      
      # Standardize deviance by number of observations in training data
      std.deviance.train <- deviance.maxpost.taxa.train/n.obs.train
      
      # tidy deviance data
      train.deviance <- data.table(null.deviance = deviance.primitive.taxa.train, 
                                   residual.deviance = deviance.maxpost.taxa.train, 
                                   std.deviance = std.deviance.train, 
                                   D2 = d.train,
                                   Type = "Training",
                                   Fold = k,
                                   Model = "jSDM",
                                   Trial = folder,
                                   stringsAsFactors = F)
      
      train.deviance$Taxon <- names(std.deviance.train)
      train.deviance$n.samples <- n.obs.train[train.deviance$Taxon]
      train.deviance$n <- n.train[train.deviance$Taxon]
      
      train.deviance <- train.deviance[, c("Taxon", "Type", "Fold", "Model", "Trial", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n")]
      
      # tidy probability data
      obs.train <- data.table(occur.taxa)
      obs.train$SiteId <- sites
      obs.train$SampId <- samples
      obs.train <- gather(obs.train, Taxon, Obs, -SiteId, -SampId)
      
      p.train <- data.table(p.maxpost, stringsAsFactors = FALSE)
      colnames(p.train) <- colnames(occur.taxa)
      p.train$SiteId <- sites
      p.train$SampId <- samples
      p.train <- gather(p.train, Taxon, Pred, -SiteId, -SampId)
      p.train <- left_join(p.train, obs.train, by = c("SiteId", "SampId", "Taxon"))
      p.train$Type <- "Training"
      p.train$Fold <- k
      p.train$Model <- "jSDM"
      p.train$Trial <- folder
      
      ### TEST prob/dev
      # Read observations and inputs for testing 
      paste(full.path, '_train', k, sep="")
      obs.test <- read.csv(paste(full.path, '_train', k, '/inputs/test', k, '.csv', sep=""), header = TRUE, stringsAsFactors = F)
      predictors <- read.csv(paste(full.path, '_train', k, '/inputs/predictors.csv', sep=""), header = TRUE, sep = ',', stringsAsFactors = F)
      
      # Drop taxa that do not occur in training data
      y.test <- obs.test[, c("SiteId", "SampId", names(n.train))]
      # Join predictors to test observations
      predictors.test <- left_join(y.test[, c("SiteId", "SampId")], predictors, by=c("SiteId", "SampId"))
      
      # Drop rows in y.test and predictors.test where the latter has any NAs
      ind <- !apply(is.na(predictors.test[,inf.fact]),1,FUN=any)
      ind <- ifelse(is.na(ind),FALSE,ind)
      predictors.test <- predictors.test[ind, ]
      y.test <- y.test[ind, ]
      
      # Keep the test sites
      test.sites <- predictors.test$SiteId
      test.samples <- predictors.test$SampId
      
      # Drop the test sites and keep the inputs
      test.inputs <- as.matrix(predictors.test[, inf.fact])
      
      # Calculate predicted probabilities based on fitted parameters and test data 
      z <- matrix(rep(alpha.taxa.maxpost,nrow(test.inputs)),nrow=nrow(test.inputs),byrow=TRUE) + 
        test.inputs%*%beta.taxa.maxpost
      p.maxpost <- 1/(1+exp(-z))
      
      y.test <- y.test[, names(n.train)]
      n.test <- apply(y.test, 2, sum, na.rm = TRUE)
      n.obs.test <- apply(y.test, 2, function(j){sum(!is.na(j))})
      p.primitive.test <- apply(y.test, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
      p.primitive.test <- ifelse(p.primitive.test==0,1e-4,p.primitive.test)
      p.primitive.test <- matrix(rep(p.primitive.test,nrow(y.test)),nrow=nrow(y.test),byrow=TRUE)
      
      dev.resid.test <- sign(y.test-0.5)*sqrt(-2*(y.test*log(p.maxpost)+(1-y.test)*log(1-p.maxpost)))
      deviance.maxpost.test   <- sum(dev.resid.test^2, na.rm = T)
      deviance.primitive.test <- -2*sum(y.test*log(p.primitive.test)+(1-y.test)*log(1-p.primitive.test), na.rm = T)
      
      ### Deviance-based statistics
      deviance.maxpost.taxa.test <- apply(dev.resid.test^2,2,sum, na.rm=T) # NAs removed
      deviance.primitive.taxa.test <- -2*apply(y.test*log(p.primitive.test)+(1-y.test)*log(1-p.primitive.test),2,sum,na.rm=T)
      
      # D2
      d.test <- 1-deviance.maxpost.taxa.test/deviance.primitive.taxa.test
      
      # Relative deviance
      std.deviance.test <- deviance.maxpost.taxa.test/n.obs.test
      
      # Tidy deviance data
      test.deviance <- data.table(null.deviance = deviance.primitive.taxa.test, 
                                  residual.deviance = deviance.maxpost.taxa.test, 
                                  std.deviance = std.deviance.test, 
                                  D2 = d.test, 
                                  Type = "Testing",
                                  Fold = k,
                                  Model = "jSDM",
                                  Trial = folder,
                                  stringsAsFactors = F)
      
      test.deviance$Taxon <- names(std.deviance.test)
      test.deviance$n.samples <- n.obs.test[test.deviance$Taxon]
      test.deviance$n <- n.test[test.deviance$Taxon]
      
      test.deviance <- test.deviance[, c("Taxon", "Type", "Fold", "Model", "Trial", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n")]
      
      # tidy the probability data
      obs.test <- setDT(obs.test)
      obs.test <- gather(obs.test, Taxon, Obs, -SiteId, -SampId)
      
      p.test <- data.table(p.maxpost, stringsAsFactors = FALSE)
      colnames(p.test) <- colnames(occur.taxa)
      p.test$SiteId <- test.sites
      p.test$SampId <- test.samples
      p.test <- gather(p.test, Taxon, Pred, -SiteId, -SampId)
      p.test <- left_join(p.test, obs.test, by = c("SiteId", "SampId", "Taxon"))
      p.test$Type <- "Testing"
      p.test$Fold <- k
      p.test$Model <- "jSDM"
      p.test$Trial <- folder
      
      ### Bind the k-fold results
      output$deviance <- bind_rows(output$deviance, train.deviance)
      output$deviance <- bind_rows(output$deviance, test.deviance)
      output$probability <- bind_rows(output$probability, p.train, p.test)
      
      cat('Trial / fold: ', folder, "/", k, '\n')
    }
  }
  return(output)
}

# Calculate predicted probabilities over the full posterior (warning: very memory intensive)
cv.jsdm.pred <- function(trials){
  pred <- list()
  # pred <- data.table()
  for (t in 1:length(trials)){
    trial <- trials[t]
    for (k in 1:3){
      model.image <- new.env()
      path <- paste('outputs/jsdm_p1/jsdm_',trial,'_train', k, sep="")
      load(paste(path, "/Inv_JSDM_D1.RData", sep = ""), envir = model.image)
      
      # Extract objects from workspace image to local environment
      sites <- model.image$sites
      samples <- model.image$samples
      occur.taxa <- model.image$occur.taxa
      env.cond <- model.image$env.cond
      inf.fact <- model.image$inf.fact
  
      jsdm <- extract(model.image$res.D1,permuted=TRUE,inc_warmup=FALSE)
      rm(model.image)
      
      # Name dimensions of beta^taxa array
      dimnames(jsdm[["beta_taxa"]])[[1]] <- 1:dim(jsdm[["beta_taxa"]][,,])[1]
      dimnames(jsdm[["beta_taxa"]])[[2]] <- inf.fact
      dimnames(jsdm[["beta_taxa"]])[[3]] <- colnames(occur.taxa)
      
      # Name dimensions of alpha^taxa array
      dimnames(jsdm[["alpha_taxa"]])[[1]] <- 1:nrow(jsdm[["alpha_taxa"]])
      dimnames(jsdm[["alpha_taxa"]])[[2]] <- colnames(occur.taxa)
      
      beta.taxa <- jsdm[["beta_taxa"]]
      alpha.taxa <- jsdm[["alpha_taxa"]]
      
      # PREDICTION: Prepare the predictive samples and inputs
      obs.test <- read.csv(paste(path, '/inputs/test', k, '.csv', sep = ''), header = TRUE, stringsAsFactors = F)
      predictors <- read.csv(paste(path, "/inputs/predictors.csv", sep = ""), header = TRUE, sep = ',', stringsAsFactors = F)
      
      # Drop taxa that do not occur in training data
      obs.test <- obs.test[, c("SiteId", "SampId", colnames(occur.taxa))]
      # Join predictors to test observations
      predictors.test <- left_join(obs.test[, c("SiteId", "SampId")], predictors, by=c("SiteId", "SampId"))
      
      # Drop rows in obs.test and predictors.test where the latter has any NAs
      ind <- !apply(is.na(predictors.test[,inf.fact]),1,FUN=any)
      ind <- ifelse(is.na(ind),FALSE,ind)
      predictors.test <- predictors.test[ind, ]
      obs.test <- obs.test[ind, ]
      rm(ind)
      
      # PREDICTION: Propogate posterior taxon-specific parameter distributions
      # through the joint model
      
      # Transpose matrix of inputs (without ID columns)
      X <- t(as.matrix(predictors.test[, inf.fact]))
      
      # Initialize array of dimensions s (samples),i (sites),j ()
      # For each taxon, compute x_ik*beta_jk
      m <- array(
        apply(beta.taxa, 3, function(beta.K){
          beta.K %*% X
        }), 
        dim=c(dim(beta.taxa)[1], nrow(predictors.test), dim(beta.taxa)[3])
        # e.g., dim=c(10000, 580, 245) for full community
      )
      rm(beta.taxa); gc() # recover memory
      
      # Calculate predicted probabilities through the link function
      # dimensions s, j, i
      pred.fold <- array(
        apply(m, 2, function(xb){
          link.function(alpha.taxa + xb)
        }),
        dim=c(nrow(alpha.taxa), ncol(alpha.taxa), nrow(predictors.test))
        # e.g., dim=c(10000, 245, 580) for full community
      )
      rm(m, alpha.taxa); gc() # recover memory
      
      dimnames(pred.fold)[[1]] <- 1:dim(pred.fold)[1]
      dimnames(pred.fold)[[2]] <- colnames(occur.taxa)
      dimnames(pred.fold)[[3]] <- predictors.test$SampId
      
      pred.fold <- melt(pred.fold)
      setDT(pred.fold)
      colnames(pred.fold) <- c("Sample", "Taxon", "SampId", "Pred")
      pred.fold <- select(pred.fold, Taxon, SampId, Pred)
      
      # Append results and print output
      pred[[paste("fold",k,sep="")]] <- pred.fold
      # pred <- bind_rows(pred, pred.fold)
      rm(pred.fold)
      cat("Fold completed:", k, "\n")
    }
  }
  gc()
  return(pred)

  # # Initialize array of dimensions s,i,k
  # X <- t(X) # transpose X prior to using it in apply()
  # # For each taxon, compute x*beta for each posterior sample, site, and variable
  # m <- array(
  #   apply(beta.taxa, 3, function(beta.K){
  #   beta.K %*% X
  #     }), 
  #   dim=c(10000, 580, 245)
  # )
  # 
  # gc()
  # 
  # # For each site, compute the probabilities
  # predictions <- array(
  #   apply(m, 2, function(x.beta){
  #     link.function(alpha.taxa + x.beta)
  #   }),
  #   dim=c(10000, 245, 580)
  # )
  
  # Parallelized implementation of probabilistic predictions in JSDM
  # library(snow)
  # cl <- makeCluster(4,type="SOCK")
  # clusterExport(cl, list="X", envir = .GlobalEnv)
  # m <- array(
  #   parApply(cl=cl, beta.taxa, MARGIN=3, function(beta.K){
  #     beta.K %*% X
  #   }), 
  #   dim=c(10000, 580, 245)
  # )
  # stopCluster(cl)
  # 
  # cl <- makeCluster(4,type="SOCK")
  # clusterExport(cl, list="alpha.taxa", envir = .GlobalEnv)
  # 
  # predictions <- array(
  #   parApply(cl=cl, m, 2, function(x.beta){
  #     link.function(alpha.taxa + x.beta)
  #   }),
  #   dim=c(10000, 245, 580)
  # )
  # stop.cluster(cl)
}

# Process models ####
# Check R workspace image filename for optional joint model extensions
# Requires path containing one workspace image as function argument
# Identifies only very specific filename extensions (use with caution!)
identify.jsdm <- function(full.path){
  # Identify the RData workspace image
  files <- list.files(full.path)
  workspace <- files[endsWith(files, ".RData")] # failure point: assumes one .RData file
  
  # Identify the individual joint model extensions
  extensions <- gsub(".RData|Inv_JSDM_", "", workspace)
  extensions <- strsplit(extensions, "_")
  
  # Check for community parameter correlations
  if ("corr" %in% extensions[[1]]){
    correlations <- 1
  } else if ("nocorr" %in% extensions[[1]]){
    correlations <- 0
  }
  
  # Check for site-specific effects
  if ("site" %in% extensions[[1]]){
    site.effect <- 1
  } else if ("nosite" %in% extensions[[1]]){
    site.effect <- 0
  }
  
  # Check for latent variables
  latvar <- strsplit(extensions[[1]][3], "latvar")
  if ("no" %in% latvar){
    n.latent <- 0
  } else  if (latvar > 0){
    n.latent <- as.integer(latvar)
  } else{
    n.latent <- 0
  }
  model.info <- list("workspace" = workspace,
                     "correlations" = correlations,
                     "site.effects" = site.effect,
                     "n.latent" = n.latent)
  return(model.info)
}

# Modified function of extract.jsdm(): Extracts results from directory containing folder or folders of JSDM results
# into one tidy output (list of data tables and the model)
# directory <- 'outputs/jsdm_p1'
# folders <- 'TT1'

# jsdm <- extract.jsdm.extensions(directory, folder)

extract.jsdm <- function(directory, folders) {
  # Store multiple model data results
  output <- list("deviance" = data.table(), "probability" = data.table(), "parameters" = data.table())
  
  # Loop through JSDMs, loading each workspace image into the local environment
  for (f in 1:length(folders)){
    cat("Model: ", folder, " -> | ")
    
    cat("Load image -> | ")
    
    folder <- folders[f]
    full.path <- paste(getwd(), directory, folder, sep="/")
    extensions <- identify.jsdm(full.path)
    model.image <- new.env()
    load(paste(full.path, extensions$workspace, sep = "/"), envir = model.image)
    
    # Stanfit object
    res <- model.image$res

    # Peter's code START:
    # correct chains for multiple local maxima of latent variables:
    res.orig <- res
    res.extracted.trace1 <- extract(res,permuted=FALSE,inc_warmup=FALSE)
   
    if ( extensions$n.latent > 0 )
    {
      name.x    <- "x_lat"
      name.beta <- "beta_lat"
      parnames <- names(res@sim$samples[[1]])
      ind.x    <- which(substring(parnames,1,nchar(name.x))==name.x)
      ind.beta <- which(substring(parnames,1,nchar(name.beta))==name.beta)
      ind.notwarmup <- round(res@sim$warmup/res@sim$thin+1):round(res@sim$iter/res@sim$thin)
      if ( n.latent == 1 )
      {
        # calculate mean of x (make it positive for chain 1 for compatibility across runs):
        s <- 0
        for ( ind in ind.x ) s <- s + sum(res@sim$samples[[1]][[ind]][ind.notwarmup])
        if ( s < 0 )  # change sign of both x_lat and beta_lat to keep results equal
        {
          for ( ind in ind.x )    res@sim$samples[[1]][[ind]] <- - res@sim$samples[[1]][[ind]]
          for ( ind in ind.beta ) res@sim$samples[[1]][[ind]] <- - res@sim$samples[[1]][[ind]]
        }
        if ( length(res@sim$samples) > 1 )  # adapt other chains to first:
        {
          for ( ch in 2:length(res@sim$samples) )
          {
            sumsq1 <- 0
            sumsq2 <- 0
            for ( ind in ind.x ) 
            {
              sumsq1 <- sumsq1 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]-res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
              sumsq2 <- sumsq2 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]+res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
            }
            if ( sumsq2 < sumsq1 )  # change sign of both x_lat and beta_lat to keep results equal
            {
              for ( ind in ind.x )    res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
              for ( ind in ind.beta ) res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
            }
          }
        }
      }
      else
      {
        if ( n.latent == 2 )
        {
          n.x.i      <- round(0.5*length(ind.x))
          ind.x.1    <- ind.x[1:n.x.i]
          ind.x.2    <- ind.x.1 + n.x.i
          ind.beta.1 <- ind.beta[rep(c(TRUE,FALSE),round(0.5*length(ind.beta)))]
          ind.beta.2 <- ind.beta.1 + 1
          # calculate mean of x (make it positive for chain 1 for compatibility across runs):
          s <- 0
          for ( ind in ind.x.1 ) s <- s + sum(res@sim$samples[[1]][[ind]][ind.notwarmup])
          if ( s < 0 )  # change sign of both x_lat and beta_lat to keep results equal
          {
            for ( ind in ind.x.1 )    res@sim$samples[[1]][[ind]] <- - res@sim$samples[[1]][[ind]]
            for ( ind in ind.beta.1 ) res@sim$samples[[1]][[ind]] <- - res@sim$samples[[1]][[ind]]
          }
          s <- 0
          for ( ind in ind.x.2 ) s <- s + sum(res@sim$samples[[1]][[ind]][ind.notwarmup])
          if ( s < 0 )  # change sign of both x_lat and beta_lat to keep results equal
          {
            for ( ind in ind.x.2 )    res@sim$samples[[1]][[ind]] <- - res@sim$samples[[1]][[ind]]
            for ( ind in ind.beta.2 ) res@sim$samples[[1]][[ind]] <- - res@sim$samples[[1]][[ind]]
          }
          if ( length(res@sim$samples) > 1 )  # adapt other chains to first:
          {
            for ( ch in 2:length(res@sim$samples) )
            {
              sumsq1 <- 0
              sumsq2 <- 0
              sumsq3 <- 0
              sumsq4 <- 0
              for ( ind in ind.x.1 ) 
              {
                sumsq1 <- sumsq1 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]      -res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
                sumsq2 <- sumsq2 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]      +res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
                sumsq3 <- sumsq3 + sum((res@sim$samples[[ch]][[ind+n.x.i]][ind.notwarmup]-res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
                sumsq4 <- sumsq4 + sum((res@sim$samples[[ch]][[ind+n.x.i]][ind.notwarmup]+res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
              }
              ind.min <- which.min(c(sumsq1,sumsq2,sumsq3,sumsq4))
              if ( ind.min < 3 )  # indices of latent variables ok, check signs
              {
                if ( ind.min == 2 )  # change sign of variable 1
                {
                  for ( ind in ind.x.1 )    res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
                  for ( ind in ind.beta.1 ) res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
                }
                sumsq1 <- 0           # check variable 2
                sumsq2 <- 0
                for ( ind in ind.x.2 )
                {
                  sumsq1 <- sumsq1 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]-res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
                  sumsq2 <- sumsq2 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]+res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
                }
                if ( sumsq2 < sumsq1 ) # change sign of variable 2 
                {
                  for ( ind in ind.x.2 )    res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
                  for ( ind in ind.beta.2 ) res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
                }
              }
              else   # change variables
              {
                if ( ind.min == 3 )  # change without changing sign
                {
                  for ( ind in ind.x.1 )
                  {
                    tmp <- res@sim$samples[[ch]][[ind]]
                    res@sim$samples[[ch]][[ind]] <- res@sim$samples[[ch]][[ind+n.x.i]]
                    res@sim$samples[[ch]][[ind+n.x.i]] <- tmp
                  }
                  for ( ind in ind.beta.1 ) 
                  {
                    tmp <- res@sim$samples[[ch]][[ind]]
                    res@sim$samples[[ch]][[ind]] <- res@sim$samples[[ch]][[ind+1]]
                    res@sim$samples[[ch]][[ind+1]] <- tmp
                  }
                }
                else   # change with changing sign
                {
                  for ( ind in ind.x.1 )
                  {
                    tmp <- res@sim$samples[[ch]][[ind]]
                    res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind+n.x.i]]
                    res@sim$samples[[ch]][[ind+n.x.i]] <- - tmp
                  }
                  for ( ind in ind.beta.1 ) 
                  {
                    tmp <- res@sim$samples[[ch]][[ind]]
                    res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind+1]]
                    res@sim$samples[[ch]][[ind+1]] <- - tmp
                  }
                }
                sumsq1 <- 0          # check variable 2
                sumsq2 <- 0
                for ( ind in ind.x.2 )
                {
                  sumsq1 <- sumsq1 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]-res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
                  sumsq2 <- sumsq2 + sum((res@sim$samples[[ch]][[ind]][ind.notwarmup]+res@sim$samples[[1]][[ind]][ind.notwarmup])^2)
                }
                if ( sumsq2 < sumsq1 ) # change sign of variable 2 
                {
                  for ( ind in ind.x.2 )    res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
                  for ( ind in ind.beta.2 ) res@sim$samples[[ch]][[ind]] <- - res@sim$samples[[ch]][[ind]]
                }
              }
            }
          }
        }
        else
        {
          print("*** chain merging not yet implemented for more than two latent variables ***\n")
        }
      }
    }
    
    # Extract the modified stanfit object
    res.extracted.trace2 <- extract(res,permuted=FALSE,inc_warmup=FALSE)
    res.extracted       <- extract(res,permuted=TRUE,inc_warmup=FALSE)
    # Peter's code END
    
    # Extract objects from workspace image and stanfit object
    sites <- model.image$sites # sites with complete predictors and non-zero taxa included in Stan model
    samples <- model.image$samples
    occur.taxa <- model.image$occur.taxa # non-zero site-species matrix included in Stan model
    env.cond <- model.image$env.cond # complete influence factors included in Stan model
    inf.fact <- model.image$inf.fact # influence factors included in Stan model

    # Name dimensions of posterior beta^taxa samples
    dimnames(res.extracted[["beta_taxa"]])[[1]] <- 1:dim(res.extracted[["beta_taxa"]][,,])[1]
    dimnames(res.extracted[["beta_taxa"]])[[2]] <- inf.fact
    dimnames(res.extracted[["beta_taxa"]])[[3]] <- colnames(occur.taxa)
    
    x <- as.matrix(env.cond[,inf.fact])
    colnames(x) <- inf.fact
    y <- as.matrix(occur.taxa)
    ind.maxpost <- which.max(res.extracted[["lp__"]])
    mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
    sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
    mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
    sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
    names(mu.beta.comm.maxpost) <- inf.fact
    names(sigma.beta.comm.maxpost) <- inf.fact
    
    dimnames(res.extracted[["beta_taxa"]][,,]) <- list(1:dim(res.extracted[["beta_taxa"]][,,])[1], inf.fact, colnames(occur.taxa))
    colnames(res.extracted[["alpha_taxa"]][,]) <- colnames(occur.taxa)
    
    alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,] # alpha at max posterior
    beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,] # "" ""
    rownames(beta.taxa.maxpost) <- inf.fact
    
    beta <- t(beta.taxa.maxpost)
    beta <- data.table(beta, stringsAsFactors = F)
    beta$Taxon <- colnames(occur.taxa)
    beta <- gather(beta, Variable, Parameter, -Taxon)
    
    # Match positions of site-samples to unique sites
    unique.sites <- unique(sites)
    siteIND <- match(sites, unique.sites) 
    
    if (extensions$site.effects){
      # Get named vector of site effects
      gamma.maxpost <- res.extracted[["gamma_site"]][ind.maxpost,]
      names(gamma.maxpost) <- unique.sites
      
      # Get named vectors of random effects applied to the samples
      gamma.samples <- gamma.maxpost[sites]
    }
    
    if (extensions$n.latent){
      x.lat <- res.extracted[["x_lat"]]
      beta.lat <- res.extracted[["beta_lat"]]
      
      if (extensions$n.latent==1){
        x.lat.maxpost <- x.lat[ind.maxpost,]
        names(x.lat.maxpost) <- unique.sites
        
        colnames(beta.lat) <- colnames(occur.taxa)
        beta.lat.maxpost <- beta.lat[ind.maxpost,]
        names(beta.lat.maxpost) <- colnames(occur.taxa)
      } else if (extensions$n.latent > 1){
        x.lat.maxpost <- x.lat[ind.maxpost,,]
        rownames(x.lat.maxpost) <- unique.sites
        colnames(x.lat.maxpost) <- 1:extensions$n.latent
        
        beta.lat.maxpost <- beta.lat[ind.maxpost,,]
        rownames(beta.lat.maxpost) <- 1:extensions$n.latent
        colnames(beta.lat.maxpost) <- colnames(occur.taxa)
      }
    }
    
    # Calculate new results based on extracted objects and optional extensions
    cat("Process -> | ")
    # Get the occurrence frequency to order labels
    n <- apply(y, 2, sum, na.rm=T)
    n <- sort(n, decreasing = T)  
    
    ### Probability
    # If site effects FALSE and latent variables ZERO (simplest joint model)
    if (!extensions$site.effects & extensions$n.latent==0){
      z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) + 
        x%*%beta.taxa.maxpost
    }
    # If site effects TRUE and latent variables FALSE
    if (extensions$site.effects & !extensions$n.latent){
      z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
        x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE)
    }
    # If site effects TRUE and latent variables == 1
    if (extensions$site.effects & extensions$n.latent){
      if (extensions$n.latent==1){
        lv <- sapply(beta.lat.maxpost, function(j){
          j * x.lat.maxpost[siteIND]
        })
        
        z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
          x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE) + lv
      } else if (extensions$n.latent > 1){
        # Multiply latent variables with latent coefficients for multiple latent variables
        x <- x.lat.maxpost[siteIND,]
        b <- t(beta.lat.maxpost)
        
        dim(x)#: 580, 2
        head(x, 3)
        dim(b)#: 245, 2
        head(b, 3)

        j <- b[1, ]
        
        # Attempt 1: generates correct output for one example taxa
        test1 <- t(apply(x, 1, function(i){
          j * i
        }))
        
        # # Attempt 2: generates matrix with correct results, wrong dimensions
        test2 <- apply(b, 1, function(j){ # vector of length: 2
          # cat(j," ", class(j),length(j),"\n")
          t(apply(x, 1, function(i){
            j * i
          })) # dim: 580 rows, 2 cols
        })
        
        # Attempt 3: generate array with correct results and dimensions
        test3 <- array(apply(b, 1, function(j){
          apply(x, 1, function(i){
            j * i
          })
        }),
        dim=c(n.latent, length(samples), ncol(occur.taxa)))

        identical(test1, test3[,,1])
        
        # Attempt 4: generate list of matrices with each matrix belonging to a taxon. Bind the matrices together into a 3-dimensional array 
        test4 <- lapply(rownames(b), function(j){
          row <- b[j,]
          t(apply(x, 1, function(i){
            row * i
          }))
        })
        
        # test <- do.call(abind, test4)
        test <- simplify2array(test4)
      }
    }
    p.maxpost <- 1/(1+exp(-z))
    
    p.primitive <- apply(y, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
    p.primitive <- ifelse(p.primitive==0,1e-4,p.primitive)
    p.primitive <- matrix(rep(p.primitive,nrow(y)),nrow=nrow(y),byrow=TRUE)
    dev.resid <- sign(y-0.5)*sqrt(-2*(y*log(p.maxpost)+(1-y)*log(1-p.maxpost)))
    deviance.maxpost   <- sum(dev.resid^2, na.rm = T)
    deviance.primitive <- -2*sum(y*log(p.primitive)+(1-y)*log(1-p.primitive), na.rm = T)
    
    ### Deviance
    deviance.maxpost.taxa <- apply(dev.resid^2,2,sum, na.rm=T) # NAs removed
    # Residual deviance
    deviance.primitive.taxa <- -2*apply(y*log(p.primitive)+(1-y)*log(1-p.primitive),2,sum,na.rm=T)
    deviance.fit <- 1-deviance.maxpost.taxa/deviance.primitive.taxa ### D2
    
    n.samples <- apply(y, 2, function(j){
      sum(!is.na(j))
    })
    
    deviance <- data.table(Taxon = names(deviance.fit), Model = "jSDM", null.dev = deviance.primitive.taxa, res.dev = deviance.maxpost.taxa, std.res.dev = deviance.maxpost.taxa/n.samples, D2 = deviance.fit, n.samples = n.samples[names(deviance.fit)], n = n[names(deviance.fit)], stringsAsFactors = F)
    
    row.names(deviance) <- 1:nrow(deviance)
    
    ### Tidy data
    # Melt the occurrence data into columns of Taxon and Obs
    obs <- data.table(occur.taxa)
    obs$SiteId <- sites
    obs$SampId <- samples
    obs <- gather(obs, Taxon, Obs, -SiteId, -SampId)
    
    # Melt the probabilities into columns of Taxon and Pred
    prob <- data.table(p.maxpost)
    colnames(prob) <- colnames(occur.taxa)
    prob$SiteId <- sites
    prob$SampId <- samples
    prob <- gather(prob, Taxon, Pred, -SiteId, -SampId)
    
    # Combine occurrence/probability data and join the site coordinates
    probability <- left_join(obs, prob, by = c("SiteId", "SampId", "Taxon"))
    probability <- na.omit(probability)
    rm(obs, prob)
    
    ### Prepare output
    probability$Model <- "jSDM"
    deviance$Model <- "jSDM"
    beta$Model <- "jSDM"
    
    probability$Trial <- folder
    deviance$Trial <- folder
    beta$Trial <- folder
    
    # Output results to a list
    output$extensions <- extensions
    output[[folder]]$jSDM <- TRUE # flag for select.jsdm() to identify list containing the res.extracted objects
    
    if (extensions$site.effects){
      output[[folder]]$gamma.maxpost <- gamma.maxpost
      output[[folder]]$gamma.samples <- gamma.samples
    }
    
    if (extensions$correlations){
      corrfact.comm <- res.extracted[["corrfact_comm"]]
      dimnames(corrfact.comm)[[1]] <- 1:dim(res.extracted[["corrfact_comm"]][,,])[1]
      dimnames(corrfact.comm)[[2]] <- inf.fact
      dimnames(corrfact.comm)[[3]] <- inf.fact
      
      corrfact.comm.maxpost <- corrfact.comm[ind.maxpost, ,]
      output[[folder]]$corrfact.comm <- corrfact.comm
      output[[folder]]$corrfact.comm.maxpost <- corrfact.comm.maxpost
    }
    
    if (extensions$n.latent==1){
      output[[folder]]$x.lat <- x.lat
      output[[folder]]$x.lat.maxpost <- x.lat.maxpost
      output[[folder]]$beta.lat <- beta.lat
      output[[folder]]$beta.lat.maxpost <- beta.lat.maxpost
      output[[folder]]$lv.samples <- lv
    }
    
    output$deviance <- bind_rows(output$deviance, deviance)
    output$probability <- bind_rows(output$probability, probability)
    output$parameters <- bind_rows(output$parameters, beta)
    
    # Store the site, samples, community, and input data
    output[[folder]]$sites <- sites
    output[[folder]]$samples <- samples
    output[[folder]]$occur.taxa <- occur.taxa
    output[[folder]]$inf.fact <- inf.fact
    output[[folder]]$env.cond <- env.cond
    
    # Store priors for community parameters
    output[[folder]]$mu.alpha.comm.pripar <- model.image$data$mu_alpha_comm_pripar
    output[[folder]]$sigma.alpha.comm.pripar <- model.image$data$sigma_alpha_comm_pripar
    output[[folder]]$mu.beta.comm.pripar <- model.image$data$mu_beta_comm_pripar
    output[[folder]]$sigma.beta.comm.pripar <- model.image$data$sigma_beta_comm_pripar
    
    # Store max. posterior community parameters
    output[[folder]]$mu.alpha.comm <- res.extracted[["mu_alpha_comm"]]
    output[[folder]]$mu.alpha.comm.maxpost <- mu.alpha.comm.maxpost
    output[[folder]]$sigma.alpha.comm <- res.extracted[["sigma_alpha_comm"]]
    output[[folder]]$sigma.alpha.comm.maxpost <- sigma.alpha.comm.maxpost
    
    colnames(res.extracted[["mu_beta_comm"]]) <- inf.fact
    output[[folder]]$mu.beta.comm <- res.extracted[["mu_beta_comm"]]
    output[[folder]]$mu.beta.comm.maxpost <- mu.beta.comm.maxpost
    
    colnames(res.extracted[["sigma_beta_comm"]]) <- inf.fact
    output[[folder]]$sigma.beta.comm <- res.extracted[["sigma_beta_comm"]]
    output[[folder]]$sigma.beta.comm.maxpost <- sigma.beta.comm.maxpost
    
    # Store posterior taxon-specific parameters
    output[[folder]]$alpha.taxa <- res.extracted[["alpha_taxa"]]
    output[[folder]]$alpha.taxa.maxpost <- alpha.taxa.maxpost
    
    output[[folder]]$beta.taxa  <- res.extracted[["beta_taxa"]]
    output[[folder]]$beta.taxa.maxpost <- beta.taxa.maxpost
    
    
    
    cat("DONE","\n")
  }
  return(output)
}

# Return tidy, complete Variable Selection Results (VSR; no summary statistics are calculated)
extract.vsr <- function(trial, sample){
  image <- new.env()
  
  load(paste('outputs/', trial, '/variable_selection.RData',sep=''), envir = image)
  output <- image$output
  n <- occur.freq(sample)
  rm(image)
  
  # Get the standardized deviance over k-folds, for all parameters, all taxa, all models
  d <- data.table()
  for (i in 1:length(output)){
    x <- output[[i]]
    
    # Each fold contains the same models run against each taxa in the fold;
    # repeat the models as many times as the number of taxa in each fold
    m <- apply(x$models.to.run, 1, function(m){
      paste(m, collapse = " ")
    })
    
    # Fill the Model column for each fold
    x$fold1$Model <- rep(m, nrow(x$fold1)/nrow(x$models.to.run))
    x$fold2$Model <- rep(m, nrow(x$fold2)/nrow(x$models.to.run))
    x$fold3$Model <- rep(m, nrow(x$fold3)/nrow(x$models.to.run))
    
    x$fold1$Fold <- 1
    x$fold2$Fold <- 2
    x$fold3$Fold <- 3
    
    # Bind the datasets by row into one
    d <- bind_rows(d, x$fold1, x$fold2, x$fold3)
  }
  n <- occur.freq(sample)
  d$n <- n[d$Taxon]
  return(d)
}

# Propogate quantiles of a subsample of the posterior through the jSDM
# to obtain predicted probabilites based on 5th and 95th quantiles (default arguments)
propogate.jsdm.pred <- function(result, get.quantiles=TRUE, quantiles=c(0.05, 0.95)){
  jsdm <- select.jsdm(result)
  
  # Extract objects from workspace image to local environment
  sites <- jsdm$sites
  samples <- jsdm$samples
  occur.taxa <- jsdm$occur.taxa
  env.cond <- jsdm$env.cond
  inf.fact <- jsdm$inf.fact
  
  # Name dimensions of beta^taxa array
  dimnames(jsdm[["beta_taxa"]])[[1]] <- 1:dim(jsdm[["beta_taxa"]][,,])[1]
  dimnames(jsdm[["beta_taxa"]])[[2]] <- inf.fact
  dimnames(jsdm[["beta_taxa"]])[[3]] <- colnames(occur.taxa)
  
  # Name dimensions of alpha^taxa array
  dimnames(jsdm[["alpha_taxa"]])[[1]] <- 1:nrow(jsdm[["alpha_taxa"]])
  dimnames(jsdm[["alpha_taxa"]])[[2]] <- colnames(occur.taxa)
  
  beta.taxa <- jsdm[["beta_taxa"]]
  alpha.taxa <- jsdm[["alpha_taxa"]]
  
  # Propogate a subsample of the full posterior through the model to obtain predicted probabilities,
  # THEN obtain the 5th and 95th quantiles of the predicted probabilities
  posterior.samples <- 1:nrow(alpha.taxa)
  subsample.size <- length(posterior.samples)*0.2
  
  set.seed(2352)# ensure a reproducible random sample
  subsample <- sample(posterior.samples, subsample.size, replace=FALSE)
  
  alpha.taxa.subsample <- alpha.taxa[subsample, ]
  beta.taxa.subsample <- beta.taxa[subsample, ,]

  # Transpose matrix of inputs (without ID columns)
  X <- t(as.matrix(env.cond[, inf.fact]))
  
  # Initialize array of dimensions s (samples),i (sites),j ()
  # For each taxon, compute x_ik*beta_jk
  m <- array(
    apply(beta.taxa.subsample, 3, function(beta.K){
      beta.K %*% X
    }), 
    dim=c(subsample.size, nrow(env.cond), dim(beta.taxa.subsample)[3])
    # e.g., dim=c(10000, 580, 245) for full community
  )
  
  # Calculate predicted probabilities through the link function
  # dimensions s, j, i
  prob <- array(
    apply(m, 2, function(xb){
      link.function(alpha.taxa.subsample + xb)
    }),
    dim=c(subsample.size, ncol(alpha.taxa.subsample), nrow(env.cond))
    # e.g., dim=c(10000, 245, 580) for full community
  )
  rm(m); gc()
  
  dimnames(prob)[[1]] <- subsample
  dimnames(prob)[[2]] <- colnames(occur.taxa)
  dimnames(prob)[[3]] <- samples
  
  prob <- melt(prob)
  
  # Format the dataset
  colnames(prob) <- c("Subsample", "Taxon", "SampId", "Pred")
  prob$Taxon <- as.character(prob$Taxon)
  prob$SampId <- as.character(prob$SampId)
  
  # Add the SiteId to SampId (not the best solution but it works)
  s <- data.table(SiteId = sites, SampId = samples)
  
  # Build and join the observations
  obs <- s %>%
    cbind(occur.taxa) %>%
    gather(Taxon, Obs, -SiteId, -SampId) %>%
    select(-SiteId)
  
  if(get.quantiles){
    prob <- prob %>%
      group_by(SampId, Taxon) %>%
      summarise(quantile.05 = quantile(Pred, 0.05), 
                quantile.95 = quantile(Pred, 0.95)) %>%
      gather(Quantile, Pred, quantile.05:quantile.95, -SampId) %>%
      mutate(Quantile=ifelse(Quantile=="quantile.05", 0.05, 0.95)) %>%
      left_join(obs, by=c("SampId", "Taxon")) %>%
      left_join(s, by="SampId")
  } else{
    prob <- prob %>%
      left_join(obs, by=c("SampId", "Taxon")) %>%
      left_join(s, by="SampId")
  }
  return(prob)
}

# # Extract results from calibrated jSDM
# extract.jsdm <- function(dir, trials, epsilon) {
#   if (missing(epsilon)){
#     epsilon <- FALSE
#   }
#   # Create a local environment for loading workspace image
#   model.image <- new.env()
#   
#   # Store multiple model data results
#   output <- list("deviance" = data.table(), "probability" = data.table(), "parameters" = data.table())
#   # Loop through JSDM trials, loading each workspace image into the local environment
#   for (t in 1:length(trials)){
#     cat("Trial: ", trials[t], " | ")
#     load(paste(dir,"/jsdm_",trials[t], "/Inv_JSDM_D1.RData", sep = ""), envir = model.image)
#     cat("Load image -> | ")
#     
#     # Extract objects from workspace image to local environment
#     sites <- model.image$sites # sites with complete predictors and non-zero taxa included in Stan model
#     samples <- model.image$samples
#     occur.taxa <- model.image$occur.taxa # non-zero site-species matrix included in Stan model
#     env.cond <- model.image$env.cond # complete influence factors included in Stan model
#     inf.fact <- model.image$inf.fact # influence factors included in Stan model
#     
#     set.ident <- model.image$set.ident
#     res <- model.image$res.D1
#     res.extracted <- extract(res,permuted=TRUE,inc_warmup=FALSE)
#     
#     # Name dimensions of posterior beta^taxa samples
#     dimnames(res.extracted[["beta_taxa"]])[[1]] <- 1:dim(res.extracted[["beta_taxa"]][,,])[1]
#     dimnames(res.extracted[["beta_taxa"]])[[2]] <- inf.fact
#     dimnames(res.extracted[["beta_taxa"]])[[3]] <- colnames(occur.taxa)
#     
#     # res.summary <- summary(res)
#     # res.table <- res.summary$summary
#     # 
#     # 
#     # # Extract parameter names
#     # tname <- colnames(occur.taxa)
#     # r <- row.names(res.table) # get the parameter names
#     # 
#     # # Replace the beta parameter names
#     # q <- r[grep("beta_taxa", r)]
#     # q <- gsub(".*\\[(.*)\\].*", "\\1", q)
#     # q <- strsplit(q, ",")
#     # q <- lapply(q, as.numeric)
#     # q <- lapply(q, function(x){
#     #   x <- paste(inf.fact[x[1]],"_", tname[x[2]], sep = "")
#     # })
#     # q <- unlist(q)
#     # r[grep("beta_taxa", r)] <- q
#     # 
#     # # Replace the intercept names
#     # q <- r[grep("alpha_taxa", r)]
#     # q <- gsub(".*\\[(.*)\\].*", "\\1", q)
#     # q <- as.numeric(q)
#     # q <- sapply(q, function(x){
#     #   x <- paste("alpha_", tname[x], sep = "")
#     # })
#     # r[grep("alpha_taxa", r)] <- q
#     # rm(tname, q)
#     # 
#     # dims <- dim(res.extracted.trace)
#     # 
#     # plot.nrow <- 10
#     # plot.ncol <-  5
#     # 
#     # pdf(paste("outputs/jsdm_cmodel/jsdm_", trials[t],"/",trials[t],"_traceplots.pdf", sep=""),width=8,height=12)
#     # par(mfrow=c(plot.nrow,plot.ncol),mar=c(2,2,2,0.5)+0.2) # c(bottom, left, top, right)
#     # for ( i in 1:ceiling(length(names(res))/(plot.nrow*plot.ncol)) )
#     # {
#     #   start <- (i-1)*plot.nrow*plot.ncol+1
#     #   end   <- min(start+plot.nrow*plot.ncol-1,length(names(res)))
#     #   for ( j in start:end )
#     #   {
#     #     plot(numeric(0),numeric(0),type="n",cex.axis=0.8,
#     #          xlim=c(0,dims[1]),ylim=range(res.extracted.trace[,,j]),xlab="",ylab="")
#     #     title(main=r[j], cex.main = 0.7)
#     #     # dimnames(res.extracted.trace)[[3]][j]
#     #     for ( k in 1:dims[2] ) lines(1:dims[1],res.extracted.trace[,k,j],col=k)
#     #   }
#     # }
#     # dev.off()
#     # cat("Traceplots -> | ")
#     
#     x <- as.matrix(env.cond[set.ident,inf.fact])
#     colnames(x) <- inf.fact
#     y <- as.matrix(occur.taxa[set.ident,])
#     ind.maxpost <- which.max(res.extracted[["lp__"]])
#     mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
#     sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
#     mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
#     sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
#     names(mu.beta.comm.maxpost) <- inf.fact
#     names(sigma.beta.comm.maxpost) <- inf.fact
#     
#     dimnames(res.extracted[["beta_taxa"]][,,]) <- list(1:dim(res.extracted[["beta_taxa"]][,,])[1], inf.fact, colnames(occur.taxa))
#     colnames(res.extracted[["alpha_taxa"]][,]) <- colnames(occur.taxa)
#     
#     alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,] # alpha at max posterior
#     beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,] # "" ""
#     rownames(beta.taxa.maxpost) <- inf.fact
#     
#     if (epsilon){
#       # Match positions of site-samples to unique sites
#       unique.sites <- unique(sites)
#       siteIND <- match(sites, unique.sites) 
#       
#       # Get named vector of site-specific random effects
#       eps.maxpost <- res.extracted[["eps"]][ind.maxpost,]
#       names(eps.maxpost) <- unique.sites
#       
#       # Get named vectors of random effects applied to the samples
#       eps.x <- eps.maxpost[sites]
#     }
#     
#     # Get the occurrence frequency to order labels
#     n <- apply(y, 2, sum, na.rm=T)
#     n <- sort(n, decreasing = T)  
#     
#     cat("Process -> | ")
#     
#     beta <- t(beta.taxa.maxpost)
#     beta <- data.table(beta, stringsAsFactors = F)
#     beta$Taxon <- colnames(occur.taxa)
#     beta <- gather(beta, Variable, Parameter, -Taxon)
#     
#     ### Probability
#     if (epsilon){
#       z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
#         x%*%beta.taxa.maxpost + matrix(rep(eps.x,ncol(y)),nrow=nrow(x),byrow=FALSE)
#     }else{
#       z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) + 
#         x%*%beta.taxa.maxpost
#     }
#     p.maxpost <- 1/(1+exp(-z))
#     # p.primitive <- apply(y,2,sum, na.rm = TRUE)/nrow(y) # NAs removed
#     p.primitive <- apply(y, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
#     p.primitive <- ifelse(p.primitive==0,1e-4,p.primitive)
#     p.primitive <- matrix(rep(p.primitive,nrow(y)),nrow=nrow(y),byrow=TRUE)
#     dev.resid <- sign(y-0.5)*sqrt(-2*(y*log(p.maxpost)+(1-y)*log(1-p.maxpost)))
#     deviance.maxpost   <- sum(dev.resid^2, na.rm = T)
#     deviance.primitive <- -2*sum(y*log(p.primitive)+(1-y)*log(1-p.primitive), na.rm = T)
#     
#     ### Deviance
#     deviance.maxpost.taxa <- apply(dev.resid^2,2,sum, na.rm=T) # NAs removed
#     # Residual deviance
#     deviance.primitive.taxa <- -2*apply(y*log(p.primitive)+(1-y)*log(1-p.primitive),2,sum,na.rm=T)
#     deviance.fit <- 1-deviance.maxpost.taxa/deviance.primitive.taxa ### D2
#     
#     n.samples <- apply(y, 2, function(j){
#       sum(!is.na(j))
#     })
#     
#     deviance <- data.table(Taxon = names(deviance.fit), Model = "jSDM", null.dev = deviance.primitive.taxa, res.dev = deviance.maxpost.taxa, std.res.dev = deviance.maxpost.taxa/n.samples, D2 = deviance.fit, n.samples = n.samples[names(deviance.fit)], n = n[names(deviance.fit)], stringsAsFactors = F)
#     
#     row.names(deviance) <- 1:nrow(deviance)
#     
#     ### Tidy data
#     # Melt the occurrence data into columns of Taxon and Obs
#     obs <- data.table(occur.taxa)
#     obs$SiteId <- sites
#     obs$SampId <- samples
#     obs <- gather(obs, Taxon, Obs, -SiteId, -SampId)
#     
#     # Melt the probabilities into columns of Taxon and Pred
#     prob <- data.table(p.maxpost)
#     colnames(prob) <- colnames(occur.taxa)
#     prob$SiteId <- sites
#     prob$SampId <- samples
#     prob <- gather(prob, Taxon, Pred, -SiteId, -SampId)
#     
#     # Combine occurrence/probability data and join the site coordinates
#     probability <- left_join(obs, prob, by = c("SiteId", "SampId", "Taxon"))
#     probability <- na.omit(probability)
#     rm(obs, prob)
#     
#     ### Prepare output
#     probability$Model <- "jSDM"
#     deviance$Model <- "jSDM"
#     beta$Model <- "jSDM"
#     
#     probability$Trial <- trials[t]
#     deviance$Trial <- trials[t]
#     beta$Trial <- trials[t]
#     
#     if (epsilon){
#       output$eps.maxpost <- eps.maxpost
#       output$eps.x <- eps.x
#     }
#     
#     output$deviance <- bind_rows(output$deviance, deviance)
#     output$probability <- bind_rows(output$probability, probability)
#     output$parameters <- bind_rows(output$parameters, beta)
#     
#     # Store the site, samples, community, and input data
#     output[[trials[t]]] <- res.extracted
#     output[[trials[t]]]$sites <- sites
#     output[[trials[t]]]$samples <- samples
#     output[[trials[t]]]$occur.taxa <- occur.taxa
#     output[[trials[t]]]$inf.fact <- inf.fact
#     output[[trials[t]]]$env.cond <- env.cond
#     
#     # Store priors for community parameters
#     output[[trials[t]]]$mu_alpha_comm_pripar <- model.image$data$mu_alpha_comm_pripar
#     output[[trials[t]]]$sigma_alpha_comm_pripar <- model.image$data$sigma_alpha_comm_pripar
#     output[[trials[t]]]$mu_beta_comm_pripar <- model.image$data$mu_beta_comm_pripar
#     output[[trials[t]]]$sigma_beta_comm_pripar <- model.image$data$sigma_beta_comm_pripar
#     
#     # Store max. posterior community parameters
#     output[[trials[t]]]$mu.alpha.comm.maxpost <- mu.alpha.comm.maxpost
#     output[[trials[t]]]$sigma.alpha.comm.maxpost <- sigma.alpha.comm.maxpost
#     output[[trials[t]]]$mu.beta.comm.maxpost <- mu.beta.comm.maxpost
#     output[[trials[t]]]$sigma.beta.comm.maxpost <- sigma.beta.comm.maxpost
#     
#     # Store max. posterior taxon parameters
#     output[[trials[t]]]$alpha.taxa.maxpost <- alpha.taxa.maxpost
#     output[[trials[t]]]$beta.taxa.maxpost <- beta.taxa.maxpost
#     
#     cat("DONE","\n")
#   }
#   return(output)
# }

# Extract significant responses from posterior taxon-specific parameters based on 5% p-values
extract.resp <- function(results){
  jsdm <- select.jsdm(results)
  n <- occur.freq(jsdm$occur.taxa)
  
  beta.samples <- jsdm$beta.taxa
  dimnames(beta.samples) <- list(1:dim(beta.samples)[1], jsdm$inf.fact, colnames(jsdm$occur.taxa)) # name the dimensions
  taxon.samples <- melt(beta.samples)
  colnames(taxon.samples) <- c("Sample", "Variable", "Taxon", "Value")
  taxon.samples <- setDT(taxon.samples)
  taxon.samples$Variable <- as.character(taxon.samples$Variable)
  taxon.samples$Taxon <- as.character(taxon.samples$Taxon)
  
  beta.taxa.response <- matrix(nrow=length(names(n)), ncol=length(jsdm$inf.fact))
  # beta.taxa.response[j, k] <- c
  for (k in 1:length(jsdm$inf.fact)){
    variable <- jsdm$inf.fact[k]
    samples <- taxon.samples[Variable == variable,]
    
    for (j in 1:length(names(n))){
      taxon <- names(n)[j]
      sample <- samples[Taxon==taxon,]
      # Fill matrix: beta.taxa.response
      significance <- quantile(sample$Value, probs = c(0.05, 0.95))
      
      # If 5th quantile greater than 0, set positive
      # If 95th quantile less than 0, set negative
      if (significance[1] > 0){ # significant positive
        c <- 1
      }
      if (significance[2] < 0){ # significant negative
        c <- -1
      }
      # If posterior is !(positive) AND !(negative), set grey
      if (!(significance[1] > 0) & !(significance[2] < 0)){
        c <- 0
      }
      beta.taxa.response[j, k] <- c
    }
  }
  colnames(beta.taxa.response) <- jsdm$inf.fact
  beta.taxa.response <- as.data.table(beta.taxa.response)
  beta.taxa.response$Taxon <- names(n) 
  return(beta.taxa.response)
}

# Given a jSDM result, extract the linear predictor (z) per variable/taxon/sample (z=beta*(x-mean(x))). Uncentered inputs are constructed within the function by calling prepare.inputs().
linear.predictor <- function(results){
  jsdm <- select.jsdm(results)
  
  inf.fact <- jsdm$inf.fact
  inf.fact <- inf.fact[!(inf.fact=="Temp2")]
  
  n <- occur.freq(jsdm$occur.taxa)

  # Failure point: untransformed inputs must be prepared based on inputs in JSDM
  
  x <- p[, c("SiteId", "SampId", inf.fact)]
  x <- x[x$SiteId %in% jsdm$sites & x$SampId %in% jsdm$samples, ]
  # x <- prepare.inputs(inf.fact, sample.bdms, center = F)
  # x <- select(x, -Temp2)
  
  x.mean <- apply(x[, inf.fact], 2, mean, na.rm = T)
  
  # Loop through taxa in each influence factor, with an exception for Temp (must include Temp^2)
  K <- lapply(inf.fact, function(k){
    if (k=="Temp"){
      J <- lapply(names(n), function(j){
        # Get taxon/variable data
        variable <- c(k, "Temp2")
        
        # Subset parameters
        bjk <- results$parameters[Taxon==j & Variable %in% variable,]
        bj.temp <- bjk[Variable=="Temp", "Parameter"]$Parameter
        bj.temp2 <- bjk[Variable=="Temp2", "Parameter"]$Parameter
        
        # Select input data
        x.ik <- x[,c("SiteId", "SampId", "Temp")]
        xk <- x.ik[["Temp"]]
        
        # Calculate z-scale
        z <- bj.temp*(xk-x.mean["Temp"]) + bj.temp2*(xk-x.mean["Temp"])^2
        
        # Bind the data
        dt <- data.table(SiteId = x.ik[["SiteId"]], SampId = x.ik[["SampId"]], z = z, x = xk, Taxon = j, Variable=variable[1])
        # slopes <- bind_rows(slopes, dt)
        return(dt)
      })
    } else{
      J <- lapply(names(n), function(j){
        # Filter data and calculate z-scale
        bjk <- results$parameters[Taxon==j & Variable==k,]
        bjk <- bjk$Parameter
        
        # Select input data
        x.ik <- x[,c("SiteId", "SampId", k)]
        xk <- x.ik[[k]]
        
        # Calculate z-scale
        z <- bjk*(xk-x.mean[k])
        
        # Bind the data
        dt <- data.table(SiteId = x.ik[["SiteId"]], SampId = x.ik[["SampId"]], z = z, x = xk, Taxon = j, Variable = k)
        # slopes <- bind_rows(slopes, dt)
        return(dt)
      })
    }
    slopes.k <- rbindlist(J)
    return(slopes.k)
  })
  
  # If a site effect is present, extract the linear predictor
  if(results$extensions$site.effects==1){
    gamma.maxpost.samples <- jsdm$gamma.samples
    
    dt <- data.table(SiteId = jsdm$sites, SampId = jsdm$samples, z = gamma.maxpost.samples, x = NA, Taxon = NA, Variable = "Site effect")
    K[[length(K)+1]] <- dt
  }
  
  # If a single latent variable is present, extract the linear predictor
  if(results$extensions$latent.variables==1){
    x.lat.maxpost <- jsdm$x.lat.maxpost
    x.lat.maxpost.samples <- x.lat.maxpost[jsdm$sites]
    
    dt <- as.tibble(jsdm$lv.samples)
    
    dt$SiteId <- jsdm$sites
    dt$SampId <- jsdm$samples
    dt$x <- x.lat.maxpost.samples
    dt$Variable <- "LV"
    dt <- gather(dt, Taxon, z, -SiteId, -SampId, -x, -Variable)
    dt <- select(dt, SiteId, SampId, z, x, Taxon, Variable)
    
    K[[length(K)+1]] <- dt
  }
  # Bind the list of data tables into one big data table
  slopes.K <- rbindlist(K)
  
  transparency <- n/max(n)
  transparency <- rescale(transparency, to = c(0.05, 0.15))
  slopes.K$alpha <- transparency[slopes.K$Taxon]
  
  return(slopes.K)
}
# # Extract residual deviance from a jSDM result; analyze dependence of residuals
# among taxa
# extract.residuals <- function(results){
# }

# Given a jSDM result (list), returns the full sample of the marginal posterior taxon-specific parameters in tidy data table
extract.beta <- function(results){
  # Get beta samples
  jsdm <- select.jsdm(results)
  beta.samples <- melt(jsdm$beta.taxa) # Transform 3d array into 2d data.table
  
  # Get trial name
  ind <- sapply(results, function(i){
    i$jSDM==TRUE
  })
  ind <- unlist(ind)
  
  beta.samples$Trial <- names(ind)
  
  # Format the dataset
  colnames(beta.samples) <- c("Sample", "Variable", "Taxon", "Value", "Trial")
  beta.samples$Variable <- as.character(beta.samples$Variable)
  beta.samples$Taxon <- as.character(beta.samples$Taxon)
  setDT(beta.samples)
  return(beta.samples)
}

# Maps ####
# Plot modelled probabilites across Switzerland given a run.isdm() output object
map.isdm <- function(results, fileName) {
  pdf(paste(fileName, ".pdf", sep=''), paper = 'special', width = 10.5, onefile = TRUE)
  taxa <- results$taxa
  for (j in 1:length(taxa)){
    # For the given taxa and set of models, get the correct model
    taxon <- names(taxa)[j]
    
    m <- select.isdm(results, taxon)
    dj <- D2(m)
    dj <- round(d, 2)
    fv <- m$fitted.values
    
    # Tidy the DF for plotting
    dt <- results$mdata[, c("SiteId", taxon)]
    dt <- na.omit(dt)
    dt <- cbind(dt, fv)
    colnames(dt) <- c("SiteId", "Obs", "Pred")
    dt <- left_join(dt, inputs$xy, by = "SiteId")
    dt$Obs <- as.factor(dt$Obs)
    
    g <- ggplot()
    g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
    g <- g + geom_point(data = dt, aes(X, Y, color = Obs, size = Pred), alpha = 0.35)
    # g <- g + scale_size_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 7))
    g <- g + scale_radius(limits = c(0,1), breaks =seq(0, 1, 0.2), range = c(2, 6))
    g <- g + labs(title = paste("Probability of occurrence vs observations of",paste(taxon)),
                  subtitle = paste("iSDM: y ~", paste(results$inf.fact, collapse=" ", sep = " "), "- page", j, "\nD2 = ", dj),
                  x = "",
                  y = "",
                  size = "Probability of\noccurrence")
    g <- g + theme_minimal(base_size = 15)
    g <- g + theme(plot.title = element_text(), 
                   plot.subtitle = element_text(),
                   panel.grid.major=element_line(colour="transparent"),
                   axis.text = element_blank())
    g <- g + scale_y_continuous(breaks=NULL)
    g <- g + scale_x_continuous(breaks=NULL)
    g <- g + guides(color = guide_legend(override.aes = list(size=6)))
    g <- g + scale_color_manual(name = "Observation", values=c("#FF0000", "#0077FF"), labels=c("Absence", "Presence"))
    cat("Plotting taxon: ", taxon, "\n")
    print(g)
    rm(d, fv, m)
  }
  dev.off()
}

# Plot modelled probabilites across Switzerland given an extact.jsdm() output
map.jsdm <- function(results, fileName){
  jsdm <- select.jsdm(results)
  
  # Get probability/observations, deviance, and occurrence frequency for all taxa
  probability <- left_join(results$probability, inputs$xy, by="SiteId")
  setDT(probability)
  
  deviance <- results$deviance[, c("Taxon", "D2")]
  d <- deviance$D2
  names(d) <- deviance$Taxon
  
  taxa <- occur.freq(jsdm$occur.taxa)
  
  pdf(paste(fileName, ".pdf", sep=''), paper = 'special', width = 10.5, onefile = TRUE)
  
  for (j in 1:length(taxa)){
    # For the given taxa, get the probabilities/observations
    taxon <- names(taxa)[j]
    dt <- probability[Taxon==taxon, ]
    dj <- round(d[taxon], 2)
    dt$Obs <- as.factor(dt$Obs)
    
    g <- ggplot()
    g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
    g <- g + geom_point(data = dt, aes(X, Y, color = Obs, size = Pred), alpha = 0.35)
    # g <- g + scale_size_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 7))
    g <- g + scale_radius(limits = c(0,1), breaks =seq(0, 1, 0.2), range = c(2, 6))
    g <- g + labs(title = paste("Probability of occurrence vs observations of",paste(taxon)),
                  subtitle = paste("jSDM: y ~", paste(jsdm$inf.fact, collapse=" ", sep = " "), "- page", j, "\nD2 = ", dj),
                  x = "",
                  y = "")
    g <- g + theme_minimal(base_size = 15)
    g <- g + theme(plot.title = element_text(), 
                   plot.subtitle = element_text(),
                   panel.grid.major = element_line(colour="transparent"),
                   axis.text = element_blank())
    g <- g + scale_y_continuous(breaks=NULL)
    g <- g + scale_x_continuous(breaks=NULL)
    g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
    g <- g + labs(size = "Probability of\noccurrence")
    g <- g + scale_color_manual(name = "Observation", values=c("0" = "#FF0000", "1" = "#0077FF"), labels=c("Absence", "Presence"))
    cat("Plotting taxon: ", taxon, "\n")
    print(g)
  }
  dev.off()
}

map.jsdm.taxon <- function(results, taxon, legend=TRUE){
  jsdm <- select.jsdm(results)
  
  # Get probability/observations, deviance, and occurrence frequency for all taxa
  probability <- left_join(results$probability, inputs$xy, by="SiteId")
  setDT(probability)

  # For the given taxa, get the probabilities/observations
  dt <- probability[Taxon==taxon, ]
  dt$Obs <- as.factor(dt$Obs)
  
  taxon.label <- sub("_", " ", taxon)
  
  g <- ggplot()
  g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
  g <- g + geom_point(data = dt, aes(X, Y, color = Obs, size = Pred), alpha = 0.35)
  g <- g + scale_radius(limits = c(0,1), breaks =seq(0, 1, 0.2), range = c(2, 6))
  
  # Adjust theme and axes/labels
  g <- g + scale_y_continuous(breaks = NULL)
  g <- g + scale_x_continuous(breaks = NULL)
  g <- g + theme_minimal(base_size = 15)
  g <- g + theme(plot.title = element_text(hjust = 0.5),
                 axis.title = element_blank(),
                 axis.text = element_blank,
                 panel.grid.major = element_line(colour="transparent"),
                 plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
  
  # Adjust legend and colors
  g <- g + scale_color_manual(name = "Observation", values=c("0" = "#FF0000", "1" = "#0077FF"), labels = c("Absence", "Presence"))
  g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
  g <- g + labs(title = taxon.label,
                size = "Probability of\noccurrence")
  if (!legend){
    g <- g + guides(colour=FALSE, size=FALSE)
  }
  return(g)
}

map.jsdm.pred <- function(results, fileName){
  cat("Propogating posterior quantiles through joint model...\n")
  dt <- extract.jsdm.pred(results, get.quantiles=TRUE) # Get predicted probabilities with default quantiles c(0.05, 0.95)
  dt <- left_join(dt, inputs$xy, by="SiteId")
  setDT(dt)
  
  # Format data for ggplot() aesthetics
  dt <- na.omit(dt)
  dt$Obs <- as.factor(dt$Obs)
  dt$Alpha <- ifelse(dt$Quantile==0.05, 0.65, 0.35)
  dt$Alpha <- as.factor(dt$Alpha)
  dt$Shape <- ifelse(dt$Quantile==0.05, 19, 21)
  dt$Stroke <- ifelse(dt$Quantile==0.05, 0, 0.75)
  
  # Get specific taxa and exp. variables
  jsdm <- select.jsdm(results)
  taxa <- occur.freq(jsdm$occur.taxa)
  inf.fact <- jsdm$inf.fact
  
  pdf(paste(fileName, ".pdf", sep=''), paper = 'special', width = 10.5, onefile = TRUE)
  
  # Loop through taxa, subset data by taxon and plot
  for (j in 1:length(taxa)){
    taxon <- names(taxa[j])
    
    plot.data <- dt[Taxon==taxon, ]
    
    # Map geometries
    g <- ggplot()
    g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
    g <- g + geom_point(data = plot.data, aes(X, Y, size = Pred, alpha = Alpha, color = Obs, stroke = Stroke, shape = Shape))
    
    # Configure themes and labels
    g <- g + labs(title = paste("Probability of occurrence vs observations of", taxon),
                  subtitle = paste("jSDM:", paste(inf.fact, collapse = " ", sep = " "), "- page", j),
                  x = "",
                  y = "",
                  size = "Probability of\noccurrence",
                  alpha = "Posterior",
                  color = "Observation")
    g <- g + theme_minimal(base_size = 15)
    g <- g + theme(plot.title = element_text(hjust = 0.5), 
                   plot.subtitle = element_text(hjust = 0.5),
                   panel.grid.major = element_line(colour="transparent"))
    
    # Configure legends and scales
    g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1),
                    alpha = guide_legend(override.aes = list(size=6, shape=c(19,21), stroke=c(0,0.75), color="black"), order=2),
                    color = guide_legend(override.aes = list(size=6, stroke=0), order=3))
    g <- g + scale_y_continuous(breaks=NULL)
    g <- g + scale_x_continuous(breaks=NULL)
    g <- g + scale_radius(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 6))
    g <- g + scale_color_manual(values=c("0" = "#FF0000", "1" = "#0077FF"), labels=c("Absence", "Presence"))
    g <- g + scale_alpha_manual(values=c("0.65"="0.65", "0.35"="0.35"), labels=c("5th quantile", "95th quantile"))
    g <- g + scale_shape_identity() # Plot the shape according to the data
    
    cat("Plotting taxon: ", taxon, "\n")
    print(g)
  }
  dev.off()
}
# # Beginning with a full model, drop one variable at a time (OVAT) and measure
# # the change in D-squared
# VarImp <- function(cm, predictors){
#   # Calculate % increase in residual deviance vs. full model
#   # Calculate % decrease in D2 when variable is dropped
#   variables <- colnames(select(predictors, -SiteId, -SampId))
#   
#   # Obtain deviance statistics for full model
#   fm <- RunModel(cm, predictors, 100, trace=F)
#   full.deviance <- fm$deviance; rm(fm)
#   # Fill/calculate columns: Removed, D2.change, res.dev.change
#   full.deviance$Removed <- "None"
#   full.deviance$D2.change <- 0
#   full.deviance$res.dev.change <- 0
#   if (class(predictors)[1]=="data.frame"){ setDT(predictors)}
#   for (k in 1:length(variables)){
#     v <- variables[k]
#     cat("Drop variable:", v,"\n")
#     # Select data
#     vss <- variables[-k]
#     x <- predictors[, c("SiteId", "SampId", vss), with=F]
#     # Run subset
#     fm <- RunModel(cm, x, 100, trace=F)
# 
#     # Extract deviance
#     dev.drop <- fm$deviance; rm(fm)
#     # Fill/calculate columns: Removed, D2.change, res.dev.change
#     dev.drop$Removed <- v
#     dev.drop$D2.change <- (full.deviance[Taxon==dev.drop$Taxon & Removed=="None", "D2"]-dev.drop$D2)/full.deviance[Taxon==dev.drop$Taxon & Removed=="None", "D2"]
#     dev.drop$res.dev.change <- (full.deviance[Taxon==dev.drop$Taxon & Removed=="None", "res.dev"]-dev.drop$res.dev)/full.deviance[Taxon==dev.drop$Taxon & Removed=="None", "res.dev"]
#     # Rbindlist() to full model result
#     full.deviance <- bind_rows(full.deviance, dev.drop)
#   }
#   return(full.deviance)
# }

# plotSDM + elevation map ###
# library(rasterVis)
# library(ggplot2)

# dem <- raster("inputs/dem/dem500m.tif")

# geoSDM <- function(results, fileName) {
#   # Get the coordinates for all distinct sites
#   xy <- select(inv, SiteId, X, Y)
#   xy <- distinct(xy, SiteId, X, Y)
#   
#   taxa <- results$taxa
#   d <- sapply(results$models, function(m){
#     (m$null.deviance-m$deviance)/m$null.deviance
#   })
#   
#   page <- 1
#   
#   pdf(paste('outputs/', fileName, ".pdf", sep=''), paper = 'special', width = 11, onefile = TRUE)
#   
#   for (t in 1:length(taxa)){
#     # For the given taxa and set of models, get the correct model
#     taxon <- names(taxa)[t]
#     index <- match(taxon, names(taxa))
#     m <- results$models[[index]]
#     d.squared <- round(d[t], 2)
#     fv <- m$fitted.values
#     
#     # Tidy the DF for plotting
#     dt <- results$mdata[, c("SiteId", taxon)]
#     dt <- na.omit(dt)
#     dt <- cbind(dt, fv)
#     colnames(dt) <- c("SiteId", "Observed", "Fitted")
#     dt <- left_join(dt, xy, by = "SiteId")
#     dt$PA <- ifelse(dt$Observed == 1, "blue", "red")
#     
#     # dt$Fitted <- dt$Predicted^3
#     # dt$Size <- rescale(dt$Size, to = c(0,1))
#     dt$Observed <- as.factor(dt$Observed)
#     dt$PA <- as.factor(dt$PA)
#     
#     g <- gplot(dem, maxpixels = 5e5)
#     g <- g + geom_tile(aes(fill = value))
#     g <- g + scale_fill_gradient(low = 'white', high = 'black', na.value = "transparent")
#     g <- g + coord_equal()
#     
#     
#     # g <- g + geom_point(data = inputs$ch, aes(POINT_X, POINT_Y))
#     # g <- g + geom_path(lineend = "round")
#     
#     g <- g + geom_point(data = dt, aes(X, Y, color = Observed, size = Fitted), alpha = 0.50)
#     g <- g + scale_size_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 7))
#     g <- g + labs(x = "",
#                   y = "")
#     g <- g + scale_colour_brewer(palette = "Set1")
#     g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#     g <- g + scale_y_continuous(breaks=NULL)
#     g <- g + scale_x_continuous(breaks=NULL)
#     g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
#     g <- g + labs(color = "Observation", size = "Predicted probability", fill = "Elevation (m)")
#     g <- g + theme_minimal(base_size = 15)
#     print(g)
#     cat("Plotting taxa: ", taxon, "\n")
#     page <- page + 1
#     print(g)
#   }
#   dev.off()
# }

# # Extract GLM statistics: significant parameters, deviance, and AIC
# ExtractGLM <- function(results, taxa){
#   # Extract deviance and AIC statistics from all models
#   d <- sapply(results$models, function(x){
#     c(x$deviance, x$null.deviance, x$aic)
#   })
#   d <- t(d)
#   colnames(d) <- c("res.dev", "null.dev", "aic")
#   
#   # Create tidy DF with taxa names
#   d <- as.data.frame(d, stringsAsFactors = FALSE)
#   d$Taxa <- names(results$taxa)
#   
#   
#   d$D2 <- ((d$null.dev - d$res.dev) / d$null.dev)
#   
#   # Select significant model statistics in (d) based on d$taxa in (n)
#   n <- results$parameters$Taxa
#   
#   d <- d[d[, "Taxa"] %in% n,]
#   
#   # Combine model fit statistics with parameters
#   modsum <- left_join(d, results$parameters, by = "Taxa")
#   colnames(modsum)[which(colnames(modsum)=="Taxa")] <- taxa
#   
#   # Add frequency of occurrence for taxa
#   modsum$n <- results$taxa[modsum[, taxa]]
#   return(modsum)
# }
# 
# # For each taxa, for each influence factor, fit taxa ~ inf.fact and extract D2
# varFit <- function(comm, pr){
#   # Only keep community sites that have complete predictors
#   pr <- na.omit(pr)
#   input_cm <- comm[comm$SiteId %in% pr$SiteId, ]
#   
#   # Convert to data table if necessary
#   if (class(input_cm)[1] == "data.frame"){ input_cm <- as.data.table(input_cm)}
#   if (class(pr)[1] == "data.frame"){ pr <- as.data.table(pr)}
#   
#   # Get taxa with sufficient frequency of occurrence
#   nsites <- uniqueN(input_cm$SiteId)
#   cutoff <- nsites*n
#   ntaxa <- ProjectStats(comm)
#   ntaxa <- filter(ntaxa, n_occurrence > cutoff)
#   taxa <- ntaxa$Taxa
#   
#   # Prepare data for model input
#   sites <- select(input_cm, SiteId)
#   input_cm <- lapply(input_cm[,-1,with=FALSE], function(x){as.factor(x)})
#   input_cm <- bind_cols(sites, input_cm)
#   mdata = left_join(input_cm, pr, by = "SiteId")
#   tc <- trainControl(method = "none")
# 
#   # Prepare output data structures
#   iff <- colnames(pr[, -1])
#   mdb <- matrix(data = NA, nrow = length(taxa), ncol = length(iff))
#   
#   for (k in 1:length(iff)){
#     for (t in 1:length(taxa)){
#       pform <- paste(taxa[t], "~", iff[k])
#       iform <- as.formula(pform)
# 
#       cat("Fitting:", pform, "\n")
#       m <- train(iform, family = binomial(link = "logit"), data = mdata, method = "glm", trControl = tc, maxit = 100, na.action = na.omit)
#       D2 <- (m$null.deviance - m$deviance) / m$null.deviance
#       # }
#       mdb[t, k] <- D2
#     }
#   }
#   
#   colnames(mdb) <- iff
#   rownames(mdb) <- taxa
#   mdb <- melt(mdb)
#   colnames(mdb) <- c("Taxa", "Variable", "D2")
#   return(mdb)
# }
# 
# # For each taxa, run full model and then minus each influence factor from the model
# varStep <- function(cm, predictors){
#   # duplicate environmental conditions
#   rownames(predictors) <- predictors[,"SiteId"]
#   predictors <- predictors[cm$SiteId,]
#   row.names(predictors) <- 1:nrow(predictors) # reset the row names
#   
#   inf.fact <- colnames(predictors[,-1])
#   
#   # exclude SITES that do not have full sets of influence factors:
#   ind <- !apply(is.na(predictors[,inf.fact]),1,FUN=any)
#   ind <- ifelse(is.na(ind),FALSE,ind)
#   
#   input_cm <- cm[ind, ]
#   predictors   <- predictors[ind, ]
#   
#   # Count and filter taxa by cutoff
#   ntaxa <- apply(input_cm[,-1], 2, sum, na.rm = TRUE)
#   ntaxa <- sort(ntaxa, decreasing = TRUE)
#   ntaxa <- ntaxa[ntaxa > 0]
#   taxa <- names(ntaxa)
#   
#   # Get the influence factors and prepare the output matrix
#   pr <- predictors[, -1]
#   iff <- colnames(pr)
#   col <- length(iff)+1
#   # Prepare input data
#   mdata <- cbind(input_cm, pr)
#   # Prepare output data matrix
#   mdb <- matrix(data = NA, nrow = length(taxa), ncol = col)
#   
#   for (t in 1:length(taxa)){
#     # Run the complete model
#     pform <- paste(taxa[t], "~", paste(iff, collapse = "+"))
#     iform <- as.formula(pform)
#     m <- glm(iform, family = binomial(link = "logit"), data = mdata, maxit = 100, na.action = na.omit)
#     m$data <- NULL
#     D2 <- (m$null.deviance - m$deviance) / m$null.deviance
#     mdb[t, col] <- D2
#     
#     # Drop each influence factor and re-run the model
#     for (k in 1:length(iff)){
#       pform <- paste(taxa[t], "~", paste(iff[-k], collapse = "+"))
#       iform <- as.formula(pform)
#       m <- glm(iform, family = binomial(link = "logit"), data = mdata, maxit = 100, na.action = na.omit)
#       m$data <- NULL
#       D2 <- (m$null.deviance - m$deviance) / m$null.deviance
#       mdb[t, k] <- D2 
#     }
#   }
#   
#   mdb <- data.frame(mdb, stringsAsFactors = F)
#   colnames(mdb) <- c(iff, "Full Model")
#   mdb$Taxa <- taxa
#   mdb <- gather(mdb, "Variable", "D2", -Taxa)
#   return(mdb)
# }
# 
# # For every taxa, do a stepwise logistic regression
# RunStepwise <- function(cm, predictors) {
#   # Only keep community sites that have complete predictors
#   input_cm <- cm[cm$SiteId %in% predictors$SiteId, ]
#   
#   # Convert the community/predictor data to data table if needed
#   if (class(input_cm)[1] == "data.frame"){ input_cm <- as.data.table(input_cm)}
#   if (class(predictors)[1] == "data.frame"){ predictors <- as.data.table(predictors)}
#   
#   # Filter taxa by cutoff
#   nsites <- length(unique(input_cm$SiteId))
#   cutoff <- nsites*0.05
#   ntaxa <- ProjectStats(cm)
#   ntaxa <- filter(ntaxa, n_occurrence > cutoff)
#   taxa <- ntaxa$Taxa
#   
#   # Prepare data for model input
#   sites <- select(input_cm, SiteId)
#   input_cm <- lapply(input_cm[,-1,with=FALSE], function(x){as.factor(x)})
#   cat("Model inputs:", nrow(input_cm), "sites with", ncol(input_cm), "taxa \n")
#   input_cm <- bind_cols(sites, input_cm)
#   mdata = left_join(input_cm, predictors, by = "SiteId")
#   
#   # Prepare the influence factors for the formulas
#   inf.factors <- paste(colnames(predictors[, -1, with=FALSE]), collapse = "+")
#   
#   # Pre-allocate matrix X with j rows and k columns
#   x <- matrix(data = 0, nrow = length(taxa), ncol = ncol(predictors[, -1, with = FALSE]))
#   colnames(x) <- colnames(predictors[, -1, with = FALSE])
#   rownames(x) <- taxa
#   
#   for (y in 1:length(taxa)){
#     printFormula <- paste(taxa[y], "~", inf.factors)
#     inputFormula <- as.formula(printFormula)
#     
#     cat("Fitting:", printFormula, "\n")
#     
#     # Fit the initial model and do backwards selection
#     m <- glm(inputFormula, family = binomial(link = "logit"), data = mdata, maxit = 100, na.action = na.pass)
#     mstep <- step(m, scale = 0, trace = FALSE, direction = "backward")
#     
#     # Insert the beta parameters into output matrix X
#     beta <- coefficients(mstep)
#     beta <- beta[names(beta) != "(Intercept)"] # Drop the intercept
#     x[taxa[y], names(beta)] <- unname(beta) # Insert k parameters for taxa j into matrix X
#   }
#   return(x)
# }
# 
# # Extract GAM statistics: significant splines, deviance, and AIC
# ExtractGAM <- function(results, taxa){
#   # Extract deviance and AIC statistics from all models
#   d <- lapply(results$models, function(x){
#     var <- c(x$deviance, x$null.deviance, x$aic)
#   })
#   
#   # Create tidy DF with taxa names
#   d <- do.call(rbind, d)
#   d <- cbind(d, results$taxa$Taxa)
#   d <- as.data.frame(d, stringsAsFactors = FALSE)
#   colnames(d) <- c("gam.res.dev", "gam.null.dev", "gam.aic", taxa)
#   d[, 1] <- as.numeric(d[, 1])
#   d[, 2] <- as.numeric(d[, 2])
#   d[, 3] <- as.numeric(d[, 3])
#   
#   d$gam.pdev <- ((d[,2]-d[,1])/d[,2])*100
#   
#   return(d)
#   # # Select significant model statistics in (d) based on d$taxa in (n)
#   # n <- row.names(results$effects)
#   # 
#   # d <- d[d[, taxa] %in% n,]
#   # 
#   # # Combine model statistics with significant effects
#   # e <- as.data.frame(results$effects)
#   # e <- add_rownames(e, var = taxa)
#   # modsum <- left_join(d, e, by = taxa)
#   # 
#   # # Add frequency of occurrence for taxa
#   # z <- as.data.frame(results$taxa)
#   # colnames(z) <- c(taxa, "n")
#   # z <- z[z[, taxa] %in% n,]
#   # modsum <- left_join(modsum, z, by = taxa)
# }
# 
# # Adds higher taxonomic levels to model summaries
# # AddTaxa(proj, proj.glm, "Family")
# AddTaxa <- function(ProjectData, ModelObject, taxa) {
#   if (taxa == "Family"){
#     dt <- select(ProjectData, Phylum, Class, Order, Family)
#     dt <- distinct(dt)
#   }
#   
#   if (taxa == "Genus") {
#     dt <- select(ProjectData, Phylum, Class, Order, Family, Genus)
#     dt <- distinct(dt)
#   }
#   
#   if (taxa == "Species") {
#     dt <- select(ProjectData, Phylum, Class, Order, Family, Genus, Species)
#     dt <- distinct(dt)
#   }
#   
#   #  dt <- dt[!duplicated(dt[, "Family"]),]
#   ModelObject$summary <- left_join(ModelObject$summary, dt, by = taxa)
# }
# 
# # Moves columns within results$summary table; function is called in Workflow()
# SortCol <- function(results){
#   # Desired order of all possible taxonomic levels, model statistics, and predictors
#   apc <- c("Phylum", "Class","Order", "Family", "Genus", "Species", "n", "res.dev", "null.dev", "aic", "pdev", "D2", 
#            "Urban", "Traffic", "HM", "LUI", "LUIw", "IAR","TU.Dm", "TU.Cr", "PPP", "PPP2","AgIDW", "p_farmland", "A10m", "A20m", "A50m", "Ag100m", "A100m", "A500m", "A1km", "F10m", "F20m", "F50m", "F100m", "F500m", "F1km", "G10m", "G20m", "G50m", "G100m", "G500m", "G1km", "P", "WW", "LU", "FRI", "ED", "WV", "BM", "RW", "WiVar", "BedMod", "EcoMorph", "Flow", "FV", "Velocity", "Temp", "Tempsq", "Temp2", "Shading", "Morph", "FF.int", "gam.res.dev", "gam.null.dev", "gam.pdev")
#   y <- colnames(results$summary)
#   z <- intersect(apc, y)
#   results$summary <- results$summary[z]
#   return(results)
# }
# 
# # Calls previously defined functions to output complete model summaries
# Workflow <- function(Project, Taxa, Model) {
#   # # Create community matrix for all taxa, all programs
#   # # for colonization potential index (CPI) influence factor
#   # p <- unique(inv$MonitoringProgram)
#   # p <- p[!p %in% "VD_old"]
#   # data <- filter(inv, MonitoringProgram %in% p & !is.na(inv[, get(Taxa)]))
#   # cc <- complete.cases(data[, Taxa, with = FALSE])
#   # wcm <- data[cc, ]
#   # wcm <- CC(wcm, Taxa, "PA")
#   # wcm <- as.data.frame(wcm)
#   # rm(p, data, cc)
#   
#   # Workflow: SelectData, CC, RunModel, ExtractGLM, AddTaxa, SortCol
#   dt <- SelectData(Project, Taxa)
#   comm <- CC(dt, Taxa, "PA")
#   if (Model == "GLM") {
#     m <- RunModel(comm, predictors, 0.05, "GLM")
#     m$summary <- ExtractGLM(m, Taxa)
#     m$summary <- AddTaxa(dt, m, Taxa)
#     m <- SortCol(m)
#   }
#   if (Model == "GAM") {
#     m <- RunModel(comm, predictors, 0.05, "GAM")
#     m$summary <- ExtractGAM(m, Taxa)
#   }
#   return(m)
# }


### Model development ###
# Selects invertebrate observations by programs and taxonomic resolution
# Select invertebrates by Project ("BDM"), Taxa ("Species")
# SelectData <- function(Project, Taxa) {
#   if (length(Project) > 1){
#     data <- filter(inv, MonitoringProgram %in% Project & !is.na(inv[, get(Taxa)]))
#     if (class(data)[1] == 'data.frame') {data <- as.data.table(data)}
#   }
#   else {
#     data <- filter(inv, MonitoringProgram == Project & !is.na(inv[, get(Taxa)]))
#     if (class(data)[1] == 'data.frame') {data <- as.data.table(data)}
#   }
#   return(data)
# }
# 
# # Create site-species matrix (only presence-absence works presently)
# CC <- function(Dataset, Taxa, Method) {
#   if (nrow(Dataset) > 0) {
#     data <- select(Dataset, SiteId, get(Taxa), get(Method))
#     # Sum PA of taxa by site
#     ct <- dcast(data, SiteId ~ get(Taxa), value.var = "PA", fun.aggregate = sum)
#     # Convert positive values to 1, else set values to 0
#     if (class(ct)[1] == 'data.frame') {ct <- as.data.table(ct)}
#     pa <- apply(ct[,-1, with = FALSE], c(1,2), function(x){if (x > 0) {x = 1} else x = 0})
#     pam <- cbind(ct[, 1, with = FALSE], pa)
#     names(pam)[1] <- "SiteId"
#     pam$SiteId <- as.character(pam$SiteId)
#     # pam[,-1, with=FALSE] <- apply(pam[,-1, with=FALSE], 2, function(x){as.numeric(x)})
#     return(pam)
#   }
# }

# Calculate: number of occurrences / probability of occurrence by site/ proportion of occurrences
# Depreciated function: replaced with apply() on post-processed community data
# ProjectStats <- function(cm) {
#   if (class(cm)[1] == "data.frame"){ cm <- as.data.table(cm)}
#   ntaxa <- apply(cm[,-1, with = FALSE], 2, sum, na.rm = TRUE)
#   ntaxa <- sort(ntaxa, decreasing = TRUE)
#   ntaxa <- as.data.frame(ntaxa)
#   ntaxa <-  rownames_to_column(ntaxa, "Taxon")
#   colnames(ntaxa) <- c("Taxa", "n_occurrence")
#   return(ntaxa)
# }


# # Extract GLM statistics: significant parameters, deviance, and AIC
# ExtractGLM <- function(results, taxa){
#   # Extract deviance and AIC statistics from all models
#   d <- sapply(results$models, function(x){
#     c(x$deviance, x$null.deviance, x$aic)
#   })
#   d <- t(d)
#   colnames(d) <- c("res.dev", "null.dev", "aic")
#   
#   # Create tidy DF with taxa names
#   d <- as.data.frame(d, stringsAsFactors = FALSE)
#   d$Taxa <- names(results$taxa)
#   
#   
#   d$D2 <- ((d$null.dev - d$res.dev) / d$null.dev)
#   
#   # Select significant model statistics in (d) based on d$taxa in (n)
#   n <- results$parameters$Taxa
#   
#   d <- d[d[, "Taxa"] %in% n,]
#   
#   # Combine model fit statistics with parameters
#   modsum <- left_join(d, results$parameters, by = "Taxa")
#   colnames(modsum)[which(colnames(modsum)=="Taxa")] <- taxa
#   
#   # Add frequency of occurrence for taxa
#   modsum$n <- results$taxa[modsum[, taxa]]
#   return(modsum)
# }


# # For each taxa, for each influence factor, fit taxa ~ inf.fact and extract D2
# varFit <- function(comm, pr){
#   # Only keep community sites that have complete predictors
#   pr <- na.omit(pr)
#   input_cm <- comm[comm$SiteId %in% pr$SiteId, ]
#   
#   # Convert to data table if necessary
#   if (class(input_cm)[1] == "data.frame"){ input_cm <- as.data.table(input_cm)}
#   if (class(pr)[1] == "data.frame"){ pr <- as.data.table(pr)}
#   
#   # Get taxa with sufficient frequency of occurrence
#   nsites <- uniqueN(input_cm$SiteId)
#   cutoff <- nsites*n
#   ntaxa <- ProjectStats(comm)
#   ntaxa <- filter(ntaxa, n_occurrence > cutoff)
#   taxa <- ntaxa$Taxa
#   
#   # Prepare data for model input
#   sites <- select(input_cm, SiteId)
#   input_cm <- lapply(input_cm[,-1,with=FALSE], function(x){as.factor(x)})
#   input_cm <- bind_cols(sites, input_cm)
#   mdata = left_join(input_cm, pr, by = "SiteId")
#   tc <- trainControl(method = "none")
# 
#   # Prepare output data structures
#   iff <- colnames(pr[, -1])
#   mdb <- matrix(data = NA, nrow = length(taxa), ncol = length(iff))
#   
#   for (k in 1:length(iff)){
#     for (t in 1:length(taxa)){
#       pform <- paste(taxa[t], "~", iff[k])
#       iform <- as.formula(pform)
# 
#       cat("Fitting:", pform, "\n")
#       m <- train(iform, family = binomial(link = "logit"), data = mdata, method = "glm", trControl = tc, maxit = 100, na.action = na.omit)
#       D2 <- (m$null.deviance - m$deviance) / m$null.deviance
#       # }
#       mdb[t, k] <- D2
#     }
#   }
#   
#   colnames(mdb) <- iff
#   rownames(mdb) <- taxa
#   mdb <- melt(mdb)
#   colnames(mdb) <- c("Taxa", "Variable", "D2")
#   return(mdb)
# }
# 
# # For each taxa, run full model and then minus each influence factor from the model
# varStep <- function(cm, predictors){
#   # duplicate environmental conditions
#   rownames(predictors) <- predictors[,"SiteId"]
#   predictors <- predictors[cm$SiteId,]
#   row.names(predictors) <- 1:nrow(predictors) # reset the row names
#   
#   inf.fact <- colnames(predictors[,-1])
#   
#   # exclude SITES that do not have full sets of influence factors:
#   ind <- !apply(is.na(predictors[,inf.fact]),1,FUN=any)
#   ind <- ifelse(is.na(ind),FALSE,ind)
#   
#   input_cm <- cm[ind, ]
#   predictors   <- predictors[ind, ]
#   
#   # Count and filter taxa by cutoff
#   ntaxa <- apply(input_cm[,-1], 2, sum, na.rm = TRUE)
#   ntaxa <- sort(ntaxa, decreasing = TRUE)
#   ntaxa <- ntaxa[ntaxa > 0]
#   taxa <- names(ntaxa)
#   
#   # Get the influence factors and prepare the output matrix
#   pr <- predictors[, -1]
#   iff <- colnames(pr)
#   col <- length(iff)+1
#   # Prepare input data
#   mdata <- cbind(input_cm, pr)
#   # Prepare output data matrix
#   mdb <- matrix(data = NA, nrow = length(taxa), ncol = col)
#   
#   for (t in 1:length(taxa)){
#     # Run the complete model
#     pform <- paste(taxa[t], "~", paste(iff, collapse = "+"))
#     iform <- as.formula(pform)
#     m <- glm(iform, family = binomial(link = "logit"), data = mdata, maxit = 100, na.action = na.omit)
#     m$data <- NULL
#     D2 <- (m$null.deviance - m$deviance) / m$null.deviance
#     mdb[t, col] <- D2
#     
#     # Drop each influence factor and re-run the model
#     for (k in 1:length(iff)){
#       pform <- paste(taxa[t], "~", paste(iff[-k], collapse = "+"))
#       iform <- as.formula(pform)
#       m <- glm(iform, family = binomial(link = "logit"), data = mdata, maxit = 100, na.action = na.omit)
#       m$data <- NULL
#       D2 <- (m$null.deviance - m$deviance) / m$null.deviance
#       mdb[t, k] <- D2 
#     }
#   }
#   
#   mdb <- data.frame(mdb, stringsAsFactors = F)
#   colnames(mdb) <- c(iff, "Full Model")
#   mdb$Taxa <- taxa
#   mdb <- gather(mdb, "Variable", "D2", -Taxa)
#   return(mdb)
# }
# 
# # For every taxa, do a stepwise logistic regression
# RunStepwise <- function(cm, predictors) {
#   # Only keep community sites that have complete predictors
#   input_cm <- cm[cm$SiteId %in% predictors$SiteId, ]
#   
#   # Convert the community/predictor data to data table if needed
#   if (class(input_cm)[1] == "data.frame"){ input_cm <- as.data.table(input_cm)}
#   if (class(predictors)[1] == "data.frame"){ predictors <- as.data.table(predictors)}
#   
#   # Filter taxa by cutoff
#   nsites <- length(unique(input_cm$SiteId))
#   cutoff <- nsites*0.05
#   ntaxa <- ProjectStats(cm)
#   ntaxa <- filter(ntaxa, n_occurrence > cutoff)
#   taxa <- ntaxa$Taxa
#   
#   # Prepare data for model input
#   sites <- select(input_cm, SiteId)
#   input_cm <- lapply(input_cm[,-1,with=FALSE], function(x){as.factor(x)})
#   cat("Model inputs:", nrow(input_cm), "sites with", ncol(input_cm), "taxa \n")
#   input_cm <- bind_cols(sites, input_cm)
#   mdata = left_join(input_cm, predictors, by = "SiteId")
#   
#   # Prepare the influence factors for the formulas
#   inf.factors <- paste(colnames(predictors[, -1, with=FALSE]), collapse = "+")
#   
#   # Pre-allocate matrix X with j rows and k columns
#   x <- matrix(data = 0, nrow = length(taxa), ncol = ncol(predictors[, -1, with = FALSE]))
#   colnames(x) <- colnames(predictors[, -1, with = FALSE])
#   rownames(x) <- taxa
#   
#   for (y in 1:length(taxa)){
#     printFormula <- paste(taxa[y], "~", inf.factors)
#     inputFormula <- as.formula(printFormula)
#     
#     cat("Fitting:", printFormula, "\n")
#     
#     # Fit the initial model and do backwards selection
#     m <- glm(inputFormula, family = binomial(link = "logit"), data = mdata, maxit = 100, na.action = na.pass)
#     mstep <- step(m, scale = 0, trace = FALSE, direction = "backward")
#     
#     # Insert the beta parameters into output matrix X
#     beta <- coefficients(mstep)
#     beta <- beta[names(beta) != "(Intercept)"] # Drop the intercept
#     x[taxa[y], names(beta)] <- unname(beta) # Insert k parameters for taxa j into matrix X
#   }
#   return(x)
# }
# 
# # Extract GAM statistics: significant splines, deviance, and AIC
# ExtractGAM <- function(results, taxa){
#   # Extract deviance and AIC statistics from all models
#   d <- lapply(results$models, function(x){
#     var <- c(x$deviance, x$null.deviance, x$aic)
#   })
#   
#   # Create tidy DF with taxa names
#   d <- do.call(rbind, d)
#   d <- cbind(d, results$taxa$Taxa)
#   d <- as.data.frame(d, stringsAsFactors = FALSE)
#   colnames(d) <- c("gam.res.dev", "gam.null.dev", "gam.aic", taxa)
#   d[, 1] <- as.numeric(d[, 1])
#   d[, 2] <- as.numeric(d[, 2])
#   d[, 3] <- as.numeric(d[, 3])
#   
#   d$gam.pdev <- ((d[,2]-d[,1])/d[,2])*100
#   
#   return(d)
#   # # Select significant model statistics in (d) based on d$taxa in (n)
#   # n <- row.names(results$effects)
#   # 
#   # d <- d[d[, taxa] %in% n,]
#   # 
#   # # Combine model statistics with significant effects
#   # e <- as.data.frame(results$effects)
#   # e <- add_rownames(e, var = taxa)
#   # modsum <- left_join(d, e, by = taxa)
#   # 
#   # # Add frequency of occurrence for taxa
#   # z <- as.data.frame(results$taxa)
#   # colnames(z) <- c(taxa, "n")
#   # z <- z[z[, taxa] %in% n,]
#   # modsum <- left_join(modsum, z, by = taxa)
# }

# Adds higher taxonomic levels to model summaries
# AddTaxa(proj, proj.glm, "Family")
# AddTaxa <- function(ProjectData, ModelObject, taxa) {
#   if (taxa == "Family"){
#     dt <- select(ProjectData, Phylum, Class, Order, Family)
#     dt <- distinct(dt)
#   }
#   
#   if (taxa == "Genus") {
#     dt <- select(ProjectData, Phylum, Class, Order, Family, Genus)
#     dt <- distinct(dt)
#   }
#   
#   if (taxa == "Species") {
#     dt <- select(ProjectData, Phylum, Class, Order, Family, Genus, Species)
#     dt <- distinct(dt)
#   }
#   
#   #  dt <- dt[!duplicated(dt[, "Family"]),]
#   ModelObject$summary <- left_join(ModelObject$summary, dt, by = taxa)
# }

# # Moves columns within results$summary table; function is called in Workflow()
# SortCol <- function(results){
#   # Desired order of all possible taxonomic levels, model statistics, and predictors
#   apc <- c("Phylum", "Class","Order", "Family", "Genus", "Species", "n", "res.dev", "null.dev", "aic", "pdev", "D2", 
#            "Urban", "Traffic", "HM", "LUI", "LUIw", "IAR","TU.Dm", "TU.Cr", "PPP", "PPP2","AgIDW", "p_farmland", "A10m", "A20m", "A50m", "Ag100m", "A100m", "A500m", "A1km", "F10m", "F20m", "F50m", "F100m", "F500m", "F1km", "G10m", "G20m", "G50m", "G100m", "G500m", "G1km", "P", "WW", "LU", "FRI", "ED", "WV", "BM", "RW", "WiVar", "BedMod", "EcoMorph", "Flow", "FV", "Velocity", "Temp", "Tempsq", "Temp2", "Shading", "Morph", "FF.int", "gam.res.dev", "gam.null.dev", "gam.pdev")
#   y <- colnames(results$summary)
#   z <- intersect(apc, y)
#   results$summary <- results$summary[z]
#   return(results)
# }
# 
# # Calls previously defined functions to output complete model summaries
# Workflow <- function(Project, Taxa, Model) {
#   # # Create community matrix for all taxa, all programs
#   # # for colonization potential index (CPI) influence factor
#   # p <- unique(inv$MonitoringProgram)
#   # p <- p[!p %in% "VD_old"]
#   # data <- filter(inv, MonitoringProgram %in% p & !is.na(inv[, get(Taxa)]))
#   # cc <- complete.cases(data[, Taxa, with = FALSE])
#   # wcm <- data[cc, ]
#   # wcm <- CC(wcm, Taxa, "PA")
#   # wcm <- as.data.frame(wcm)
#   # rm(p, data, cc)
#   
#   # Workflow: SelectData, CC, RunModel, ExtractGLM, AddTaxa, SortCol
#   dt <- SelectData(Project, Taxa)
#   comm <- CC(dt, Taxa, "PA")
#   if (Model == "GLM") {
#     m <- RunModel(comm, predictors, 0.05, "GLM")
#     m$summary <- ExtractGLM(m, Taxa)
#     m$summary <- AddTaxa(dt, m, Taxa)
#     m <- SortCol(m)
#   }
#   if (Model == "GAM") {
#     m <- RunModel(comm, predictors, 0.05, "GAM")
#     m$summary <- ExtractGAM(m, Taxa)
#   }
#   return(m)
# }

# # Plot overall probability of taxa by presence-absence
# plotGLM <- function(results, Taxa, fileName) {
#   pdf(paste("outputs/", fileName, ".pdf", sep=''), onefile = TRUE)
#   
#   mdt <- data.table()
#   taxa <- results$summary[, Taxa]
#   for (i in 1:length(taxa)){
#     # For the given taxa and set of models, get the correct model
#     pos <- match(taxa[i], results$taxa$Taxa)
#     m <- results$models[[pos]]
#     
#     fv <- m$fitted.values
#     
#     # Tidy the DF for plotting
#     dt <- results$mdata
#     dt <- na.omit(dt)
#     dt <- cbind(dt, fv)
#     dt <- select(dt, SiteId, eval(as.symbol(taxa[i])), fv)
#     colnames(dt) <- c("SiteId", "Observed", "Predicted")
#     mdt <- rbind(mdt, dt)
#   }
#   title <- "Overall probability versus observation"
#   yaxis <- "Overall probability"
#   xaxis <- "Observed presence-absence"
#   
#   # par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, mar = c(5,5,5,5))
#   # par(mar = c(5,5,5,5))
#   boxplot(Predicted ~ Observed, data = mdt, 
#           main = title, ylab = yaxis, xlab = xaxis, cex.names = 2, ylim=c(0,1))
#   dev.off()
# }
# 
# # Based on plotSDM(): prints boxplots of predicted probabilities for significant taxa, one taxa per page
# plotGLMs <- function(results, Taxa, fileName) {
#   
#   taxa <- results$summary[, Taxa]
#   page <- 1
#   
#   pdf(paste("output/", fileName, ".pdf", sep=''), onefile = TRUE)
#   
#   for (i in 1:length(taxa)){
#     # For the given taxa and set of models, get the correct model
#     pos <- match(taxa[i], results$taxa$Taxa)
#     m <- results$models[[pos]]
#     
#     fv <- m$fitted.values
#     
#     # Tidy the DF for plotting
#     dt <- results$mdata
#     dt <- na.omit(dt)
#     dt <- cbind(dt, fv)
#     dt <- select(dt, SiteId, eval(as.symbol(taxa[i])), fv)
#     colnames(dt) <- c("SiteId", "Observed", "Predicted")
#     
#     title <- paste(page, " - Modelled presence-absence of",paste(taxa[i]))
#     yaxis <- "Model probability"
#     xaxis <- "Observed presence-absence"
#     cat("Plotting taxa: ", taxa[i], "\n")
#     page <- page + 1
#     par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, mar = c(5,5,5,5))
#     boxplot(Predicted ~ Observed, data = dt, 
#             main = title, ylab = yaxis, xlab = xaxis, cex.names = 2, ylim=c(0,1))
#   }
#   dev.off()
# }
# 
# # Plot GAM splines for all taxa, one taxa per page
# plotGAMs <- function(results, fileName) {
#   t <- results$taxa
#   taxa <- t$Taxa
#   page <- 1
#   
#   pdf(paste("outputs/", fileName, ".pdf", sep=''), paper = 'special', width = 12, height = 9, onefile = TRUE)
#   
#   for (i in 1:length(taxa)){
#     x <- which(results$taxa$Taxa == taxa[i])
#     m <- results$models[[x]]
#     page <- page + 1
#     cat("Plotting taxa: ", taxa[i], "\n")
#     title <- paste("Spline", taxa[i])
#     plot.gam(m, pages = 1, main = title)
#   }
#   dev.off()
# }
# 
# # Plot GAM splines for one taxa
# plotGAM <- function(results, Organism) {
#   x <- which(results$taxa$Taxa == Organism)
#   m <- results$models[[x]]
#   plot.gam(m, pages = 1)
# }
# 

map.jsdm.pred.taxon <- function(results, taxon, legend=TRUE){
  # Get predicted probabilities with default quantiles c(0.05, 0.95)
  dt <- propogate.jsdm.pred(results, get.quantiles = TRUE) 
  dt <- left_join(dt, inputs$xy, by="SiteId")
  setDT(dt)
  
  # Format data for ggplot() aesthetics
  dt <- na.omit(dt)
  dt$Obs <- as.factor(dt$Obs)
  dt$Alpha <- ifelse(dt$Quantile==0.05, 0.65, 0.35)
  dt$Alpha <- as.factor(dt$Alpha)
  dt$Shape <- ifelse(dt$Quantile==0.05, 19, 21)
  dt$Stroke <- ifelse(dt$Quantile==0.05, 0, 0.75)
  
  taxon.label <- sub("_", " ", taxon)
  
  plot.data <- dt[Taxon==taxon, ]
  # Map geometries
  g <- ggplot()
  g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
  g <- g + geom_point(data = plot.data, aes(X, Y, size = Pred, alpha = Alpha, color = Obs, stroke = Stroke, shape = Shape))
  
  # Configure theme and labels
  g <- g + theme_minimal(base_size = 15)
  g <- g + theme(plot.title = element_text(hjust = 0.5),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 panel.grid.major = element_line(colour="transparent"),
                 plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
  
  g <- g + labs(title = taxon.label,
                x = "",
                y = "",
                size = "Probability of\noccurrence",
                alpha = "Posterior",
                color = "Observation")
  
  # Configure legends and scales
  g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1),
                  alpha = guide_legend(override.aes = list(size=6, shape=c(19,21), stroke=c(0,0.75), color="black"), order=2),
                  color = guide_legend(override.aes = list(size=6, stroke=0), order=3))
  g <- g + scale_y_continuous(breaks=NULL)
  g <- g + scale_x_continuous(breaks=NULL)
  g <- g + scale_radius(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 6))
  g <- g + scale_color_manual(values=c("0" = "#FF0000", "1" = "#0077FF"), labels=c("Absence", "Presence"))
  g <- g + scale_alpha_manual(values=c("0.65"="0.65", "0.35"="0.35"), labels=c("5th quantile", "95th quantile"))
  g <- g + scale_shape_identity()
  cat("Constructing ggplot for:", taxon, "\n")
  if (!legend){
    g <- g + guides(size=FALSE, alpha=FALSE, color=FALSE)
  }
  return(g)
}
# # Plot map of influence factors (0-100)
# # Requires: predictors, inputs$ch, inputs$xy, inv
map.inputs <- function(x, fileName) {
  # Get predictor names without "SiteId"
  variables <- unique(x$Variable)
  variables <- variables[variables != "Temp2"]
  x <- left_join(x, inputs$xy, by = "SiteId")
  setDT(x)
  pdf(paste(fileName, ".pdf", sep=''), paper = 'special', width = 11, onefile = TRUE)
  trials <- unique(x$Trial)
  for (t in 1:length(trials)){
    trial <- trials[t]
    for (k in 1:length(variables)){
      variable <- variables[k]
      plot.data <- x[Trial == trial & Variable == variable,]
      
      # Set up scales
      k.min <- round(min(plot.data[["Value"]], na.rm = T), 1)
      k.max <- round(max(plot.data[["Value"]], na.rm = T), 1)
      k.int <- (k.max - k.min)/5 # ; if.int <- round(if.int)
      
      g <- ggplot()
      g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
      g <- g + geom_point(data = plot.data, aes(X, Y, size = Value), alpha = 0.35)
      # g <- g + scale_size_continuous(name = variable, limits = c(k.min, k.max), breaks = seq(k.min, k.max, k.int), range = c(2, 7))
      g <- g + scale_radius(name = variable, limits = c(k.min, k.max), breaks = seq(k.min, k.max, k.int), range = c(2, 7))
      g <- g + scale_colour_brewer(palette = "Set1")
      
      g <- g + theme(panel.grid.major = element_line(colour="transparent"))
      # Plot titles and labels
      g <- g + scale_y_continuous(breaks=NULL)
      g <- g + scale_x_continuous(breaks=NULL)
      g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
      g <- g + theme_minimal()
      g <- g + labs(main = paste("Trial:", trial,"| Variable:", variable),
                    x = "",
                    y = "")
      cat("Plotting: ", trial,"/ ", variable, "\n")
      print(g)
    }
  }
  
  dev.off()
}

# Other plots ####
# Given output from extract.jsdm(), plot community and taxon-specific parameters (and priors)
plot.comm <- function(results, filename){
  jsdm <- select.jsdm(results)
  
  # ALPHA #
  alpha.taxa <- jsdm$alpha_taxa
  colnames(alpha.taxa) <- colnames(jsdm$occur.taxa)
  # alpha.taxa <- melt(alpha.taxa, variable.factor = FALSE)
  # colnames(alpha.taxa) <- c("Taxon", "Value")
  
  # Taxon with lowest SD should have highest density
  alpha.sd <- apply(alpha.taxa, 2, function(j){
    sd(j, na.rm=T)
  })
  alpha.taxon <- alpha.taxa[, names(which.min(alpha.sd))]
  ymax <- max(density(alpha.taxon)$y)
  # mu_alpha_comm ~ normal(mu_alpha_comm_pripar,sigma_alpha_comm_pripar)
  prior <- density(rnorm(length(jsdm$lp__), mean=jsdm$mu_alpha_comm_pripar, sd=jsdm$sigma_alpha_comm_pripar))
  
  mu <- jsdm$mu.alpha.comm.maxpost
  sd <- jsdm$sigma.alpha.comm.maxpost
  x <- seq(mu-4*sd,mu+4*sd,length.out=201)
  alpha.comm <- dnorm(x, mu, sd)
  
  # Set transparency based on occurrence frequency
  n <- apply(jsdm$occur.taxa,2,sum,na.rm=T)
  n <- sort(n, decreasing=T)
  transparency <- n/length(n)
  transparency <- rescale(transparency, to=c(0.2, 0.6))
  
  pdf(paste(filename,'.pdf',sep=''), onefile = TRUE)
  par(cex=1.25)
  plot(numeric(0), numeric(0), 
       xlim = range(x), ylim=c(0, ymax),
       main = "Prior/posterior alpha_comm, posterior alpha_taxa",
       ylab = "Density", xlab = "Value")
  
  for (j in 1:length(n)){
    lines(density(alpha.taxa[,j]), col=alpha("black", transparency[j]), lwd=0.5)
  }
  # Prior alpha_comm
  lines(prior, lwd=3, col="black")
  # Posterior alpha_comm
  lines(x, alpha.comm, type="l",  col = "grey50", lwd = 5)
  
  # BETA #
  beta.samples <- jsdm$beta_taxa
  dimnames(beta.samples) <- list(1:dim(beta.samples)[1], jsdm$inf.fact, colnames(jsdm$occur.taxa)) # name the dimensions
  taxon.samples <- melt(beta.samples)
  colnames(taxon.samples) <- c("Sample", "Variable", "Taxon", "Value")
  taxon.samples <- setDT(taxon.samples)
  taxon.samples$Variable <- as.character(taxon.samples$Variable)
  taxon.samples$Taxon <- as.character(taxon.samples$Taxon)
  
  # Get SD by variable, taxon
  beta.taxa.sd <- apply(beta.samples, c(2,3), sd)
  
  for (k in 1:length(jsdm$inf.fact)){
    variable <- jsdm$inf.fact[k]
    samples <- taxon.samples[Variable == variable,]
    
    # Taxon with lowest SD should have highest density
    taxon <- names(which.min(beta.taxa.sd[variable, ]))
    beta.taxon <- samples[Taxon == taxon,]
    ymax <- max(density(beta.taxon$Value)$y)
    rm(taxon, beta.taxon)
    
    # Community parameter distribution
    mu <- jsdm$mu.beta.comm.maxpost[variable]
    sd <- jsdm$sigma.beta.comm.maxpost[variable]
    x <- seq(mu-4*sd,mu+4*sd,length.out=201)
    beta.comm <- dnorm(x, mu, sd)
    # Distribution of max. posterior taxon parameters
    beta.taxa.maxpost.density <- density(jsdm$beta.taxa.maxpost[variable, ])
    
    # Match expressions to influence factors
    labels <- c("A10m" = expression(paste(beta["A10m"], " (1/%)")),
                "Temp" = expression(paste(beta["Temp"], " (1/", degree, "C)")),
                "Temp2" = expression(paste(beta["Temp"^2], " (1/", degree, "C"^2,")")),
                "FV" = expression(paste(beta["FV"], " (s/m)")),
                "bFRI" = expression(paste(beta["bFRI"], " (1/%)")),
                "FRI" = expression(paste(beta["FRI"], " (1/%)")),
                "F10m" = expression(paste(beta["F10m"], " (1/%)")),
                "IAR" = expression(paste(beta["IAR"], " 1/(spray treatments * fraction cropland)")),
                "Urban" = expression(paste(beta["Urban"], " (1/%)")),
                "LUD" = expression(paste(beta["LUD"], " (km"^2,"/CE)"))
    )
    
    plot(numeric(0), numeric(0),
         xlim = c(min(x), max(x)),
         ylim=c(0, ymax),
         xlab = labels[variable],
         ylab = "Density")
    abline(v=0)
    # Plot the taxon-specific parameters
    for (j in 1:length(names(n))){
      taxon <- names(n)[j]
      sample <- samples[Taxon==taxon,]
      sample <- sample$Value
      # Fill matrix: beta.taxa.response
      significance <- quantile(sample, probs = c(0.05, 0.95))
      
      # If 5th quantile greater than 0, set positive
      # If 95th quantile less than 0, set negative
      if (significance[1] > 0){ # significant positive
        c <- "blue"
      }
      if (significance[2] < 0){ # significant negative
        c <- "red"
      }
      # If posterior is !(positive) AND !(negative), set grey
      if (!(significance[1] > 0) & !(significance[2] < 0)){
        c <- "grey55"
      }
      
      # beta.taxa.response[j, k] <- c
      
      lines(density(sample), type="l", col=alpha(c, 0.30), lwd=1)
    }
    # Plot community parameter distribution
    lines(x, beta.comm, type="l", col = "grey50", lwd = 5)
    # Plot maximum posterior values over all taxa
    lines(beta.taxa.maxpost.density, type="l", col="black", lwd=2, lty='longdash')
  }
  dev.off()
}

# Plot probabilities versus influence factors for each taxon, given extract.jsdm() or run.isdm() output object
plot.prob <- function(results, filename){
  # Detect whether "results" argument is a collection of iSDMs or a jSDM,
  # create objects for preparation of inputs 
  if ("models" %in% names(results)){
    K <- results$inf.fact
    y <- data.table(SiteId = results$mdata$SiteId, SampId = results$mdata$SampId)
  }else{
    jsdm <- select.jsdm(results)
    K <- jsdm$inf.fact
    y <- data.table(SiteId = jsdm$sites, SampId = jsdm$samples)
  }
  
  inputs <- prepare.inputs(K, y, center = FALSE)
  inputs <- gather(inputs, Variable, Value, -SiteId, -SampId)
  inputs <- filter(inputs, Variable != "Temp2")
  
  # labeller() returns expressions for nice plot labels for each influence factor
  inputs$Label <- factor(inputs$Variable, levels = K[K !="Temp2"])
  levels(inputs$Label) <- labeller(levels(inputs$Label))
  
  # Join input data to probabilities
  prob <- left_join(results$probability, inputs, by=c("SiteId", "SampId"))
  prob$Obs <- as.factor(prob$Obs)
  setDT(prob); setkey(prob, Taxon)
  
  # Prepare taxa/D2 statistics
  results$deviance <- arrange(results$deviance, -n)
  taxa <- results$deviance$Taxon
  d <- results$deviance$D2; names(d) <- taxa
  
  pdf(paste(filename, '.pdf', sep=''), paper = 'special', height = 9, width = 12, onefile = TRUE)
  for (j in 1:length(taxa)){
    taxon <- taxa[j]
    d.squared <- round(d[taxon], 2)
    taxon.prob <- prob[Taxon == taxon, ]
    
    g <- ggplot(data = taxon.prob, aes(x = Value, y = Pred, color = Obs))
    g <- g + geom_point(alpha = 0.25)
    g <- g + theme_bw(base_size=15)
    g <- g + facet_wrap(~ Label, scales = "free_x", labeller=label_parsed, strip.position="bottom")
    g <- g + labs(title = paste("Probability of occurrence (based on the maximum marginal posterior distribution of the taxon-specific parameter) vs explanatory variables for",paste(taxon)),
                  x = "Explanatory variable",
                  y = "Probability of occurrence",
                  color = "Observation")
    g <- g + theme(strip.background = element_blank(), 
                   strip.placement = "outside",
                   plot.title = element_text(size=10))
    g <- g + scale_color_manual(name = "Observation", values=c("#FF0000", "#0077FF"), labels=c("Absence", "Presence"))
    g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
    cat("Plotting taxon: ", taxon, "\n")
    print(g)
  }
  dev.off()
}

# Plot probabilities versus influence factors for specific taxon, given extract.jsdm() or run.isdm() output object
plot.prob.taxon <- function(results, taxon, legend=TRUE){
  # Detect whether "results" argument is a collection of iSDMs or a jSDM,
  # create objects for preparation of inputs 
  if ("models" %in% names(results)){
    K <- results$inf.fact
    y <- data.table(SiteId = results$mdata$SiteId, SampId = results$mdata$SampId)
  }else{
    jsdm <- select.jsdm(results)
    K <- jsdm$inf.fact
    y <- data.table(SiteId = jsdm$sites, SampId = jsdm$samples)
  }
  
  # Prepare inputs and reshape into tidy data
  inputs <- prepare.inputs(K, y, center = FALSE)
  inputs <- gather(inputs, Variable, Value, -SiteId, -SampId)
  inputs <- filter(inputs, Variable != "Temp2")
  
  # labeller() returns expressions for nice plot labels for each influence factor
  inputs$Label <- factor(inputs$Variable, levels = K[K !="Temp2"])
  levels(inputs$Label) <- labeller(levels(inputs$Label))
  
  # Join inputs to probabilities
  prob <- left_join(results$probability, inputs, by=c("SiteId", "SampId"))
  prob$Obs <- as.factor(prob$Obs)
  setDT(prob)
  
  taxon.prob <- prob[Taxon == taxon, ]
  taxon.label <- sub("_", " ", taxon)
  
  g <- ggplot(data = taxon.prob, aes(x = Value, y = Pred, color = Obs))
  g <- g + geom_point(alpha = 0.25)
  g <- g + theme_bw(base_size=15)
  g <- g + facet_wrap(~ Label, scales = "free_x", labeller=label_parsed, strip.position="bottom")
  g <- g + labs(title = taxon.label,
    y = "Probability of occurrence",
    color = "Observation")
  g <- g + theme(strip.background = element_blank(), 
                 strip.placement = "outside",
                 axis.title.x = element_blank(),
                 plot.title = element_text(hjust = 0.5))
  g <- g + scale_color_manual(name = "Observation", values=c("#FF0000", "#0077FF"), labels=c("Absence", "Presence"))
  g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
  if (!legend){
    g <- g + guides(colour=FALSE, size=FALSE)
  }
  return(g)
}


# # Plot training and testing error for GLM / GAM
# plotError <- function(results, taxa){
#   m <- SelectModel(results, taxa)
#   
#   # Training
#   p <- m$fitted.values
#   x <- m$train[, "Epeorus_assimilis"]
#   df <- data.frame(p, x)
#   df$TT <- "Train"
#   colnames(df) <- c("Predicted", "Observed", "TT")
#   
#   # Testing
#   p <- predict(m, newdata = m$test, type = "response")
#   x <- m$test[, "Epeorus_assimilis"]
#   dfp <- data.frame(p, x)
#   dfp$TT <- "Test"
#   colnames(dfp) <- c("Predicted", "Observed", "TT")
#   
#   # Prepare data for ggplot()
#   df <- rbind(df, dfp)
#   levels(df$Observed)[levels(df$Observed)=="1"] <- "Present"
#   levels(df$Observed)[levels(df$Observed)=="0"] <- "Absent"
#   
#   # Plot data
#   t <- paste("Testing and training error for", "Epeorus_assimilis")
#   ggplot(data = df, aes(x = Observed, y = Predicted)) + 
#     geom_boxplot() +
#     facet_grid(.~TT) +
#     # ggtitle(t) +
#     labs(title = t, subtitle = s) +
#     theme_bw()
#   
#   # Density plot
#   # ggplot(data = df, aes(x = Predicted, color = Observed)) + 
#   #   scale_colour_manual(values = c("grey", "black")) +
#   #   geom_density()
# }
