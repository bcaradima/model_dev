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

# Alias functions for base R
s <- base::summary
h <- utils::head
n <- base::names
d <- base::dim

# Print total NAs per column
check.na <- function(df){
  apply(df, 2, function(column){
    sum.na <- sum(is.na(column))
  })
}

check.uniqueN <- function(df){
  require(data.table)
  apply(df, 2, function(column){
    uniqueN(column)
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
  # sum(-2*(obs*log(pred)+(1-obs)*log(1-pred)),na.rm=TRUE)
  sum(-2*log(ifelse(test.y==1, test.p.null, 1-test.p.null)), na.rm = TRUE)
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
  dir.top <- trial
  dir.input <- paste(trial,'/inputs',sep='')
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
run.isdm <- function(cm, predictors, trace=FALSE, messages=TRUE) {
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
  
  models <- lapply(taxa, function(j) {
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
    if (messages==TRUE){
      cat("iSDM:", print.formula, "\n")
    }
    return(m)
  })
  names(models) <- taxa
  
  # Get parameters
  parameters <- sapply(models, function(m) {
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
  parameters <- as.tibble(parameters)
  
  # Spread the estimates to column by predictor
  parameters <- gather(parameters, Taxa, Estimate, -Variable)
  colnames(parameters) <- c("Variable", "Taxon", "Parameter")
  parameters$Model <- "iSDM"
  parameters <- parameters[, c("Taxon", "Variable", "Parameter", "Model")]
  
  # Get deviance statistics
  deviance <- lapply(taxa, function(j){
    m <- models[[j]]
    
    # Number of samples
    n.samples <- length(na.omit(mdata[, j]))
    
    null.dev <- m$null.deviance
    res.dev <- m$deviance
    std.res.dev <- m$deviance/n.samples
    D2 <- (null.dev-res.dev)/null.dev
    tibble(Taxon = j, null.deviance = null.dev, residual.deviance = res.dev, std.deviance = std.res.dev, n.samples = n.samples, n.present = n[j], D2 = D2, Model = "iSDM")
  })
  deviance <- bind_rows(deviance)
  
  # Get probabilities
  probability <- lapply(taxa, function(j){
    m <- models[[j]]

    prob <- m$fitted.values
    
    dt <- as.tibble(mdata[, c("SiteId", "SampId", j)])
    colnames(dt) <- c("SiteId", "SampId", "Obs")
    dt <- na.omit(dt)
    dt$Taxon <- j
    dt$Pred <- m$fitted.values
    dt$Model <- "iSDM"
    dt <- dt[, c("SiteId", "SampId", "Taxon", "Obs", "Pred", "Model")]
    return(dt)
  })
  probability <- bind_rows(probability)
  
  # Calculate Tjur's R-squared
  tjur <- extract.tjur(probability)
  tjur$Model <- "iSDM"
  
  output <- list("mdata" = mdata, "community" = cm, "predictors" = model.inputs, "inf.fact" = inf.fact, "taxa" = n, "models" = models, "parameters" = parameters, "deviance" = deviance, "tjur" = tjur, "probability" = probability)
  return(output)
}

# Select models ####
# Returns glm object of specified taxon from run.isdm() output (models are stored in a named list)
select.isdm <- function(results, taxon) {
  # x <- which(names(results$taxa) == taxon)
  m <- results$models[[taxon]]
  return(m)
}

# Cross-validation ####
# Returns list of deviances and probabilities from GLMs under k-fold cross-validation (k=3); requires only predictors as argument, with sample.bdms and train/test data assumed to be in the workspace
cv.isdm <- function(predictors){
  cat("Starting k-fold cross-validation of GLMs:\n")
  output <- list("deviance" = list(), "probability" = list())
  
  folds <- 3
  for (fold in 1:folds) {
    # Get the observations
    train.data <- as.tibble(get(paste("train",fold,sep="")))
    test.data <- as.tibble(get(paste("test",fold,sep="")))
    
    # TRAINING phase
    train.isdm <- run.isdm(train.data, predictors, messages = FALSE)
    
    train.probability <- train.isdm$probability
    train.probability$Fold <- fold
    train.probability$Type <- "Training"
    
    train.deviance <- train.isdm$deviance
    train.deviance$Fold <- fold
    train.deviance$Type <- "Training"
    train.deviance <- train.deviance[, c("Taxon", "Type", "Fold", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]
    
    # TESTING phase
    train.taxa <- names(train.isdm$taxa) # get names of taxa in training data
    test.data <- test.data[, c("SiteId", "SampId", train.taxa)] # only keep taxa found in training
    inf.fact <- train.isdm$inf.fact

    # Prepare parameters for null model
    test.n.present <- occur.freq(test.data)
    test.n.samples <- apply(select(test.data, -SiteId, -SampId), 2, function(j){
      sum(!is.na(j), na.rm = TRUE)
    })

    freq <- test.n.present/test.n.samples
    freq <- ifelse(freq==0, 1e-04, freq)
    alpha.taxa <- -log(1/freq-1)

    # get predictions to calculate testing deviance for each taxa
    test.deviance <- list()
    test.probability <- list()
    
    for (j in 1:length(train.taxa)){
      taxon <- train.taxa[j]
      m <- train.isdm$models[[taxon]]
      null.parameters <- c(alpha.taxa[taxon], rep(0, length(inf.fact)))

      # get the observations
      obs <- test.data[, c("SiteId", "SampId", taxon)]
      colnames(obs) <- c("SiteId", "SampId", "Obs")
      obs <- na.omit(obs)

      # join the inputs to the observations
      test.model.data <- left_join(obs, predictors, by = c("SiteId", "SampId"))
      test.model.data <- na.omit(test.model.data)

      test.y <- test.model.data[["Obs"]]
      test.x <- as.matrix(test.model.data[, inf.fact])

      # null deviance
      z <- null.parameters[1] + test.x%*%null.parameters[-1]
      test.p.null <- 1/(1+exp(-z))
      test.deviance.null.taxon <- sum(-2*log(ifelse(test.y==1, test.p.null, 1-test.p.null)), na.rm = TRUE)

      # residual deviance
      # calculate predicted probabilities
      if (m$optim){
        z <- m$par["Intercept"] + test.x%*%m$par[-1]
      }else{
        z <- coef(m)["(Intercept)"] + test.x %*% coef(m)[inf.fact]
      }
      z <- as.vector(z)
      test.p.resid <- 1/(1+exp(-z))

      # calculate deviance statistics
      test.deviance.resid.taxon <- sum(-2*log(ifelse(test.y == 1, test.p.resid, 1-test.p.resid)), na.rm = TRUE)

      test.deviance.taxon <- tibble(# Group information
                              Taxon = taxon,
                              Type = "Testing",
                              Fold = fold,
                              Model = "iSDM",
                              # Performance information
                              null.deviance = test.deviance.null.taxon,
                              residual.deviance = test.deviance.resid.taxon,
                              std.deviance = test.deviance.resid.taxon/length(test.y),
                              D2 = (test.deviance.null.taxon - test.deviance.resid.taxon)/test.deviance.resid.taxon,
                              # Sample information
                              n.samples = length(test.y),
                              n.present = test.n.present[taxon])

      test.deviance.taxon <- test.deviance.taxon[, c("Taxon", "Type", "Fold", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]

      test.probability.taxon <- tibble(SiteId = test.model.data$SiteId,
                                 SampId = test.model.data$SampId,
                                 Taxon = taxon,
                                 Obs = test.model.data$Obs,
                                 Pred = test.p.resid,
                                 Model = "iSDM",
                                 Fold = fold,
                                 Type = "Testing")
      
      test.deviance[[j]] <- test.deviance.taxon
      test.probability[[j]] <- test.probability.taxon
    }
    # Combine training data and list of testing data frames into one data frame
    output$deviance[[fold]] <- bind_rows(train.deviance, test.deviance)
    output$probability[[fold]] <- bind_rows(train.probability, test.probability)
    
    cat("Fold", fold, "complete\n")
  }
  output$deviance <- bind_rows(output$deviance)
  output$probability <- bind_rows(output$probability)
  
  cat("DONE\n")
  return(output)
}

# cv.isdm <- function(predictors){
#   # Fit the iSDMs to the training data
#   cat("Training iSDMs -> | ")
#   bdm.1 <- run.isdm(train1, predictors, messages = FALSE)
#   bdm.2 <- run.isdm(train2, predictors, messages = FALSE)
#   bdm.3 <- run.isdm(train3, predictors, messages = FALSE)
#   
#   # Extract the predictions for the training data
#   bdm.1$probability$Fold <- 1
#   bdm.2$probability$Fold <- 2
#   bdm.3$probability$Fold <- 3
#   
#   train.probability <- bind_rows(bdm.1$probability, bdm.2$probability, bdm.3$probability)
#   train.probability$Type <- "Training"
#   
#   # Perform predictions on corresponding independent data; obtain predicted outcomes and deviance statistics
#   cat("Extract -> | ")
#   train.deviance <- lapply(1:3, function(fold){
#     isdm <- get(paste("bdm.",fold,sep=""))
#     
#     # TRAINING
#     # Tidy training statistics for output
#     train.deviance <- tibble(Taxon = names(isdm$taxa),
#                              # Statistics
#                              null.deviance = isdm$deviance$null.dev,
#                              residual.deviance = isdm$deviance$res.dev,
#                              std.deviance = isdm$deviance$res.dev / isdm$deviance$n.samples,
#                              D2 = isdm$deviance$D2,
#                              # Sample information
#                              n.samples = isdm$deviance$n.samples,
#                              n.present = isdm$taxa,
#                              # Validation information
#                              Type = "Training",
#                              Fold = fold,
#                              Model = "iSDM",
#                              stringsAsFactors = F)
#     
#     train.deviance <- train.deviance[, c("Taxon", "Type", "Fold", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]
#   })
#   train.deviance <- bind_rows(train.deviance)
#   
#   # TESTING
#   cat("Testing -> |")
#   cv <- lapply(1:3, function(fold){
#     # Prepare input/output data
#     test.data <- as.tibble(get(paste("test", fold, sep = ""))) # get test data corresponding to fold
#     taxa <- names(isdm$taxa) # get names of taxa in training data
#     test.data <- test.data[, c("SiteId", "SampId", taxa)] # only keep taxa found in training
#     inf.fact <- isdm$inf.fact
#     
#     # Calculate intercepts from test data
#     n.present <- occur.freq(test.data)
#     n.total <- apply(select(test.data, -SiteId, -SampId), 2, function(j){
#       sum(!is.na(j), na.rm = TRUE)
#     })
#     
#     freq <- n.present/n.total
#     freq <- ifelse(freq==0, 1e-04, freq)
#     alpha.taxa <- -log(1/freq-1)
#     
#     # get predictions to calculate testing deviance for each taxa
#     J <- lapply(taxa, function(j){
#       m <- isdm$models[[j]]
#       
#       # prepare parameters for null model
#       null.parameters <- c(alpha.taxa[j], rep(0, length(inf.fact)))
#       
#       # get the observations
#       obs <- test.data[, c("SiteId", "SampId", j)]
#       colnames(obs) <- c("SiteId", "SampId", "Obs")
#       obs <- na.omit(obs)
#       
#       # join the inputs to the observations
#       test.model.data <- left_join(obs, predictors, by = c("SiteId", "SampId"))
#       test.model.data <- na.omit(test.model.data)
#       
#       test.y <- test.model.data[["Obs"]]
#       test.x <- as.matrix(test.model.data[, inf.fact])
#       
#       # null deviance
#       z <- null.parameters[1] + test.x%*%null.parameters[-1]
#       test.p.null <- 1/(1+exp(-z))
#       test.deviance.null.taxon <- sum(-2*log(ifelse(test.y==1, test.p.null, 1-test.p.null)), na.rm = TRUE)
#       
#       # residual deviance
#       # calculate predicted probabilities
#       if (m$optim){
#         z <- m$par["Intercept"] + test.x%*%m$par[-1]
#       }else{
#         z <- coef(m)["(Intercept)"] + test.x %*% coef(m)[inf.fact]
#       }
#       z <- as.vector(z)
#       test.p.resid <- 1/(1+exp(-z))
#       
#       # calculate deviance statistics
#       test.deviance.resid.taxon <- sum(-2*log(ifelse(test.y == 1, test.p.resid, 1-test.p.resid)), na.rm = TRUE)
#       
#       test.deviance <- tibble(Taxon = j,
#                    # Statistics
#                    null.deviance = test.deviance.null.taxon,
#                    residual.deviance = test.deviance.resid.taxon,
#                    std.deviance = test.deviance.resid.taxon/length(test.y),
#                    D2 = (test.deviance.null.taxon - test.deviance.resid.taxon)/test.deviance.resid.taxon,
#                    # Sample information
#                    n.samples = length(test.y),
#                    n.present = n.present[j],
#                    # Group information
#                    Type = "Testing",
#                    Fold = fold,
#                    Model = "iSDM")
#       
#       test.deviance <- test.deviance[, c("Taxon", "Type", "Fold", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]
# 
#       test.probability <- tibble(SiteId = test.model.data$SiteId, 
#                    SampId = test.model.data$SampId, 
#                    Pred = test.p.resid, 
#                    Taxon = j, 
#                    Model = "iSDM",
#                    Fold = fold,
#                    Type = "Testing")
#       
#       test.probability <- left_join(test.probability, obs, by = c("SiteId", "SampId"))
#       test.probability <- test.probability[, colnames(train.probability)]
#       
#       output <- list("probability" = test.probability, "deviance" = test.deviance)
#       return(output)
#     })
#     
#     # Bind results within each fold
#     deviance <- lapply(J, function(i){i[[2]]})
#     deviance <- bind_rows(deviance)
#     
#     probability <- lapply(J, function(i){i[[1]]})
#     probability <- bind_rows(probability)
#     
#     output <- list("probability" = probability, "deviance" = deviance)
#     return(output)
#   })
#   
#   # Bind results across the folds
#   test.deviance <- lapply(cv, function(i){i[[2]]})
#   test.deviance <- bind_rows(test.deviance)
#   
#   test.probability <- lapply(cv, function(i){i[[1]]})
#   test.probability <- bind_rows(test.probability)
#   
#   # Bind training and testing results
#   deviance <- bind_rows(train.deviance, test.deviance)
#   probability <- bind_rows(train.probability, test.probability)
#   
#   output <- list("probability" = probability, "deviance" = deviance)
#   return(output)
#   cat("| DONE\n")
# }

# Extracts tidy calibration (training) data (probability, parameters, deviance), then performs prediction (testing)
# against k-fold datasets. Warning: this function assumes a certain folder structure to read jSDM k-fold CV results
cv.jsdm <- function(directory = "outputs/paper 1 extensions", folder) {
  require(mvtnorm)
  output = list("deviance" = data.table(), "probability" = data.table())
  
  for (fold in 1:3){
    model.image <- new.env()
    cat("Load image: ", folder, " -> | ")
    full.path <- paste(getwd(), directory, paste(folder,'_train', fold, sep=""), sep="/")
    extensions <- identify.jsdm(full.path)
    load(paste(full.path, extensions$workspace, sep="/"), envir = model.image)
    
    # Read additional data for testing (assumed to be in /inputs)
    test.obs <- read.csv(paste(full.path, '/inputs/test', fold, '.csv', sep=""), header = TRUE, stringsAsFactors = F)
    predictors <- read.csv(paste(full.path, '/inputs/predictors.csv', sep=""), header = TRUE, sep = ',', stringsAsFactors = F)
    
    cat("Extract -> | ")
    
    # Stanfit object
    res <- model.image$res
    
    # Peter's code START:
    # correct chains for multiple local maxima of latent variables:
    res.orig <- res
    res.extracted.trace1 <- extract(res,permuted=FALSE,inc_warmup=FALSE)
    
    n.latent <- extensions$n.latent
    n.chain <- model.image$n.chain
    if ( extensions$n.latent > 1 ){
      name.x        <- "x_lat"
      name.beta     <- "beta_lat"
      parnames      <- names(res@sim$samples[[1]])
      ind.x         <- which(substring(parnames,1,nchar(name.x))==name.x)
      ind.beta      <- which(substring(parnames,1,nchar(name.beta))==name.beta)
      ind.notwarmup <- round(res@sim$warmup/res@sim$thin+1):round(res@sim$iter/res@sim$thin)
      print("")
      print(paste("n.latent",n.latent))
      if ( n.latent == 1 )
      {
        n.x    <- length(ind.x)
        
        # calculate means of chain 1 at all sites
        x.mean.1 <- rep(NA,length=n.x)
        for ( i in 1:n.x ) x.mean.1[i] <- mean(res@sim$samples[[1]][[ind.x[i]]][ind.notwarmup])
        if ( n.chain > 1 )  # adapt other chains to first:
        {
          for ( ch in 2:n.chain )  # match chains to first chain
          {
            print(paste("chain",ch))
            # calculate means of chain ch at all sites
            x.mean.ch <- rep(NA,length=n.x)
            for ( i in 1:n.x ) x.mean.ch[i] <- mean(res@sim$samples[[ch]][[ind.x[i]]][ind.notwarmup])
            dx.min.plus  <- sum(abs(x.mean.1-x.mean.ch))
            dx.min.minus <- sum(abs(x.mean.1+x.mean.ch))
            if ( dx.min.minus < dx.min.plus )
            {
              # change sign:
              x.mean.ch <- -x.mean.ch
              res@sim$samples[[ch]][[ind.x]]    <- -res@sim$samples[[ch]][[ind.x]]
              res@sim$samples[[ch]][[ind.beta]] <- -res@sim$samples[[ch]][[ind.beta]]
              print("sign changed")
            }
          }
        }
      }
      else  # n.latent > 1
      {
        n.x.i    <- round(length(ind.x)/n.latent)
        ind.x    <- matrix(ind.x,ncol=n.latent,byrow=FALSE)    # reformat ind.x
        ind.beta <- matrix(ind.beta,ncol=n.latent,byrow=TRUE)  # reformat ind.beta
        
        # calculate means of chain 1 at all sites and for all latent variables
        x.mean.1 <- matrix(NA,nrow=n.x.i,ncol=n.latent)
        for ( i in 1:n.x.i )
        {
          for ( k in 1:n.latent ) x.mean.1[i,k] <- mean(res@sim$samples[[1]][[ind.x[i,k]]][ind.notwarmup])
        }
        
        if ( n.chain > 1 )
        {
          for ( ch in 2:n.chain )  # match chains to first chain
          {
            print(paste("chain",ch))
            # calculate means of chain ch at all sites and for all latent variables
            x.mean.ch <- matrix(NA,nrow=n.x.i,ncol=n.latent)
            for ( i in 1:n.x.i )
            {
              for ( k in 1:n.latent ) x.mean.ch[i,k] <- mean(res@sim$samples[[ch]][[ind.x[i,k]]][ind.notwarmup])
            }
            for ( k in 1:n.latent )  # match latent variables and signs
            {
              k.min <- k
              plus.min <- TRUE 
              dx.min.plus  <- sum(abs(x.mean.1[,k]-x.mean.ch[,k]))
              dx.min.minus <- sum(abs(x.mean.1[,k]+x.mean.ch[,k]))
              dx.min <- dx.min.plus
              if ( dx.min.minus < dx.min.plus ) { plus.min <- FALSE; dx.min <- dx.min.minus }
              if ( n.latent > k )
              {
                for ( kk in (k+1):n.latent )
                {
                  dx.min.plus  <- sum(abs(x.mean.1[,kk]-x.mean.ch[,kk]))
                  dx.min.minus <- sum(abs(x.mean.1[,kk]+x.mean.ch[,kk]))
                  if ( dx.min.plus  < dx.min ) { dx.min <- dx.min.plus;  plus.min <- TRUE; k.min <- kk }
                  if ( dx.min.minus < dx.min ) { dx.min <- dx.min.minus; plus.min <- FALSE; k.min <- kk }
                }
              }
              print(paste("k",k))
              print(paste("k.min",k.min))
              # exchange variables:
              tmp <- x.mean.ch[,k.min]
              x.mean.ch[,k.min] <- x.mean.ch[,k]
              if ( plus.min ) x.mean.ch[,k] <-  tmp
              else            x.mean.ch[,k] <- -tmp
              for ( i in 1:nrow(ind.x) )
              {
                tmp <- res@sim$samples[[ch]][[ind.x[i,k.min]]]
                res@sim$samples[[ch]][[ind.x[i,k.min]]] <- res@sim$samples[[ch]][[ind.x[i,k]]]
                if ( plus.min ) res@sim$samples[[ch]][[ind.x[i,k]]] <-  tmp
                else            res@sim$samples[[ch]][[ind.x[i,k]]] <- -tmp
              }
              print("exchange x completed")
              for ( j in 1:nrow(ind.beta) )
              {
                tmp <- res@sim$samples[[ch]][[ind.beta[j,k.min]]]
                res@sim$samples[[ch]][[ind.beta[j,k.min]]] <- res@sim$samples[[ch]][[ind.beta[j,k]]]
                if ( plus.min ) res@sim$samples[[ch]][[ind.beta[j,k]]] <-  tmp
                else            res@sim$samples[[ch]][[ind.beta[j,k]]] <- -tmp
              }
              print("exchange beta completed")
            }
          }
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
    
    # Name dimensions of parameters within stanfit object
    dimnames(res.extracted[["beta_taxa"]]) <- list(1:dim(res.extracted[["beta_taxa"]][,,])[1], inf.fact, colnames(occur.taxa))
    colnames(res.extracted[["alpha_taxa"]]) <- colnames(occur.taxa)
    
    colnames(res.extracted[["mu_beta_comm"]]) <- inf.fact
    colnames(res.extracted[["sigma_beta_comm"]]) <- inf.fact
    
    # Extract inputs (x), observations (y), and parameters at maximum posterior
    x <- as.matrix(env.cond[,inf.fact])
    colnames(x) <- inf.fact
    y <- as.matrix(occur.taxa)
    ind.maxpost <- which.max(res.extracted[["lp__"]])
    mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
    sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
    mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
    sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
    alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,]
    beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,]
    
    # Name dimensions of maximum posterior community parameters
    names(mu.beta.comm.maxpost) <- inf.fact
    names(sigma.beta.comm.maxpost) <- inf.fact
    rownames(beta.taxa.maxpost) <- inf.fact
    
    # # Get the occurrence frequency to order labels
    # n <- occur.freq(y)
    
    # Match positions of site-samples to unique sites
    unique.sites <- unique(sites)
    siteIND <- match(sites, unique.sites) 
    
    # Additional extraction of JSDM results based on identify.jsdm() output
    # Check extensions for site effects
    if (extensions$site.effects){
      # Get full sample for site effects
      gamma <- res.extracted[["gamma_site"]]
      
      # Get named vector of site effects
      gamma.maxpost <- gamma[ind.maxpost,]
      names(gamma.maxpost) <- unique.sites
      
      # Get named vectors of random effects applied to the samples
      gamma.samples <- gamma.maxpost[sites]
    }
    
    # Check extensions for latent variables
    if (extensions$n.latent){
      x.lat <- res.extracted[["x_lat"]]
      beta.lat <- res.extracted[["beta_lat"]]
      
      # Check for single latent variable
      if (extensions$n.latent==1){
        x.lat.maxpost <- x.lat[ind.maxpost,]
        beta.lat.maxpost <- beta.lat[ind.maxpost,]
        
        rownames(beta.lat) <- 1:nrow(beta.lat)
        colnames(beta.lat) <- colnames(occur.taxa)
        
        names(beta.lat.maxpost) <- colnames(occur.taxa)
        
        # Check for single site- or sample-specific latent variable
        if (extensions$lat.site==1){
          rownames(x.lat) <- 1:nrow(x.lat)
          colnames(x.lat) <- unique.sites
          
          names(x.lat.maxpost) <- unique.sites
        } else if (extensions$lat.site==0){
          rownames(x.lat) <- 1:nrow(x.lat)
          colnames(x.lat) <- samples
          
          names(x.lat.maxpost) <- samples
        }
        # Check for multiple latent variables  
      } else if (extensions$n.latent > 1){
        x.lat.maxpost <- x.lat[ind.maxpost,,]
        beta.lat.maxpost <- beta.lat[ind.maxpost,,]
        
        dimnames(beta.lat) <- list(1:dim(res.extracted[["beta_lat"]][,,])[1], 1:extensions$n.latent, colnames(occur.taxa))
        
        rownames(beta.lat.maxpost) <- 1:extensions$n.latent
        colnames(beta.lat.maxpost) <- colnames(occur.taxa)
        
        # Check for site- or sample-specific latent variables
        if (extensions$lat.site==1){ # Sites
          dimnames(x.lat) <- list(1:dim(res.extracted[["x_lat"]][,,])[1], unique.sites, 1:extensions$n.latent)
          
          rownames(x.lat.maxpost) <- unique.sites
          colnames(x.lat.maxpost) <- 1:extensions$n.latent
        } else if (extensions$lat.site==0){ # Samples
          dimnames(x.lat) <- list(1:dim(res.extracted[["x_lat"]][,,])[1], samples, 1:extensions$n.latent)
          
          rownames(x.lat.maxpost) <- samples
          colnames(x.lat.maxpost) <- 1:extensions$n.latent
        }
      }
    }
    
    # Calculate model outcomes and statistics based on extracted objects and optional extensions
    cat("Process -> | ")
    
    ### Calibration results
    # Check if site effects AND latent variables are disabled (FF0)
    if (!extensions$site.effects & extensions$n.latent==0){
      z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) + 
        x%*%beta.taxa.maxpost
    }
    
    # Check if site effects are enabled AND latent variables are disabled (TT0)
    if (extensions$site.effects & !extensions$n.latent){
      z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
        x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE)
    }
    
    # Check if site effects are enabled AND latent variables are enabled (TT1-TTx)
    if (extensions$site.effects & extensions$n.latent){
      # Calculate linear predictor for single latent variable
      if (extensions$n.latent==1){
        # Check if latent variable is site- or sample-specific
        if (extensions$lat.site==1){
          # Extend site-specific LV to the samples to calculate probabilities
          lv <- sapply(beta.lat.maxpost, function(j){
            j * x.lat.maxpost[siteIND]
          })
        } 
        else if (extensions$lat.site==0){
          lv <- sapply(beta.lat.maxpost, function(j){
            j * x.lat.maxpost
          })
        }
        # Calculate linear predictor for all terms (fixed+site+latent terms)
        z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
          x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE) + lv
      } 
      # Calculate linear predictor for multiple latent variables
      else if (extensions$n.latent > 1){
        # Check if latent variables are site- or sample-specific
        if (extensions$lat.site==1){
          x.lat.maxpost.samples <- x.lat.maxpost[siteIND,]
          beta.lat.t <- t(beta.lat.maxpost)
          
          # Prepare x.lat*beta.lat as a 3-dimensional array
          lv <- lapply(rownames(beta.lat.t), function(j){
            row <- beta.lat.t[j,]
            t(apply(x.lat.maxpost.samples, 1, function(i){
              row * i
            }))
          })
          # Bind list of matrices into array
          lv <- simplify2array(lv)
          rm(x.lat.maxpost.samples, beta.lat.t)
        } 
        else if (extensions$lat.site==0){
          beta.lat.t <- t(beta.lat.maxpost)
          
          # Prepare x.lat*beta.lat as a 3-dimensional array
          lv <- lapply(rownames(beta.lat.t), function(j){
            row <- beta.lat.t[j,]
            t(apply(x.lat.maxpost, 1, function(i){
              row * i
            }))
          })
          # Bind list of matrices into array
          lv <- simplify2array(lv)
          rm(beta.lat.t)
        }
        
        # Latent terms are summed by dimensions 1 and 3 (i.e., by "sites" or samples and taxa)
        z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
          x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE) + apply(lv, c(1,3), sum)
      }
    }
    
    p.maxpost <- 1/(1+exp(-z))
    train.y <- occur.taxa
    train.n.present <- apply(train.y, 2, sum, na.rm = TRUE)
    train.n.samples <- apply(train.y, 2, function(j){sum(!is.na(j))})
    
    train.p.null <- apply(train.y, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
    train.p.null <- ifelse(train.p.null==0,1e-4,train.p.null)
    train.p.null <- matrix(rep(train.p.null,nrow(train.y)),nrow=nrow(train.y),byrow=TRUE)
    
    # Calculate deviance over the samples
    train.deviance.resid <- sign(train.y-0.5)*sqrt(-2*(train.y*log(p.maxpost)+(1-train.y)*log(1-p.maxpost)))
    
    # Calculate the sum of squared residual deviance for one taxa
    # deviance.maxpost.train   <- sum(dev.resid.train^2, na.rm = T)
    # Equivalent to deviance calculated for iSDMs under cross-validation:
    # sum(-2*log(ifelse(isdm$mdata$Gammaridae == 1, m$fitted.values, 1-m$fitted.values)), na.rm = TRUE)
    
    # Calculate the null ("primitive") deviance
    train.deviance.null <- -2*sum(train.y*log(train.p.null)+(1-train.y)*log(1-train.p.null), na.rm = T)
    
    # Calculate the sum of squared residual deviance for all taxa
    train.deviance.resid.taxa <- apply(train.deviance.resid^2,2,sum, na.rm=T) # NAs removed
    train.deviance.null.taxa <- -2*apply(train.y*log(train.p.null)+(1-train.y)*log(1-train.p.null),2,sum,na.rm=T)
    
    train.d <- 1-train.deviance.resid.taxa/train.deviance.null.taxa ###
    
    # Standardize deviance by number of observations in training data
    train.std.deviance <- train.deviance.resid.taxa/train.n.samples
    
    # tidy deviance data
    train.deviance <- tibble(null.deviance = train.deviance.null.taxa, 
                                 residual.deviance = train.deviance.resid.taxa, 
                                 std.deviance = train.std.deviance, 
                                 D2 = train.d,
                                 Type = "Training",
                                 Fold = fold,
                                 Model = "jSDM",
                                 Trial = folder,
                                 stringsAsFactors = F)
    
    train.deviance$Taxon <- names(train.std.deviance)
    train.deviance$n.samples <- train.n.samples[train.deviance$Taxon]
    train.deviance$n.present <- train.n.present[train.deviance$Taxon]
    
    train.deviance <- train.deviance[, c("Taxon", "Type", "Fold", "Model", "Trial", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]
    
    # tidy probability data
    train.obs <- as.tibble(occur.taxa)
    train.obs$SiteId <- sites
    train.obs$SampId <- samples
    train.obs <- gather(train.obs, Taxon, Obs, -SiteId, -SampId)
    
    train.p <- as.tibble(p.maxpost, stringsAsFactors = FALSE)
    colnames(train.p) <- colnames(occur.taxa)
    train.p$SiteId <- sites
    train.p$SampId <- samples
    train.p <- gather(train.p, Taxon, Pred, -SiteId, -SampId)
    train.p <- left_join(train.p, train.obs, by = c("SiteId", "SampId", "Taxon"))
    train.p$Type <- "Training"
    train.p$Fold <- fold
    train.p$Model <- "jSDM"
    train.p$Trial <- folder
    
    ### TEST results
    # Drop taxa that do not occur in training data
    test.y <- test.obs[, c("SiteId", "SampId", names(train.n.present))]
    # Join predictors to test observations
    test.predictors <- left_join(test.y[, c("SiteId", "SampId")], predictors, by=c("SiteId", "SampId"))
    
    # Drop rows in y.test and predictors.test where the latter has any NAs
    ind <- !apply(is.na(test.predictors[,inf.fact]),1,FUN=any)
    ind <- ifelse(is.na(ind),FALSE,ind)
    test.predictors <- test.predictors[ind, ]
    test.y <- test.y[ind, ]
    
    # Keep the test sites
    test.sites <- test.predictors$SiteId
    test.samples <- test.predictors$SampId
    names(test.sites) <- test.samples
    
    # Drop the test sites and keep the inputs
    test.inputs <- as.matrix(test.predictors[, inf.fact])
    
    # Calculate predicted probabilities based on fitted parameters and test data 
    # Check if site effects AND latent variables are disabled (FF0)
    if (!extensions$site.effects & extensions$n.latent==0){
      z <- matrix(rep(alpha.taxa.maxpost,nrow(test.inputs)),nrow=nrow(test.inputs),byrow=TRUE) + 
        test.inputs%*%beta.taxa.maxpost 
    }
    
    # Check for site effect and sample the site effect for new sites
    if (extensions$site.effects){
      # For new sites, sample from a normal distribution with mean zero and standard deviation equal to the maximum posterior of the standard deviation of the site effect (gamma)
      sigma.gamma <- res.extracted[["sigma_gamma"]]
      sigma.gamma.maxpost <- sigma.gamma[ind.maxpost]
      
      # Note that site effect is randomly sampled only once for each unique site and then duplicated across the samples
      # of the test data
      test.gamma.sites <- rnorm(length(unique(test.sites)), 0, sigma.gamma.maxpost)
      names(test.gamma.sites) <- unique(test.sites)
      
      test.gamma <- test.gamma.sites[test.sites]
      names(test.gamma) <- test.samples
    }
    
    # Check if site effects are enabled AND latent variables are disabled (TT0)
    if (extensions$site.effects & !extensions$n.latent){
      z <- matrix(rep(alpha.taxa.maxpost,nrow(test.inputs)),nrow=nrow(test.inputs),byrow=TRUE) +
        test.inputs%*%beta.taxa.maxpost + test.gamma
    }
    
    # Check if site effects are enabled AND latent variables are enabled (TT1-TTx)
    if (extensions$site.effects & extensions$n.latent){
      # Calculate linear predictor for single latent variable
      # if n.latent == 1, then residual covariance matrix = beta_lat * t(beta_lat)
      # if n.latent > 1, then residual covariance matrix = t(beta_lat) * beta_lat
      
      if (extensions$n.latent==1){
        lv.sigma <- beta.lat.maxpost * t(beta.lat.maxpost)
        # Check if latent variable is site- or sample-specific, sample accordingly
        if (extensions$lat.site==1){
          test.lv.sites <- rnorm(length(unique(test.sites)), mean = 0, sd = lv.sigma)
          names(test.lv.sites) <- unique(test.lv.sites)
          test.lv <- test.lv.sites[test.sites]
        }
        else if (extensions$lat.site==0){
          test.lv <- rnorm(length(test.samples), mean = 0, sd = lv.sigma)
        }
      } 
      # Calculate linear predictor for multiple latent variables
      else if (extensions$n.latent > 1){
        lv.sigma <- t(beta.lat.maxpost) %*% beta.lat.maxpost
        # Check if latent variables are site- or sample-specific
        if (extensions$lat.site==1){
          test.lv.sites <- rmvnorm(length(unique(test.sites)), mean = rep(0, nrow(lv.sigma)), sigma = lv.sigma)
          rownames(test.lv.sites) <- unique(test.sites)
          test.lv <- test.lv.sites[test.sites,]
        } 
        else if (extensions$lat.site==0){
          test.lv <- rmvnorm(length(test.samples), mean = rep(0, nrow(lv.sigma)), sigma = lv.sigma)
        }
      }
      # Calculate linear predictor for all terms (fixed+site+latent terms)
      z <- matrix(rep(alpha.taxa.maxpost,nrow(test.inputs)),nrow=nrow(test.inputs),byrow=TRUE) +
        test.inputs%*%beta.taxa.maxpost + test.gamma + test.lv
    }
    
    p.maxpost <- 1/(1+exp(-z))
    
    # Keep only taxa present in training data
    test.y <- test.y[, names(train.n.present)]
    test.n.present <- apply(test.y, 2, sum, na.rm = TRUE)
    test.n.samples <- apply(test.y, 2, function(j){sum(!is.na(j))})
    
    # Calculate null predicted probabilities
    test.p.null <- apply(test.y, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
    test.p.null <- ifelse(test.p.null==0,1e-4,test.p.null)
    test.p.null <- matrix(rep(test.p.null,nrow(test.y)),nrow=nrow(test.y),byrow=TRUE)
    
    # Calculate null and residual deviances
    test.deviance.resid <- sign(test.y-0.5)*sqrt(-2*(test.y*log(p.maxpost)+(1-test.y)*log(1-p.maxpost)))
    # deviance.maxpost.test   <- sum(test.deviance.resid^2, na.rm = T)
    test.deviance.null <- -2*sum(test.y*log(test.p.null)+(1-test.y)*log(1-test.p.null), na.rm = T)
    
    test.deviance.resid.taxa <- apply(test.deviance.resid^2,2,sum, na.rm=T) # NAs removed
    test.deviance.null.taxa <- -2*apply(test.y*log(test.p.null)+(1-test.y)*log(1-test.p.null),2,sum,na.rm=T)
    
    # Calculate explanatory power: D2
    test.d <- 1-test.deviance.resid.taxa/test.deviance.null.taxa
    
    # Relative deviance
    test.std.deviance <- test.deviance.resid.taxa/test.n.samples
    
    # Tidy deviance data
    test.deviance <- tibble(null.deviance = test.deviance.null.taxa, 
                                residual.deviance = test.deviance.resid.taxa, 
                                std.deviance = test.std.deviance, 
                                D2 = test.d, 
                                Type = "Testing",
                                Fold = fold,
                                Model = "jSDM",
                                Trial = folder,
                                stringsAsFactors = F)
    
    test.deviance$Taxon <- names(test.std.deviance)
    test.deviance$n.samples <- test.n.samples[test.deviance$Taxon]
    test.deviance$n.present <- test.n.present[test.deviance$Taxon]
    
    test.deviance <- test.deviance[, c("Taxon", "Type", "Fold", "Model", "Trial", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]
    
    # tidy the probability data
    test.obs <- as.tibble(test.obs)
    test.obs <- gather(test.obs, Taxon, Obs, -SiteId, -SampId)
    
    test.p <- as.tibble(p.maxpost, stringsAsFactors = FALSE)
    colnames(test.p) <- colnames(occur.taxa)
    test.p$SiteId <- test.sites
    test.p$SampId <- test.samples
    test.p <- gather(test.p, Taxon, Pred, -SiteId, -SampId)
    test.p <- left_join(test.p, test.obs, by = c("SiteId", "SampId", "Taxon"))
    test.p$Type <- "Testing"
    test.p$Fold <- fold
    test.p$Model <- "jSDM"
    test.p$Trial <- folder
    
    ### Bind the k-fold results
    output$deviance <- bind_rows(output$deviance, train.deviance)
    output$deviance <- bind_rows(output$deviance, test.deviance)
    output$probability <- bind_rows(output$probability, train.p, test.p)
    
    cat('Trial / fold: ', folder, "/", fold, '\n')
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
      test.obs <- test.obs[ind, ]
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
# Identifies  extensions using only very specific filename syntax (use with caution!)
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
  
  # Check for latent variables and number of latent variables
  latvar <- strsplit(extensions[[1]][3], "latvar")
  if ("no" %in% latvar){
    n.latent <- 0
  } else  if (latvar[[1]][1] > 0){
    n.latent <- as.integer(latvar[[1]][1])
  }
  
  # Check if latent variables are by site or by sample
  if ("samp" %in% latvar[[1]][2]){
    lat.site <- 0
  } else {
    lat.site <- 1
  }
  
  model.info <- list("workspace" = workspace,
                     "correlations" = correlations,
                     "site.effects" = site.effect,
                     "n.latent" = n.latent,
                     "lat.site" = lat.site)
  return(model.info)
}

# Extracts results of jSDM from folder containing workspace image (identified with identify.jsdm())
# into one output (list of tidy data tables and model results)
extract.jsdm <- function(directory = 'outputs/paper 1 extensions', folder) {
  # Store multiple model data results
  output <- list("deviance" = data.table(), "probability" = data.table(), "parameters" = data.table())
  
  cat("Load image: ", folder, " -> | ")
  
  full.path <- paste(getwd(), directory, folder, sep="/")
  extensions <- identify.jsdm(full.path)
  model.image <- new.env()
  load(paste(full.path, extensions$workspace, sep = "/"), envir = model.image)
  
  cat("Extract -> | ")
  # Stanfit object
  res <- model.image$res
  
  # Peter's code START:
  # correct chains for multiple local maxima of latent variables:
  res.orig <- res
  res.extracted.trace1 <- extract(res,permuted=FALSE,inc_warmup=FALSE)
  
  n.latent <- extensions$n.latent
  n.chain <- model.image$n.chain
  if ( extensions$n.latent > 0 ){
    name.x        <- "x_lat"
    name.beta     <- "beta_lat"
    parnames      <- names(res@sim$samples[[1]])
    ind.x         <- which(substring(parnames,1,nchar(name.x))==name.x)
    ind.beta      <- which(substring(parnames,1,nchar(name.beta))==name.beta)
    ind.notwarmup <- round(res@sim$warmup/res@sim$thin+1):round(res@sim$iter/res@sim$thin)
    print("")
    print(paste("n.latent",n.latent))
    if ( n.latent == 1 )
    {
      n.x    <- length(ind.x)
      
      # calculate means of chain 1 at all sites
      x.mean.1 <- rep(NA,length=n.x)
      for ( i in 1:n.x ) x.mean.1[i] <- mean(res@sim$samples[[1]][[ind.x[i]]][ind.notwarmup])
      if ( n.chain > 1 )  # adapt other chains to first:
      {
        for ( ch in 2:n.chain )  # match chains to first chain
        {
          print(paste("chain",ch))
          # calculate means of chain ch at all sites
          x.mean.ch <- rep(NA,length=n.x)
          for ( i in 1:n.x ) x.mean.ch[i] <- mean(res@sim$samples[[ch]][[ind.x[i]]][ind.notwarmup])
          dx.min.plus  <- sum(abs(x.mean.1-x.mean.ch))
          dx.min.minus <- sum(abs(x.mean.1+x.mean.ch))
          if ( dx.min.minus < dx.min.plus )
          {
            # change sign:
            x.mean.ch <- -x.mean.ch
            res@sim$samples[[ch]][[ind.x]]    <- -res@sim$samples[[ch]][[ind.x]]
            res@sim$samples[[ch]][[ind.beta]] <- -res@sim$samples[[ch]][[ind.beta]]
            print("sign changed")
          }
        }
      }
    }
    else  # n.latent > 1
    {
      n.x.i    <- round(length(ind.x)/n.latent)
      ind.x    <- matrix(ind.x,ncol=n.latent,byrow=FALSE)    # reformat ind.x
      ind.beta <- matrix(ind.beta,ncol=n.latent,byrow=TRUE)  # reformat ind.beta
      
      # calculate means of chain 1 at all sites and for all latent variables
      x.mean.1 <- matrix(NA,nrow=n.x.i,ncol=n.latent)
      for ( i in 1:n.x.i )
      {
        for ( k in 1:n.latent ) x.mean.1[i,k] <- mean(res@sim$samples[[1]][[ind.x[i,k]]][ind.notwarmup])
      }
      
      if ( n.chain > 1 )
      {
        for ( ch in 2:n.chain )  # match chains to first chain
        {
          print(paste("chain",ch))
          # calculate means of chain ch at all sites and for all latent variables
          x.mean.ch <- matrix(NA,nrow=n.x.i,ncol=n.latent)
          for ( i in 1:n.x.i )
          {
            for ( k in 1:n.latent ) x.mean.ch[i,k] <- mean(res@sim$samples[[ch]][[ind.x[i,k]]][ind.notwarmup])
          }
          for ( k in 1:n.latent )  # match latent variables and signs
          {
            k.min <- k
            plus.min <- TRUE 
            dx.min.plus  <- sum(abs(x.mean.1[,k]-x.mean.ch[,k]))
            dx.min.minus <- sum(abs(x.mean.1[,k]+x.mean.ch[,k]))
            dx.min <- dx.min.plus
            if ( dx.min.minus < dx.min.plus ) { plus.min <- FALSE; dx.min <- dx.min.minus }
            if ( n.latent > k )
            {
              for ( kk in (k+1):n.latent )
              {
                dx.min.plus  <- sum(abs(x.mean.1[,kk]-x.mean.ch[,kk]))
                dx.min.minus <- sum(abs(x.mean.1[,kk]+x.mean.ch[,kk]))
                if ( dx.min.plus  < dx.min ) { dx.min <- dx.min.plus;  plus.min <- TRUE; k.min <- kk }
                if ( dx.min.minus < dx.min ) { dx.min <- dx.min.minus; plus.min <- FALSE; k.min <- kk }
              }
            }
            print(paste("k",k))
            print(paste("k.min",k.min))
            # exchange variables:
            tmp <- x.mean.ch[,k.min]
            x.mean.ch[,k.min] <- x.mean.ch[,k]
            if ( plus.min ) x.mean.ch[,k] <-  tmp
            else            x.mean.ch[,k] <- -tmp
            for ( i in 1:nrow(ind.x) )
            {
              tmp <- res@sim$samples[[ch]][[ind.x[i,k.min]]]
              res@sim$samples[[ch]][[ind.x[i,k.min]]] <- res@sim$samples[[ch]][[ind.x[i,k]]]
              if ( plus.min ) res@sim$samples[[ch]][[ind.x[i,k]]] <-  tmp
              else            res@sim$samples[[ch]][[ind.x[i,k]]] <- -tmp
            }
            print("exchange x completed")
            for ( j in 1:nrow(ind.beta) )
            {
              tmp <- res@sim$samples[[ch]][[ind.beta[j,k.min]]]
              res@sim$samples[[ch]][[ind.beta[j,k.min]]] <- res@sim$samples[[ch]][[ind.beta[j,k]]]
              if ( plus.min ) res@sim$samples[[ch]][[ind.beta[j,k]]] <-  tmp
              else            res@sim$samples[[ch]][[ind.beta[j,k]]] <- -tmp
            }
            print("exchange beta completed")
          }
        }
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
  
  # Name dimensions of parameters within stanfit object
  dimnames(res.extracted[["beta_taxa"]]) <- list(1:dim(res.extracted[["beta_taxa"]][,,])[1], inf.fact, colnames(occur.taxa))
  colnames(res.extracted[["alpha_taxa"]]) <- colnames(occur.taxa)
  
  colnames(res.extracted[["mu_beta_comm"]]) <- inf.fact
  colnames(res.extracted[["sigma_beta_comm"]]) <- inf.fact
  
  # Extract inputs (x), observations (y), and parameters at maximum posterior
  x <- as.matrix(env.cond[,inf.fact])
  colnames(x) <- inf.fact
  y <- as.matrix(occur.taxa)
  ind.maxpost <- which.max(res.extracted[["lp__"]])
  mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
  sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
  mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
  sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
  alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,]
  beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,]
  
  # Name dimensions of maximum posterior community parameters
  names(mu.beta.comm.maxpost) <- inf.fact
  names(sigma.beta.comm.maxpost) <- inf.fact
  rownames(beta.taxa.maxpost) <- inf.fact
  
  # Get the occurrence frequency to order labels
  n <- occur.freq(y)
  
  # Prepare tidy data from maximum posterior beta_taxa
  beta <- t(beta.taxa.maxpost)
  beta <- data.table(beta, stringsAsFactors = F)
  beta$Taxon <- colnames(occur.taxa)
  beta <- gather(beta, Variable, Parameter, -Taxon)
  
  # Match positions of site-samples to unique sites
  unique.sites <- unique(sites)
  siteIND <- match(sites, unique.sites) 
  
  # Additional extraction of JSDM results based on identify.jsdm() output
  # Check extensions for correlated community parameters
  if (extensions$correlations){
    corrfact.comm <- res.extracted[["corrfact_comm"]]
    dimnames(corrfact.comm)[[1]] <- 1:dim(res.extracted[["corrfact_comm"]][,,])[1]
    dimnames(corrfact.comm)[[2]] <- inf.fact
    dimnames(corrfact.comm)[[3]] <- inf.fact
    
    corrfact.comm.maxpost <- corrfact.comm[ind.maxpost, ,]
    
    iter.eff <- round((res@sim$iter-res@sim$warmup)/res@sim$thin)*res@sim$chains
    n.taxa <- length(n)
    if ( extensions$n.latent > 0 )
    {
      if ( extensions$n.latent == 1 )
      {
        cov.resid  <- array(0,dim=c(iter.eff,n.taxa,n.taxa))
        corr.resid <- array(0,dim=c(iter.eff,n.taxa,n.taxa))
        for ( i in 1:iter.eff ) 
        {
          cov.resid[i,,]  <- res.extracted$beta_lat[i,] %*% t(res.extracted$beta_lat[i,])
          corr.resid[i,,] <- cov2cor(cov.resid[i,,])
        }
      }
      else
      {
        cov.resid  <- array(0,dim=c(iter.eff,n.taxa,n.taxa))
        corr.resid <- array(0,dim=c(iter.eff,n.taxa,n.taxa))
        for ( i in 1:iter.eff ) 
        {
          cov.resid[i,,]  <- t(res.extracted$beta_lat[i,,]) %*% res.extracted$beta_lat[i,,]
          corr.resid[i,,] <- cov2cor(cov.resid[i,,])
        }
      }
      cov.resid.maxpost  <- cov.resid[ind.maxpost,,]
      corr.resid.maxpost <- corr.resid[ind.maxpost,,]
    }
    rm(n.taxa)
  }
  
  # Check extensions for site effects
  if (extensions$site.effects){
    # Get named vector of site effects
    gamma.maxpost <- res.extracted[["gamma_site"]][ind.maxpost,]
    names(gamma.maxpost) <- unique.sites
    
    # Get named vectors of random effects applied to the samples
    gamma.samples <- gamma.maxpost[sites]
  }
  
  # Check extensions for latent variables
  if (extensions$n.latent){
    x.lat <- res.extracted[["x_lat"]]
    beta.lat <- res.extracted[["beta_lat"]]
    
    # Check for single latent variable
    if (extensions$n.latent==1){
      x.lat.maxpost <- x.lat[ind.maxpost,]
      beta.lat.maxpost <- beta.lat[ind.maxpost,]
      
      rownames(beta.lat) <- 1:nrow(beta.lat)
      colnames(beta.lat) <- colnames(occur.taxa)
      
      names(beta.lat.maxpost) <- colnames(occur.taxa)
      
      # Check for single site- or sample-specific latent variable
      if (extensions$lat.site==1){
        rownames(x.lat) <- 1:nrow(x.lat)
        colnames(x.lat) <- unique.sites
        
        names(x.lat.maxpost) <- unique.sites
      } else if (extensions$lat.site==0){
        rownames(x.lat) <- 1:nrow(x.lat)
        colnames(x.lat) <- samples
        
        names(x.lat.maxpost) <- samples
      }
      # Check for multiple latent variables  
    } else if (extensions$n.latent > 1){
      x.lat.maxpost <- x.lat[ind.maxpost,,]
      beta.lat.maxpost <- beta.lat[ind.maxpost,,]
      
      dimnames(beta.lat) <- list(1:dim(res.extracted[["beta_lat"]][,,])[1], 1:extensions$n.latent, colnames(occur.taxa))
      
      rownames(beta.lat.maxpost) <- 1:extensions$n.latent
      colnames(beta.lat.maxpost) <- colnames(occur.taxa)
      
      # Check for site- or sample-specific latent variables
      if (extensions$lat.site==1){ # Sites
        dimnames(x.lat) <- list(1:dim(res.extracted[["x_lat"]][,,])[1], unique.sites, 1:extensions$n.latent)
        
        rownames(x.lat.maxpost) <- unique.sites
        colnames(x.lat.maxpost) <- 1:extensions$n.latent
      } else if (extensions$lat.site==0){ # Samples
        dimnames(x.lat) <- list(1:dim(res.extracted[["x_lat"]][,,])[1], samples, 1:extensions$n.latent)
        
        rownames(x.lat.maxpost) <- samples
        colnames(x.lat.maxpost) <- 1:extensions$n.latent
      }
    }
  }
  
  # Calculate model outcomes and statistics based on extracted objects and optional extensions
  cat("Process -> | ")
  
  ### Probability
  cat("Calculate statistics -> | ")
  # Check if site effects AND latent variables are disabled (FF0)
  if (!extensions$site.effects & extensions$n.latent==0){
    z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) + 
      x%*%beta.taxa.maxpost
  }
  
  # Check if site effects are enabled AND latent variables are disabled (TT0)
  if (extensions$site.effects & !extensions$n.latent){
    z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
      x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE)
  }
  
  # Check if site effects are enabled AND latent variables are enabled (TT1-TTx)
  if (extensions$site.effects & extensions$n.latent){
    # Calculate linear predictor for single latent variable
    if (extensions$n.latent==1){
      # Check if latent variable is site- or sample-specific
      if (extensions$lat.site==1){
        # Duplicate the site-specific LV to the samples to calculate probabilities
        lv <- sapply(beta.lat.maxpost, function(j){
          j * x.lat.maxpost[siteIND]
        })
      } 
      else if (extensions$lat.site==0){
        lv <- sapply(beta.lat.maxpost, function(j){
          j * x.lat.maxpost
        })
      }
      # Calculate linear predictor for all terms (fixed+site+latent terms)
      z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
        x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE) + lv
    } 
    # Calculate linear predictor for multiple latent variables
    else if (extensions$n.latent > 1){
      # Check if latent variables are site- or sample-specific
      if (extensions$lat.site==1){
        x.lat.maxpost.samples <- x.lat.maxpost[siteIND,]
        beta.lat.t <- t(beta.lat.maxpost)
        
        # Prepare x.lat*beta.lat as a 3-dimensional array
        lv <- lapply(rownames(beta.lat.t), function(j){
          row <- beta.lat.t[j,]
          t(apply(x.lat.maxpost.samples, 1, function(i){
            row * i
          }))
        })
        # Bind list of matrices into array
        lv <- simplify2array(lv)
        rm(x.lat.maxpost.samples, beta.lat.t)
      } 
      else if (extensions$lat.site==0){
        beta.lat.t <- t(beta.lat.maxpost)
        
        # Prepare x.lat*beta.lat as a 3-dimensional array
        lv <- lapply(rownames(beta.lat.t), function(j){
          row <- beta.lat.t[j,]
          t(apply(x.lat.maxpost, 1, function(i){
            row * i
          }))
        })
        # Bind list of matrices into array
        lv <- simplify2array(lv)
        rm(beta.lat.t)
      }
      
      # Latent terms are summed by dimensions 1 and 3 (i.e., by "sites" or samples and taxa)
      z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) +
        x%*%beta.taxa.maxpost + matrix(rep(gamma.samples,ncol(y)),nrow=nrow(x),byrow=FALSE) + apply(lv, c(1,3), sum)
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
  # Residual deviance
  deviance.maxpost.taxa <- apply(dev.resid^2,2,sum, na.rm=T) # NAs removed
  # Null deviance
  deviance.primitive.taxa <- -2*apply(y*log(p.primitive)+(1-y)*log(1-p.primitive),2,sum,na.rm=T)
  # D^2 statistic
  deviance.fit <- 1-deviance.maxpost.taxa/deviance.primitive.taxa ### D2
  # Number of samples
  n.samples <- apply(y, 2, function(j){ # Number of presence-absence observations
    sum(!is.na(j))
  })
  
  deviance <- data.table(Taxon = names(deviance.fit), Model = "jSDM", null.dev = deviance.primitive.taxa, res.dev = deviance.maxpost.taxa, std.res.dev = deviance.maxpost.taxa/n.samples, D2 = deviance.fit, n.samples = n.samples[names(deviance.fit)], n = n[names(deviance.fit)], stringsAsFactors = F)
  
  row.names(deviance) <- 1:nrow(deviance)
  
  ### Prepare tidy data for output
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
  
  tjur <- extract.tjur(probability)
  tjur$Model <- folder
  
  ### Prepare output
  probability$Model <- "jSDM"
  deviance$Model <- "jSDM"
  beta$Model <- "jSDM"
  
  probability$Trial <- folder
  deviance$Trial <- folder
  beta$Trial <- folder
  
  cat("Save -> || ")
  # Output results to a list
  output$extensions <- extensions
  # output[[folder]]$jSDM <- TRUE # flag for select.jsdm() to identify list containing the res.extracted objects
  
  if (extensions$site.effects){
    output$gamma.maxpost <- gamma.maxpost
    output$gamma.samples <- gamma.samples
  }
  
  if (extensions$correlations){
    output$corrfact.comm <- corrfact.comm
    output$corrfact.comm.maxpost <- corrfact.comm.maxpost
    
  }
  
  if (extensions$n.latent > 0){
    output$x.lat <- x.lat
    output$x.lat.maxpost <- x.lat.maxpost
    output$beta.lat <- beta.lat
    output$beta.lat.maxpost <- beta.lat.maxpost
    output$xb.lat.samples <- lv
    output$corr.resid.maxpost <- corr.resid.maxpost
  }
  
  # Store tidy data
  output$deviance <- bind_rows(output$deviance, deviance)
  output$deviance.residual   <- deviance.maxpost
  output$deviance.null <- deviance.primitive
  
  output$probability <- bind_rows(output$probability, probability)
  output$parameters <- bind_rows(output$parameters, beta)
  output$tjur <- tjur
  # Store the site, samples, community, and input data
  output$sites <- sites
  output$samples <- samples
  output$occur.taxa <- occur.taxa
  output$inf.fact <- inf.fact
  output$env.cond <- env.cond
  
  # Store priors for community parameters
  output$mu.alpha.comm.pripar <- model.image$data$mu_alpha_comm_pripar
  output$sigma.alpha.comm.pripar <- model.image$data$sigma_alpha_comm_pripar
  output$mu.beta.comm.pripar <- model.image$data$mu_beta_comm_pripar
  output$sigma.beta.comm.pripar <- model.image$data$sigma_beta_comm_pripar
  
  # Store community parameters
  output$mu.alpha.comm <- res.extracted[["mu_alpha_comm"]]
  output$mu.alpha.comm.maxpost <- mu.alpha.comm.maxpost
  
  output$sigma.alpha.comm <- res.extracted[["sigma_alpha_comm"]]
  output$sigma.alpha.comm.maxpost <- sigma.alpha.comm.maxpost
  
  output$mu.beta.comm <- res.extracted[["mu_beta_comm"]]
  output$mu.beta.comm.maxpost <- mu.beta.comm.maxpost
  
  output$sigma.beta.comm <- res.extracted[["sigma_beta_comm"]]
  output$sigma.beta.comm.maxpost <- sigma.beta.comm.maxpost
  
  # Store posterior taxon-specific parameters
  output$alpha.taxa <- res.extracted[["alpha_taxa"]]
  output$alpha.taxa.maxpost <- alpha.taxa.maxpost
  
  output$beta.taxa  <- res.extracted[["beta_taxa"]]
  output$beta.taxa.maxpost <- beta.taxa.maxpost
  
  cat("DONE","\n")
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
propogate.jsdm.pred <- function(jsdm, get.quantiles=TRUE, quantiles=c(0.05, 0.95)){
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

# Extract significant responses from jSDM based on posterior taxon-specific parameters and 5% probabilities
extract.resp <- function(jsdm){
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
linear.predictor <- function(jsdm){
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
        bjk <- jsdm$parameters[Taxon==j & Variable %in% variable,]
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
        bjk <- jsdm$parameters[Taxon==j & Variable==k,]
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
  if(jsdm$extensions$site.effects==1){
    gamma.maxpost.samples <- jsdm$gamma.samples
    
    dt <- data.table(SiteId = jsdm$sites, SampId = jsdm$samples, z = gamma.maxpost.samples, x = NA, Taxon = NA, Variable = "Site effect")
    K[[length(K)+1]] <- dt
  }
  
  # If a single latent variable is present, extract the linear predictor
  if(jsdm$extensions$n.latent==1){
    x.lat.maxpost <- jsdm$x.lat.maxpost
    x.lat.maxpost.samples <- x.lat.maxpost[jsdm$sites]
    
    dt <- as.tibble(jsdm$xb.lat.samples)
    
    dt$SiteId <- jsdm$sites
    dt$SampId <- jsdm$samples
    dt$x <- x.lat.maxpost.samples
    dt$Variable <- "TT1"
    dt <- gather(dt, Taxon, z, -SiteId, -SampId, -x, -Variable)
    dt <- select(dt, SiteId, SampId, z, x, Taxon, Variable)
    
    K[[length(K)+1]] <- dt
  } 
  # If multiple latent variables are present, extract their linear predictors
  else if(jsdm$extensions$n.latent > 1){
    # Extract the z-values directly from the jSDM result
    lv <- jsdm$xb.lat.samples
    
    # Get a list of fully named data frames (one per latent variable)
    dt <- lapply(1:jsdm$extensions$n.latent, function(k){
      lv.k <- lv[,k,]
      lv.k <- as.tibble(lv.k)
      colnames(lv.k) <- colnames(jsdm$occur.taxa)
      lv.k$SiteId <- jsdm$sites
      lv.k$SampId <- jsdm$samples
      lv.k
    })
    
    # Add a column to each data frame in the list (to specify the latent variable)
    dt <- lapply(1:jsdm$extensions$n.latent, function(k){
      dt[[k]] %>% mutate(Variable = paste("TT",k,sep=""))
    })
    
    # Build a tidy dataset
    dt <- dt %>%
      rbindlist() %>% # Bind list of data frames together
      as.tibble() %>%
      gather(Taxon, z, -SiteId, -SampId, -Variable) %>% # Consolidate the data frame
      mutate(x = NA) %>%
      select(SiteId, SampId, z, x, Taxon, Variable) # reorder the columns
    
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
extract.beta <- function(jsdm){
  # Get beta samples
  beta.samples <- melt(jsdm$beta.taxa) # Transform 3d array into 2d data.table
  
  # Get trial name
  ind <- sapply(jsdm, function(i){
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

extract.tjur <- function(probability){
  mean.prob <- probability %>%
    group_by(Taxon, Obs) %>%
    summarise(MeanPred = mean(Pred, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Obs = ifelse(Obs == 1, "Present", "Absent")) %>%
    spread(Obs, MeanPred) %>%
    mutate(D = Present - Absent)
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
map.jsdm <- function(jsdm, fileName){
  # Get probability/observations, deviance, and occurrence frequency for all taxa
  probability <- left_join(jsdm$probability, inputs$xy, by="SiteId")
  setDT(probability)
  
  deviance <- jsdm$deviance[, c("Taxon", "D2")]
  d <- deviance$D2
  names(d) <- jsdm$Taxon
  
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
  # Get probability/observations, deviance, and occurrence frequency for all taxa
  probability <- left_join(jsdm$probability, inputs$xy, by="SiteId")
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

map.jsdm.pred <- function(jsdm, fileName){
  cat("Propogating posterior quantiles through joint model...\n")
  dt <- extract.jsdm.pred(jsdm, get.quantiles=TRUE) # Get predicted probabilities with default quantiles c(0.05, 0.95)
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

plot.trace <- function(directory = 'outputs/paper 1 extensions', folder){
  cat("Load image: ", folder, " -> | ")
  
  full.path <- paste(getwd(), directory, folder, sep="/")
  extensions <- identify.jsdm(full.path)
  model.image <- new.env()
  load(paste(full.path, extensions$workspace, sep = "/"), envir = model.image)
  
  files <- list.files(full.path)
  workspace <- files[endsWith(files, ".RData")] # failure point: assumes one .RData file
  filename <- gsub(".RData", "", workspace)
  
  cat("Extract -> | ")
  # Stanfit object
  res <- model.image$res
  
  # Peter's code START:
  # correct chains for multiple local maxima of latent variables:
  res.orig <- res
  res.extracted.trace1 <- extract(res,permuted=FALSE,inc_warmup=FALSE)
  
  n.latent <- extensions$n.latent
  n.chain <- model.image$n.chain
  if ( extensions$n.latent > 0 ){
    name.x        <- "x_lat"
    name.beta     <- "beta_lat"
    parnames      <- names(res@sim$samples[[1]])
    ind.x         <- which(substring(parnames,1,nchar(name.x))==name.x)
    ind.beta      <- which(substring(parnames,1,nchar(name.beta))==name.beta)
    ind.notwarmup <- round(res@sim$warmup/res@sim$thin+1):round(res@sim$iter/res@sim$thin)
    print("")
    print(paste("n.latent",n.latent))
    if ( n.latent == 1 )
    {
      n.x    <- length(ind.x)
      
      # calculate means of chain 1 at all sites
      x.mean.1 <- rep(NA,length=n.x)
      for ( i in 1:n.x ) x.mean.1[i] <- mean(res@sim$samples[[1]][[ind.x[i]]][ind.notwarmup])
      if ( n.chain > 1 )  # adapt other chains to first:
      {
        for ( ch in 2:n.chain )  # match chains to first chain
        {
          print(paste("chain",ch))
          # calculate means of chain ch at all sites
          x.mean.ch <- rep(NA,length=n.x)
          for ( i in 1:n.x ) x.mean.ch[i] <- mean(res@sim$samples[[ch]][[ind.x[i]]][ind.notwarmup])
          dx.min.plus  <- sum(abs(x.mean.1-x.mean.ch))
          dx.min.minus <- sum(abs(x.mean.1+x.mean.ch))
          if ( dx.min.minus < dx.min.plus )
          {
            # change sign:
            x.mean.ch <- -x.mean.ch
            res@sim$samples[[ch]][[ind.x]]    <- -res@sim$samples[[ch]][[ind.x]]
            res@sim$samples[[ch]][[ind.beta]] <- -res@sim$samples[[ch]][[ind.beta]]
            print("sign changed")
          }
        }
      }
    }
    else  # n.latent > 1
    {
      n.x.i    <- round(length(ind.x)/n.latent)
      ind.x    <- matrix(ind.x,ncol=n.latent,byrow=FALSE)    # reformat ind.x
      ind.beta <- matrix(ind.beta,ncol=n.latent,byrow=TRUE)  # reformat ind.beta
      
      # calculate means of chain 1 at all sites and for all latent variables
      x.mean.1 <- matrix(NA,nrow=n.x.i,ncol=n.latent)
      for ( i in 1:n.x.i )
      {
        for ( k in 1:n.latent ) x.mean.1[i,k] <- mean(res@sim$samples[[1]][[ind.x[i,k]]][ind.notwarmup])
      }
      
      if ( n.chain > 1 )
      {
        for ( ch in 2:n.chain )  # match chains to first chain
        {
          print(paste("chain",ch))
          # calculate means of chain ch at all sites and for all latent variables
          x.mean.ch <- matrix(NA,nrow=n.x.i,ncol=n.latent)
          for ( i in 1:n.x.i )
          {
            for ( k in 1:n.latent ) x.mean.ch[i,k] <- mean(res@sim$samples[[ch]][[ind.x[i,k]]][ind.notwarmup])
          }
          for ( k in 1:n.latent )  # match latent variables and signs
          {
            k.min <- k
            plus.min <- TRUE 
            dx.min.plus  <- sum(abs(x.mean.1[,k]-x.mean.ch[,k]))
            dx.min.minus <- sum(abs(x.mean.1[,k]+x.mean.ch[,k]))
            dx.min <- dx.min.plus
            if ( dx.min.minus < dx.min.plus ) { plus.min <- FALSE; dx.min <- dx.min.minus }
            if ( n.latent > k )
            {
              for ( kk in (k+1):n.latent )
              {
                dx.min.plus  <- sum(abs(x.mean.1[,kk]-x.mean.ch[,kk]))
                dx.min.minus <- sum(abs(x.mean.1[,kk]+x.mean.ch[,kk]))
                if ( dx.min.plus  < dx.min ) { dx.min <- dx.min.plus;  plus.min <- TRUE; k.min <- kk }
                if ( dx.min.minus < dx.min ) { dx.min <- dx.min.minus; plus.min <- FALSE; k.min <- kk }
              }
            }
            print(paste("k",k))
            print(paste("k.min",k.min))
            # exchange variables:
            tmp <- x.mean.ch[,k.min]
            x.mean.ch[,k.min] <- x.mean.ch[,k]
            if ( plus.min ) x.mean.ch[,k] <-  tmp
            else            x.mean.ch[,k] <- -tmp
            for ( i in 1:nrow(ind.x) )
            {
              tmp <- res@sim$samples[[ch]][[ind.x[i,k.min]]]
              res@sim$samples[[ch]][[ind.x[i,k.min]]] <- res@sim$samples[[ch]][[ind.x[i,k]]]
              if ( plus.min ) res@sim$samples[[ch]][[ind.x[i,k]]] <-  tmp
              else            res@sim$samples[[ch]][[ind.x[i,k]]] <- -tmp
            }
            print("exchange x completed")
            for ( j in 1:nrow(ind.beta) )
            {
              tmp <- res@sim$samples[[ch]][[ind.beta[j,k.min]]]
              res@sim$samples[[ch]][[ind.beta[j,k.min]]] <- res@sim$samples[[ch]][[ind.beta[j,k]]]
              if ( plus.min ) res@sim$samples[[ch]][[ind.beta[j,k]]] <-  tmp
              else            res@sim$samples[[ch]][[ind.beta[j,k]]] <- -tmp
            }
            print("exchange beta completed")
          }
        }
      }
    }
  }
  
  # Extract the modified stanfit object
  res.extracted.trace2 <- extract(res,permuted=FALSE,inc_warmup=FALSE)
  res.extracted       <- extract(res,permuted=TRUE,inc_warmup=FALSE)
  # Peter's code END
  
  # traceplot
  
  dims <- dim(res.extracted.trace1)
  pdf(paste(directory,"/traceplot_",filename,"_traces.pdf",sep=""),width=8,height=12)
  par(mfrow=c(traceplot.nrow,traceplot.ncol),mar=c(2,2,2,0.5)+0.2) # c(bottom, left, top, right)
  for ( i in 1:ceiling(length(names(res))/(traceplot.nrow*traceplot.ncol)) )
  {
    start <- (i-1)*traceplot.nrow*traceplot.ncol+1
    end   <- min(start+traceplot.nrow*traceplot.ncol-1,length(names(res)))
    for ( j in start:end )
    {
      plot(numeric(0),numeric(0),type="n",cex.axis=0.8,
           xlim=c(0,dims[1]),ylim=range(res.extracted.trace1[,,j]),xlab="",ylab="",
           main=dimnames(res.extracted.trace1)[[3]][j])
      for ( k in 1:dims[2] ) lines(1:dims[1],res.extracted.trace1[,k,j],col=k)
    }
  }
  dev.off()
  
  # traceplot 2
  
  if ( n.latent > 0 )
  {
    dims <- dim(res.extracted.trace2)
    pdf(paste(directory,"/traceplot_",filename,"_traces2.pdf",sep=""),width=8,height=12)
    par(mfrow=c(traceplot.nrow,traceplot.ncol),mar=c(2,2,2,0.5)+0.2) # c(bottom, left, top, right)
    for ( i in 1:ceiling(length(names(res))/(traceplot.nrow*traceplot.ncol)) )
    {
      start <- (i-1)*traceplot.nrow*traceplot.ncol+1
      end   <- min(start+traceplot.nrow*traceplot.ncol-1,length(names(res)))
      for ( j in start:end )
      {
        plot(numeric(0),numeric(0),type="n",cex.axis=0.8,
             xlim=c(0,dims[1]),ylim=range(res.extracted.trace2[,,j]),xlab="",ylab="",
             main=dimnames(res.extracted.trace2)[[3]][j])
        for ( k in 1:dims[2] ) lines(1:dims[1],res.extracted.trace2[,k,j],col=k)
      }
    }
    dev.off()
  }
  
}
plot.commcorr <- function(jsdm){
  require(ellipse)
  if(jsdm$extensions$correlations){
    mu.beta.maxpost <- jsdm$mu.beta.comm.maxpost
    sigma.beta.maxpost <- jsdm$sigma.beta.comm.maxpost
    alpha.taxa.maxpost <- jsdm$alpha.taxa.maxpost
    beta.taxa.maxpost  <- jsdm$beta.taxa.maxpost
    
    corr.maxpost <- diag(length(mu.beta.maxpost))
    
    # labels <- jsdm.TF0$TF0$inf.fact
    # labels <- replace(labels, labels=="Temp2", ":Temp^2")
    # corr.maxpost <- jsdm$corrfact.comm.maxpost %*% t(jsdm$corrfact.comm.maxpost)
    # colnames(corr.maxpost) <- labels
    # rownames(corr.maxpost) <- labels
    # corrplot(corr.maxpost)
    
    
    cov.maxpost <- diag(sigma.beta.maxpost) %*% corr.maxpost %*% diag(sigma.beta.maxpost)
    
    panel.ellipse <- function(x,y,...){
      # decode dimension, variable index, mean and covariance information from the scatterplot data
      # (encoded with a factor to avoid scale distortions of the plots as all elements are used to 
      # determine the plot extension)
      fact <- x[1]
      n <- x[2]/fact
      ind.x <- round(x[3]/fact)
      ind.y <- round(y[3]/fact)
      mu.x  <- x[4]/fact
      mu.y  <- y[4]/fact
      sd.x  <- x[5]/fact
      sd.y  <- y[5]/fact
      corr  <- x[5+ind.y]/fact
      points(x[-(1:(5+n))],y[-(1:(5+n))],...)
      abline(h=0,v=0)
      e <- ellipse(x=matrix(c(1,corr,corr,1),nrow=2),scale=c(sd.x,sd.y),centre=c(mu.x,mu.y),level=0.9)
      lines(x=e[,1],y=e[,2])
    }
    # encode dimension, variable index, mean and covariance information into the scatterplot data
    # (use a factor to avoid scale distortions of the plots as all elements are used to determine
    # the plot extension)
    fact <- 0.001
    cov.info <- rbind(rep(fact,length(mu.beta.maxpost)),
                      fact*rep(length(mu.beta.maxpost),length(mu.beta.maxpost)),
                      fact*(1:length(mu.beta.maxpost)),
                      fact*mu.beta.maxpost,
                      fact*sigma.beta.maxpost,
                      fact*corr.maxpost)
    plot.data <- rbind(cov.info,t(beta.taxa.maxpost))
    pairs(plot.data,
          pch=19,cex=0.5,main="",
          lower.panel=panel.ellipse,upper.panel=panel.ellipse)
  } else{
    stop("Correlated community parameters not found in joint model.")
  }
}

# Given output from extract.jsdm(), plot community and taxon-specific parameters (and priors)
plot.comm <- function(jsdm, filename){
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
  # Detect whether jSDM argument is a collection of iSDMs or a jSDM,
  # create objects for preparation of inputs 
  if ("models" %in% names(results)){
    K <- results$inf.fact
    y <- data.table(SiteId = results$mdata$SiteId, SampId = results$mdata$SampId)
  }else{
    K <- results$inf.fact
    y <- data.table(SiteId = results$sites, SampId = results$samples)
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
  if ("models" %in% names(jsdm)){
    K <- jsdm$inf.fact
    y <- data.table(SiteId = jsdm$mdata$SiteId, SampId = jsdm$mdata$SampId)
  }else{
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
  prob <- left_join(jsdm$probability, inputs, by=c("SiteId", "SampId"))
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

### DEPRECIATED ###

# # Locates and returns extracted Stan object from joint model results
# select.jsdm <- function(results){
#   ind <- sapply(results, function(i){
#     # class(i)[1]=="list"
#     i["jSDM"]==TRUE
#   })
#   ind <- unlist(ind)
#   image <- results[[names(ind[ind==TRUE])]]
# }

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
