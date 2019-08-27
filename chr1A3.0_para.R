library(vcfR)
require(xgboost)
library(dplyr) ## for mutate function
library(parallel)
# 

## Function to compute classification error ########################################
classification_error <- function(conf_mat) {
  conf_mat = as.matrix(conf_mat)
  
  error = 1 - sum(diag(conf_mat)) / sum(conf_mat)
  
  return (error)
}

## vcf to gt function ##############################################################
vcf2gt <- function(filename) {
  df <- read.table(filename)
  
  df_gt <- t(df[,10:ncol(df)]) 
  df_gt <- gsub(":.*","",df_gt)
  df_gt <- gsub("\\|","/",df_gt)
  df_gt <- gsub("/"," ",df_gt)
  write.table(df_gt,"in.txt",sep=" ", quote=FALSE,row.names=FALSE,col.names=FALSE)
  df_gt=read.table("in.txt")
  df_gt=as.matrix(df_gt)
  
  class(df_gt) <- "numeric"
  
  table_df_gt <- apply(df_gt, 2, table)
  levels_table_df_gt <- unlist(lapply(table_df_gt, function(x) sum(as.numeric(names(x))>1)))
  #df_gt_sub <- df_gt[,levels_table_df_gt==0]
  check <- unique(ceiling(which(levels_table_df_gt>0)/2))

  n <- ncol(df_gt)
  df_out <- df_gt[,seq(1,n,by=2)] + df_gt[,seq(2,n,by=2)]
  colnames(df_out) <- df$V3
  df_out <- df_out[, -check]
  return(df_out)
}

## Generate xgboost database function ###############################################
Generate_XGBoost_Dataset <- function(df, size) {
  ## To save time, we only impute SNP columns with NAs in the dataset.
  NA_cols <- which(apply(df, 2, function(x) sum(is.na(x))> 0))
  df_sub <- df[,NA_cols]
  df_rest <- df[, -NA_cols]
  N <- ncol(df_sub)
  
  #df_noNA_sub <- df_noNA[, NA_cols]
  set.seed(101)
  XGBoost_dataset <- list()
  for (i in 1:N) {
    XGBoost_dataset[[i]] <- list()
    snp_NA <- which(is.na(df_sub[,i]))
    snp_notNA <- which(!is.na(df_sub[,i]))
    
    samp <- sample(length(snp_notNA), size = length(snp_notNA)*0.2)
    train <- snp_notNA[-samp]
    test <- snp_notNA[samp]
    if (sum(snp_NA) != 0) { ## double check that this column has NAs to impute
      ## Create train and test data
      if (i == 1)  range = c((i+1): size)
      else if (i == N)  range = c((i-size): N)
      else if (i <= size/2) range = c(1:(i-1), (i+1): size)
      else if ((i + size/2) >= N)  range = c((i-size):(i-1), (i+1): N)
      else range = c((i-size/2):(i-1), (i+1): (i + size/2))
      
      train_df <- as.matrix(df_sub[train, range])
      test_df <- as.matrix(df_sub[test, range])
      pred_df <- as.matrix(df_sub[snp_NA, range])
      #pred_df <- as.matrix(df_sub[, range])
      
      if(ncol(test_df)==1) test_df <- t(test_df)
      if(ncol(pred_df)==1) pred_df <- t(pred_df)
      
      ## Create labels
      train_df_label <- as.numeric(df_sub[train, i])
      test_df_label <- as.numeric(df_sub[test, i])
      pred_df_label <- as.numeric(df_sub[snp_NA, i]) # change from df_noNA_sub to df_sub
      
      ## Prepare matrices
      XGBoost_dataset[[i]][[1]] <- xgb.DMatrix(data = train_df, label = train_df_label)
      XGBoost_dataset[[i]][[2]] <- xgb.DMatrix(data = test_df, label = test_df_label)
      XGBoost_dataset[[i]][[3]] <- test_df_label
      XGBoost_dataset[[i]][[4]] <- snp_NA
      XGBoost_dataset[[i]][[5]] <- list()
      XGBoost_dataset[[i]][[6]] <- list()
      XGBoost_dataset[[i]][[7]] <- xgb.DMatrix(data = pred_df, label = pred_df_label)
      XGBoost_dataset[[i]][[8]] <- colnames(df_sub)[i]
      
      names(XGBoost_dataset[[i]]) <- c("xgb_train", "xgb_test", "test_df_label", "snp_NA","pred","error", "xgb_pred","SNP")
    }
  }
  names(XGBoost_dataset) <- colnames(df_sub)
  out <- list()
  out[[1]] <- XGBoost_dataset
  out[[2]] <- df_sub
  out[[3]] <- df_rest
  names(out) <- c("XGBoost_dataset","df_NA","df_noNA")
  return(out)
}

## Imputation function ####################################################################
Impute_GenoType_XGBoost <- function(XGBoost_dataset, test = length(XGBoost_dataset[[1]])) {
  
  params <- list(booster = "gbtree", objective = "multi:softprob", num_class = 3, eval_metric = "mlogloss")
  df_sub <- XGBoost_dataset[[2]]
  ## Model and prediction
  XGBoost_dataset_sub <- mclapply(XGBoost_dataset[[1]][c(1:test)], function(x) {
    xgb_model <- xgb.train(params = params, data =x[[1]], nrounds = 100)
    # Predict for validation set
    xgb_val_tests <- predict(xgb_model, newdata = x[[2]])
    xgb_val_out <- matrix(xgb_val_tests, nrow = 3, ncol = length(xgb_val_tests) / 3) %>% 
      t() %>%
      data.frame() %>%
      mutate(max = max.col(., ties.method = "last"), label = x[[3]]) 
    
    x[[5]] <- xgb_val_out$max
    
    # Confustion Matrix
    xgb_val_conf <- table(true = x[[3]], pred = x[[5]])
    x[[6]] <- classification_error(xgb_val_conf)
    #df_sub[x[[4]], i] <- x[[5]]-1
    
    if (x[[6]] > 0) {
      cat("XGB Validation Classification Error Rate of", x[[8]], " is: ", x[[6]], "\n")
      cat("Test should be: ",x[[3]],"\n")
      cat("But turns out to be: ", x[[5]]-1,"\n")
    }
    else {
      cat("Passed the test. \n")
      cat("XGB Validation Classification Error Rate of", x[[8]], " is: ", x[[6]], "\n")
      cat("Test should be: ",x[[3]],"\n")
      cat("But turns out to be: ", x[[5]]-1,"\n")
      xgb_preds <- predict(xgb_model, newdata = x[[7]])
      xgb_preds_out <- matrix(xgb_preds, nrow = 3, ncol = length(xgb_preds) / 3) %>% 
        t() %>%
        data.frame() %>%
        mutate(max = max.col(., ties.method = "last")) 
      df_sub[x[[4]], x[[8]]] <- xgb_preds_out$max-1
      cat("Prediction is ",df_sub[x[[4]], x[[8]]],"\n")
    }
    cat("\n")
    return(x)
  })
  
  ## Show the results: 
  bind <- cbind(df_sub, XGBoost_dataset[[3]])
  mean_error <- mean(as.numeric(unlist(lapply(XGBoost_dataset_sub, function(x) x[[6]]))))
  cat(paste0("The mean classification error is ", mean_error, ".\n"))
  
  return(df_sub)
}

############## main ######################################################################
                                              
start_time <- Sys.time()
chr1A <- vcf2gt(filename = "chr1A.vcf")
check1_time <- Sys.time() 

cat("Time to convert vcf file is ", check1_time - start_time, "\n")

#chr1A <- read.table("chr1A_gt.txt")
XGBoost_datasets <- Generate_XGBoost_Dataset(chr1A, 100)
#filename <- "chr1A.vcf"
cat("Time to generate xgboost dataset is ", Sys.time() - check1_time, "\n")

check2_time = Sys.time()
#cl <- makeCluster(detectCores())
out <- Impute_GenoType_XGBoost(XGBoost_datasets)
#stopCluster(cl)
end_time <- Sys.time()
cat("Time to build model and finish prediction is ", end_time - check2_time, "\n")

cat(paste0("The total runing time is ", end_time - start_time))

