#' ITDRMS.fit
#'
#' Fits the scaled data.
#' @param data Data frame scaled mass spec data
#' @param ncores Integer: How many cores to use for fitting. Default is 1.
#' @param fit.length Integer: How many points should be used for fitting curves. Default is 100 which is sufficient for plotting.
#' @param control.slope Boolean: Can the control temperature lines be fitted with y=ax+b (TRUE) or just a horizontal line (y=b)?
#' @param outlier.removal Logical: If TRUE (default is TRUE), outliers will be identified as points that significantly worsen the fit and removed.
#' 
#' @import tidyverse
#' @import magrittr
#' @importFrom gtools mixedsort
#' @import drc
#' 
#' @return A list with three elements. $data is a data frame with the fitted data, $fits is a list with fit objects and $curvy_controls is a data frame with those proteins, that show dose-response behaviour at the lowest temperature (baseline has a sigmoid shape).
#' @examples 
#' data_fitted <- ITDRMS.fit(data_scaled)
#' @export

ITDRMS.fit = function(
    data=NULL,
    ncores=1,
    fit.length=100,
    control.slope=FALSE,
    outlier.removal=TRUE)
{

  if(is.null(data)) {
    stop("Please include data")
  } else {
    data=as.data.frame(data)
  }
  
  # calculate average dilution factor
  dils <- names(data) %>% as.numeric() %>% na.omit() %>% .[(.)!=0] %>% sort(decreasing=TRUE)
  dil.factor <- mean(dils[-length(dils)]/dils[-1]) %>% round(2)
  
  
  controlcond <- grep("^[[:digit:]]*C$", unique(data$condition), value=TRUE)
  conditions <- unique(data$condition)[-grep(controlcond,unique(data$condition))]
  
  data <- data  %>%
    mutate(condition=factor(condition,levels=c(controlcond,conditions))) %>%
    arrange(condition,id) %>%
    mutate(condition=as.character(condition))
  
  suppressWarnings(
    ratio_columns <- try(names(data) %>% as.numeric() %>% .[!is.na(.)] %>% as.character(),silent=TRUE)
  )
  
  ratio_data <- data %>% 
    remove_rownames() %>% unite("rowname",id,condition,sep=";") %>% column_to_rownames("rowname") %>%
    dplyr::select(all_of(ratio_columns))
  
  curvy_controls=data.frame()
  
  topconc <- max(as.numeric(ratio_columns))
  divfactor <- mynthroot(topconc/min(as.numeric(ratio_columns)[-which(ratio_columns=="0")]),(fit.length-1))
  fakedata <- data.frame("x"=topconc/divfactor^(0:(fit.length-1)) )
  curves=data.frame()
  # iterate fitting through all data
  cat("\nFitting curves...\n")
  
  
  fitted_data<-data.frame()
  fits=list()
  
  ratio_data_control <- ratio_data[grep(controlcond,row.names(ratio_data)),]
  column_names <- c(paste0("conf.int_", ratio_columns), 
                    paste0("fit_", ratio_columns), 
                    "fit", "R2", "Slope", "EC50")
  
  # Create an empty data frame with these column names
  control_fitresults <- data.frame(matrix(ncol = length(column_names), nrow = 0)) %>%
    setNames(column_names) %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(fit=as.character(fit))
  
  cat("\nFitting control lines... \n")
  
  zz<-file("temp.txt",open="wt")
  sink(zz, type="message")
  pb <- txtProgressBar(min=0, max=nrow(ratio_data_control), style=3, initial="")
  for(i in 1:nrow(ratio_data_control)) {
    
    
    condition <- ratio_data_control %>% slice(i) %>% rownames() %>% str_split(.,pattern=";") %>% unlist %>% .[2]
    id <- ratio_data_control %>% slice(i) %>% rownames() %>% str_split(.,pattern=";") %>% unlist %>% .[1]
    fitname <- ratio_data_control %>% slice(i) %>% rownames()
    fit_data_orig <- ratio_data_control %>% slice(i) %>% t() %>% as.data.frame() %>% rownames_to_column() %>%
      mutate(across(everything(), as.double)) %>% setNames(c("x","y"))
    fit_data <- fit_data_orig
    fit_data$x[1] <- fit_data_orig$x[2]/dil.factor
    y=fit_data$y
    prev_rem <- which(is.na(y))
    if(length(prev_rem)==0) {
      prev_rem=length(y)+1
    } 
    cur_ratio_names <- ratio_columns[-prev_rem]
    
    x=fit_data$x[-prev_rem]
    y=fit_data$y[-prev_rem]
    cur_ratio_columns = ratio_columns[-prev_rem]
    cur_fit_data<-data.frame("x"=x,"y"=y)
    
    if(control.slope) {
      linear <- try(lm(formula = y ~ log(x,dil.factor), na.action=na.omit), silent=TRUE)
    } else {
      linear <- try(lm(formula = y ~ 1, na.action=na.omit), silent=TRUE)
    }
    if(class(linear)!="try-error") {
      R2linear <- summary(linear)$r.squared
      
      #calculate linear prediction intervals
      control_fitresults <- bind_rows(control_fitresults,
                                      predict(linear, newdata=cur_fit_data%>%dplyr::select(x), interval="confidence") %>% as.data.frame() %>% 
                                        mutate(half=upr-fit) %>% dplyr::select(fit, half) %>% setNames(c("fit","conf.int")) %>% 
                                        mutate(x=cur_ratio_names) %>% pivot_longer(cols=!x) %>% unite("rowname",name,x) %>%
                                        mutate(rowname=str_remove(rowname,"\\.$")) %>% 
                                        mutate(rowname=factor(rowname,levels=gtools::mixedsort(unique(rowname)))) %>% arrange(rowname) %>% 
                                        column_to_rownames("rowname") %>% t() %>% as.data.frame() %>%
                                        mutate(fit="lm", "R2"=R2linear,"Slope"=coefficients(linear)[2],"EC50"=NA)
      )
    } 
    
    sigmoid <- try(fit_any_sigmoid(cur_fit_data),silent=TRUE,outFile="zzz.txt")
    if(class(sigmoid)!="try-error"&!is.na(sigmoid[1])) {
      R2sigmoid <- 1 - sum((residuals(sigmoid)^2))/sum((y-mean(y))^2)
      if( R2sigmoid>=0.6 & abs(coefficients(sigmoid)[2]-coefficients(sigmoid)[3])>=0.2 ) {
        curvy_controls <- bind_rows(curvy_controls,
                                    data.frame(
                                      "id.condition"=rownames(ratio_data_control)[i],
                                      "fit"="nls",
                                      "R2"=R2sigmoid,
                                      "Slope"=coefficients(sigmoid)[1],
                                      "EC50"=coefficients(sigmoid)[4]
                                    ) %>% remove_rownames() %>% separate_wider_delim("id.condition",";",names=c("id","condition"))
        )
      }
      
    }
    

    
    fits[[fitname]] <- linear
    
    
    setTxtProgressBar(pb, i)
  }
  sink(NULL, type="message")
  close(pb) # to close the progress bar
  
  ratio_data_conds <- ratio_data[-grep(controlcond,row.names(ratio_data)),]
  
  # cat("Subtracting baselines... \n")
  # pb <- txtProgressBar(min=0, max=nrow(ratio_data_conds), style=3, initial="")
  # 
  # for(i in 1:nrow(ratio_data_conds)) {
  #   protein=str_extract(row.names(ratio_data_conds)[i], "^.*(?=;)")
  #   # subrow <- ratio_data_control[grep(protein,row.names(ratio_data_control)),]
  #   subrow <- control_fitresults[grep(paste0(protein,";"),row.names(ratio_data_control)),] %>% dplyr::select(starts_with("fit")) %>% dplyr::select(!fit)
  #   if(nrow(subrow)==0) {
  #     next
  #   }
  #   ratio_data_conds[i,] = ratio_data_conds[i,] - subrow + 1
  #   setTxtProgressBar(pb, i)
  #   
  # }
  # close(pb)
  
  conds_fitresults_fits <- list()
  cat("Fitting curves with log-logistic function... \n")
  
  pb <- txtProgressBar(min=0, max=nrow(ratio_data_conds), style=3, initial="")
  if(ncores<=1) {
    
    conds_fitresults_fits <- progress_lapply(1:nrow(ratio_data_conds), 
                                        function(xx) ITDRMS_sub.fit(data=ratio_data_conds,i=xx,outlier.removal=outlier.removal),
                                        pb)
  } else {
    require(parabar)
    require(parallel)
    set_option("progress_track", TRUE)
    set_option("stop_forceful", TRUE)
    backend <- start_backend(cores = ncores, cluster_type = "psock", backend_type = "async")
    export(backend, c("pb","ITDRMS_sub.fit","fit_sigmoid", "outlier.removal","ratio_data_conds"), envir = environment())
    
    conds_fitresults_fits <- par_lapply(backend, 1:nrow(ratio_data_conds), function(x)
      ITDRMS_sub.fit(data=ratio_data_conds,i=x,outlier.removal=outlier.removal)
    )
    stop_backend(backend)
  }
  
  names(conds_fitresults_fits) <- rownames(ratio_data_conds)
  conds_fitresults <- lapply(conds_fitresults_fits, function(x) x$data) %>%
    rbindlist(use.names=TRUE,fill=TRUE)
  fits <- c(fits, lapply(conds_fitresults_fits, function(x) x$fit))
  
  fitresults <- bind_rows(control_fitresults,conds_fitresults)
  # if the calculated EC50 is lower than lowest drug concentration, then it might as well just be an outlier in the
  # zero-drug concentration. The experiment should be repeated with lower drug doses.
  
  all_data <- bind_cols(data,fitresults) %>%
    mutate(EC50=ifelse(EC50<as.numeric(ratio_columns[2]), paste0("<",as.numeric(ratio_columns[2])), EC50))
  
  # calculate fit and prediction values at EC50 with propagate package
  
  output <- list("data"=all_data,"fits"=fits, "curvy_controls"=curvy_controls)
  return(output)
}


