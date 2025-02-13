#' ITDRMS_sub.fit
ITDRMS_sub.fit = function(
    data=NULL,
    i=1,
    outlier.removal=TRUE
) {
  
  require(tidyverse)
  require(data.table)
  
  dils <- names(data) %>% as.numeric() %>% na.omit() %>% .[(.)!=0] %>% sort(decreasing=TRUE)
  dil.factor <- mean(dils[-length(dils)]/dils[-1]) %>% round(2)
  
  suppressWarnings(
    ratio_columns <- try(names(data) %>% as.numeric() %>% .[!is.na(.)] %>% as.character(),silent=TRUE)
  )
  conds_fitresults <- data.frame() 
  condition <- data %>% slice(i) %>% rownames() %>% str_split(.,pattern=";") %>% unlist %>% .[2]
  id <- data %>% slice(i) %>% rownames() %>% str_split(.,pattern=";") %>% unlist %>% .[1]
  fitname <- data %>% slice(i) %>% rownames()
  fit_data <- data %>% slice(i) %>% t() %>% as.data.frame() %>% rownames_to_column() %>%
    mutate(across(everything(), as.double)) %>% setNames(c("x","y"))
  lowest.x <- fit_data$x[1] <- fit_data$x[2]/dil.factor
  y=fit_data$y
  prev_rem <- which(is.na(y))
  if(length(prev_rem)==0) {
    prev_rem=length(y)+1
  } 
  x=fit_data$x[-prev_rem]
  y=fit_data$y[-prev_rem]
  cur_ratio_columns = ratio_columns[-prev_rem]
  cur_fit_data<-data.frame("x"=x,"y"=y)
  # try fitting LL.4 - for non-control conditions
  if(!str_detect(condition,"[[:digit:]]C$")|length(prev_rem>1)) {
    sigmoid <- try(fit_sigmoid(cur_fit_data),silent=TRUE)
    if(!class(sigmoid)=="list"&!class(sigmoid)=="try-error"&any(!is.na(sigmoid))) {
      R2sigmoid <- 1 - sum((residuals(sigmoid)^2))/sum((y-mean(y))^2)
      R2sigmoid_orig <- R2sigmoid
      if(outlier.removal) {
        outR2s <- vector()
        sigmoids_out <- list()
        for(out in seq_along(x)) {
          x_out <- x[-out]
          y_out <- y[-out]
          out_data <- data.frame(x=x_out,y=y_out)
          sigmoids_out[[out]] <- fit_sigmoid(out_data)
          if(class(sigmoids_out[[out]])=="list"|any(is.na(sigmoids_out[[out]]))) {
            outR2s[out] <- NA
          } else {
            outR2s[out] <- 1 - sum((residuals(sigmoids_out[[out]])^2))/sum((y_out-mean(y_out))^2)
          }
        }
        outR2s <- c(outR2s,R2sigmoid)
        outlier <- which(outR2s>(mean(outR2s,na.rm=TRUE)+2*sd(outR2s,na.rm=TRUE)))
        outlier <- which(outR2s==max(outR2s[outlier]))
        if(length(outlier)>0) {
          if(outlier!=length(outR2s)) {
            data[i,as.character(format(x[outlier], scientific=FALSE,trim=TRUE))] <- NA
            R2sigmoid <- outR2s[outlier]
            sigmoid <- sigmoids_out[[outlier]]
            cur_fit_data <- cur_fit_data %>% slice(-outlier)
          }
        } else {
          outlier=100
        }
      }
    } else {
      sigmoid <- list()
      R2sigmoid = -100
      conf.int.sig=data.frame("x"=fit_data$x,"half"=NA) #rep(NA,10))
      outlier=100
    }
  } else {
    sigmoid <- list()
    R2sigmoid=-99
    outlier=100
  }
  
  if(outlier.removal) {
    cur_ratio_columns <- cur_ratio_columns[-outlier]
  }
  

  if(class(sigmoid)!="list"&!is.na(class(sigmoid))) {
    suppressWarnings(R2 <- 1 - sum((residuals(sigmoid)^2))/sum((cur_fit_data$y-mean(cur_fit_data$y))^2))
    
    suppressWarnings(
      conds_fitresults <- predict(sigmoid, newdata=cur_fit_data%>%dplyr::select(x), interval="confidence") %>% as.data.frame() %>% 
        mutate(half=Upper-Prediction) %>% dplyr::select(Prediction, half) %>% setNames(c("fit","conf.int")) %>% 
        mutate(x=cur_ratio_columns) %>% pivot_longer(cols=!x) %>% unite("rowname",name,x) %>%
        mutate(rowname=str_remove(rowname,"\\.$")) %>% 
        mutate(rowname=factor(rowname,levels=gtools::mixedsort(unique(rowname)))) %>% arrange(rowname) %>% 
        column_to_rownames("rowname") %>% t() %>% as.data.frame() %>%
        mutate("fit"="nls","R2"=R2sigmoid,"R2orig"=R2sigmoid_orig, Slope=coefficients(sigmoid)[1],"EC50"=coefficients(sigmoid)[3])
    )
    
  } else {
    conds_fitresults <- data.frame("rows"=c(cur_ratio_columns,"fit","R2","R2orig","Slope","EC50"), "values"=rep(NA,5+length(cur_ratio_columns))) %>%
      column_to_rownames("rows") %>% setNames(row.names(data)[i]) %>% t() %>% as.data.frame() %>%
      setNames(c(if(length(cur_ratio_columns)==0){NULL} else {paste0("fit_",cur_ratio_columns)}, "fit","R2","R2orig","Slope","EC50"))
  }
  output <- list(
    data=conds_fitresults,
    fit=sigmoid
  )
  return(output)
}
