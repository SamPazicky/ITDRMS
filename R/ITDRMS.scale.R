#' ITDRMS.scale
#'
#' Scaling and normalization of the data.
#' @param data Data frame with the loaded and ideally cleaned up mass spec data
#' @param remove.control.outliers Logical. If TRUE (default), the program will remove baseline points that significantly deviate from a straight line fit.
#' @param remove.hightemp.outliers Logical: If TRUE (default is FALSE), the program will remove high temperature data points that significantly deviate from a straight line fit. This can reduce false positives but also true positives. Use with caution. 
#' @param normalization.points Integer: How many points should be used for baseline normalization of high temperature data points? Default is 1, but 2 are recommended for 10 doses and even more for more doses, depending on the thermal behaviour of the studied system.
#' @param normalize.baseline Logical: If TRUE (default), all points of baseline data (lowest temperature) will be used to normalize the baseline and the average position of all points will be 1.
#' @param normalize.selection Vector of length 2: If the data should be normalized only on a subset of data, what column and value should be used? For example, 
#' 'normalize.selection=c("Organism","Pf")' will normalized only based on the data that have value 'Pf' in the column 'Organism'. Default is NULL.
#'
#' @import tidyverse
#' @import magrittr
#' @importFrom gtools mixedsort
#'
#' @return A list with three elements. $data is a data frame with  the scaled data, $all_data is a data.frame with the input and output data and $plot is a ggplot with scaling boxplots.
#' @examples 
#' data_scaled <- ITDRMS.scale(data_cleaned)
#' @export
 


ITDRMS.scale = function(
    data=NULL,
    remove.control.outliers=TRUE,
    remove.hightemp.outliers=FALSE,
    normalization.points=1, # how many points to use
    normalize.baseline=TRUE,
    normalize.selection=NULL # can be c("column","value") to normalize only based on a selection from column
    )
{

  customPlot <- list(
    theme_bw(base_size = 12),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=ggplot2::margin(5,5,5,5, "pt")
    )
  )
  
  
  if(is.null(data)) {
    stop("Please include data")
  } else {
    data=as.data.frame(data)
  }
  
  # calculate average dilution factor
  dils <- names(data) %>% as.numeric() %>% na.omit() %>% .[(.)!=0]%>% sort(decreasing=TRUE)
  dil.factor <- mean(dils[-length(dils)]/dils[-1]) %>% round(2)
  
  controlcond <- grep("^[[:digit:]]*C$", unique(data$condition), value=TRUE)
  conditions <- unique(data$condition)[-grep(controlcond,unique(data$condition))]
  
  cat("Scaling...\n")
  
  # calculate abundance normalization factors
  
  if(is.null(normalize.selection)) {
    
    ab_norm_factors <- data %>% dplyr::select(id,condition, starts_with("Abundance")) %>%
      group_by(condition) %>%
      summarise(across(starts_with("Abundance"), ~ sum(.x, na.rm=TRUE))) %>%
      column_to_rownames("condition") %>%
      dplyr::select(starts_with("Abundance"))
    abd_mean <- ab_norm_factors %>% unlist() %>% unname() %>%
      mean(na.rm=TRUE)
    ab_norm_factors <- abd_mean/ab_norm_factors
    ab_norm_factors <- ab_norm_factors %>% as.data.frame() %>% 
      rownames_to_column(var="condition")
    
  } else {
    
    ab_norm_factors <- data %>% 
      filter(!!sym(normalize.selection[1])==normalize.selection[2]) %>%
      dplyr::select(id,condition, starts_with("Abundance")) %>%
      group_by(condition) %>%
      summarise(across(starts_with("Abundance"), ~ sum(.x, na.rm=TRUE))) %>%
      column_to_rownames("condition") %>%
      dplyr::select(starts_with("Abundance"))
    abd_mean <- ab_norm_factors %>% unlist() %>% unname() %>%
      mean(na.rm=TRUE)
    ab_norm_factors <- abd_mean/ab_norm_factors
    ab_norm_factors <- ab_norm_factors %>% as.data.frame() %>% 
      rownames_to_column(var="condition")
    
  }
  
  # apply the abundance normalization factors to ratio columns
  ratio_columns <- data %>% dplyr::select(starts_with("Abundance")) %>% names() %>% str_remove("Abundance_")
  ratio_data <- data %>% dplyr::select(id,condition,all_of(ratio_columns))
  plotdata <- ratio_data %>%
    group_by(condition) %>%
    summarise(across(all_of(ratio_columns), ~ median(.x, na.rm=TRUE))) %>%
    mutate(set="Pre-abundance adjustment")
  
  ratio_data_abadj <- ratio_data %>% 
    dplyr::select(id,condition) %>%
    left_join(ab_norm_factors,by="condition") %>%
    rename_with(~str_remove(.x,"Abundance_")) %>%
    unite("rowname",c("id","condition"),sep=";") %>%
    column_to_rownames("rowname")
  
  ratio_data_abadj <- ratio_data_abadj*(ratio_data%>%dplyr::select(all_of(ratio_columns))) %>%
    as.data.frame()
  all_data <- data %>% bind_cols(ratio_data_abadj%>%rename_with(~str_c("AbAdj_",.x)))
  
  plotdata <- bind_rows(plotdata,
                        ratio_data_abadj %>% 
                          bind_cols(ratio_data%>%dplyr::select(c(id,condition))) %>%
                          group_by(condition) %>%
                          summarise(across(all_of(ratio_columns), ~ median(.x, na.rm=TRUE))) %>%
                          mutate(set="Post-abundance adjustment")
  )
  
  # normalize blank treatment to 1
  ratio_data_abadj <- ratio_data_abadj %>% dplyr::select(gtools::mixedsort(ratio_columns,decreasing=FALSE))
  norm1columns <- gtools::mixedsort(ratio_columns)[1:normalization.points]
  for(i in 1:nrow(ratio_data_abadj)) {
    ratio_data_abadj[i,] <- ratio_data_abadj[i,]/mean(unlist(ratio_data_abadj[i,norm1columns]))

  }
  
  # calculate normalization factors
  norm_factors <- ratio_data_abadj %>%
    rownames_to_column("rowname") %>%
    separate("rowname",c("id","condition"),sep=";") %>%
    group_by(condition) %>%
    summarise(across(all_of(ratio_columns), ~ median(.x, na.rm=TRUE))) %>%
    column_to_rownames("condition") %>%
    mutate(across(everything(), ~ 1/.x)) %>%
    dplyr::select(gtools::mixedsort(ratio_columns,decreasing=FALSE))
  
  ratio_data_norm <- ratio_data %>% 
    dplyr::select(id,condition) %>%
    left_join(norm_factors%>%rownames_to_column("condition"),by="condition") %>%
    dplyr::select(-c(id,condition))
  ratio_data_norm <- ratio_data_abadj*ratio_data_norm %>%
    as.data.frame()
  
  
  
  # remove control outliers
  if (remove.control.outliers) {
    cat("Removing control outliers...\n")
    controldata <- ratio_data_norm[grep(controlcond,rownames(ratio_data_norm)),]
    after.out.intercepts <- c()
    
    x=names(controldata) %>% as.numeric()
    x_adj = x
    x_adj[1]=x_adj[2]/dil.factor
    
    pb <- txtProgressBar(min=0, max=nrow(controldata), style=3, initial="")
    for (i in 1:nrow(controldata)) {
      protein <- rownames(controldata)[i] %>% str_split_i(";",1)
      y=controldata %>% slice(i) %>% unlist() %>% unname() 
      trialfit <- try(lm(formula = y ~ 1, na.action=na.omit), silent=TRUE)
      excl_fit_Rs <- list()
      for (j in 1:length(x)) {
        y_dif <- y[-j]
        x_dif <- x_adj[-j]
        # y_dif <- y_dif + (1-mean(y_dif))
        exfit <- try(lm(formula = y_dif ~ 1, na.action=na.omit), silent=TRUE)
        if(class(exfit)!="try-error") {
          excl_fit_Rs[[j]] <- abs(exfit$coefficients[1]) #summary(exfit)$sigma
        }
      }
      
      outlier.table <- excl_fit_Rs %>% setNames(x) %>% stack() %>% setNames(c("R","x")) %>% 
        add_row(data.frame(R=abs(trialfit$coefficients[1]),x=NA)) %>%
        mutate(z=scale(R)) 
      outlier <- outlier.table %>%
        filter(abs(z)>=0) %>%
        slice_max(abs(z)) %>% 
        pull(x) %>% .[1] %>% as.character()

      if(length(outlier)>0 ) {
        if(!is.na(outlier)) {
          controldata[i,outlier] <- NA
          after.out.intercepts[protein] <- outlier.table %>% filter(x==outlier) %>% pull(R)
        }
      } else {
        after.out.intercepts[protein] <- outlier.table %>% filter(is.na(x)) %>% pull(R)
      }
      
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    after.out.intercepts <- stack(after.out.intercepts) %>% setNames(c("correction","id"))
    # normalize baseline
    if(normalize.baseline) {
      cat("Normalizing baseline to 1...\n")
      controldata <- controldata %>%
        rownames_to_column("id.condition") %>%
        separate_wider_delim(cols=id.condition,delim=";",names=c("id","condition")) %>%
        left_join(after.out.intercepts,by="id") %>%
        mutate(across(!c(id,condition), ~ .x + 1 - correction)) %>%
        mutate(rowname=paste0(id,";",condition)) %>%
        column_to_rownames("rowname") %>%
        dplyr::select(!c(id,condition,correction))
    }
    
    ratio_data_norm[grep(controlcond,rownames(ratio_data_norm)),] <- controldata
  }
  
  
  # check high-temp conditions - I always remove one outlier and fit lm. If a removal of an outlier causes the
  # curve to completely flatten (slope close to zero) then the entire sigmoid fit depends just on that point thus
  # it is not a confident hit. Then this point should be removed and if it is a first point, second point needs to be
  # scaled to 1
  
  
  if (remove.hightemp.outliers) {
    cat("Removing sigmoid fit outliers...\n")
    hightempdata <- ratio_data_norm[grep(controlcond,rownames(ratio_data_norm), invert=TRUE),]
    
    x=names(hightempdata) %>% as.numeric()
    x_adj = x
    x_adj[1]=x_adj[2]/dil.factor
    
    pb <- txtProgressBar(min=0, max=nrow(hightempdata), style=3, initial="")
    for (i in 1:nrow(hightempdata)) {
      y=hightempdata %>% slice(i) %>% unlist() %>% unname() 
      trialfit <- try(lm(formula = y ~ 1, na.action=na.omit), silent=TRUE)
      excl_fit_Rs <- list()
      for (j in 1:length(x)) {
        y_dif <- y[-j]
        x_dif <- x_adj[-j]
        # y_dif <- y_dif + (1-mean(y_dif))
        exfit <- try(lm(formula = y_dif ~ 1, na.action=na.omit), silent=TRUE) 
        excl_fit_Rs[[j]] <- abs(exfit$coefficients[1]) #summary(exfit)$r.squared
      }
      
      outlier <- excl_fit_Rs %>% setNames(x) %>% stack() %>% setNames(c("R","x")) %>% 
        add_row(data.frame(R=abs(trialfit$coefficients[1]),x=NA)) %>%
        mutate(z=scale(R)) %>%
        filter(abs(z)>=2) %>%
        filter(abs(R-1)<=0.1) %>% pull(x) %>% as.character()
      # outlier <- if(length(outlier)>1){ outlier <- c()}
      # filter(y<(mean(y)-2*sd(y)) | y>(mean(y)+2*sd(y))) %>% filter(y==min(abs(y))) %>% pull(x) %>% as.character()
      
      if(length(outlier)>0 ) {
        if(!is.na(outlier[1])) {
          hightempdata[i,outlier] <- NA
        }
      } 
      
      if(x[1] %in% outlier) {
        hightempdata[i,as.character(x)] <- hightempdata[i,as.character(x)]+1-hightempdata[i,as.character(x[2])]
      }
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    ratio_data_norm[grep(controlcond,rownames(ratio_data_norm), invert=TRUE),] <- hightempdata
  }
 
  all_data <- bind_cols(all_data, ratio_data_norm%>%rename_with(~str_c("Norm_",.x)))
  
  plotdata <- bind_rows(plotdata,
                        ratio_data_norm %>% 
                          bind_cols(ratio_data%>%dplyr::select(c(id,condition))) %>%
                          group_by(condition) %>%
                          summarise(across(all_of(ratio_columns), ~ median(.x, na.rm=TRUE))) %>%
                          mutate(set="Post-normalization adjustment")
  )
  
  plot <- all_data %>% 
    dplyr::select(id,condition,all_of(ratio_columns),contains("_")) %>%
    dplyr::select(!starts_with("Abundance")) %>% 
    rename_with(~str_c("Raw_",.x),all_of(ratio_columns)) %>% 
    pivot_longer(cols=!c(id,condition), names_to=c(".value","Dose"), values_to=c("Raw","AbAdj","Norm"),names_sep="_") %>% 
    pivot_longer(-c(id,condition,Dose), values_to="Fraction_soluble",names_to="Normalization_step") %>% 
    mutate(Dose=factor(Dose, levels=gtools::mixedsort(unique(Dose,decreasing=FALSE)))) %>%
    mutate(Normalization_step=factor(Normalization_step,levels=c("Raw","AbAdj","Norm"))) %>%
    ggplot(aes(x=Dose,y=Fraction_soluble,group=Dose)) +
    facet_wrap(vars(condition,Normalization_step), ncol=3) +
    geom_boxplot(outlier.size=0.001) +
    customPlot + theme(axis.text.x=element_text(angle=90, hjust=1,vjust=0.7)) + 
    scale_y_continuous(trans="log",name="Fraction soluble", labels=function(x) sprintf("%.2f", x), breaks=c(0.01,0.1,1,10,100))
  
  output_data <- all_data %>% dplyr::select(-all_of(ratio_columns)) %>% dplyr::select(intersect(names(.),names(data)), starts_with("Norm")) %>%
    dplyr::select(-starts_with("Abundance")) %>% rename_at(vars(starts_with("Norm")), ~str_remove(.x, "Norm_")) %>% 
    dplyr::select(id,description,condition,all_of(gtools::mixedsort(ratio_columns)),everything())
  
  return(list("data"=output_data,"all_data"=all_data,"plot"=plot))
}
