#' ITDRMS.hit
#'
#' Identifies hits from fitted mass spec data.
#' @param data Data frame: Scaled data with removed outliers and fitting statistics, ideally $data element from ITDRMS.fit output.
#' @param weighted Logical: If TRUE (default), the resulting adjusted p values will be weighted by R-squared measure of goodness of fit.
#' @param R2w Integer: How much should the R-squared weight in? Default is 10.
#' @param nMAD Integer: How much higher or lower does the area under dose-response curve (dAUC) be than mean absolute deviation (MAD)? Default is 3.
#' @param minresponse Integer: Minimal change in the solubility (e.g. 0.1 for 10% change). Default is NA.
#' @param mindAUC Integer: Minimum absolute dAUC value to be considered a hit. If both nMAD and mindAUC are given, the larger of the two will be used as the cutoff. Default is 0.1
#' @param R2line Double: Soft R2 cut-off. Proteins with R-squared value below this value will not be considered as hits despite favourable dAUC and p-value.
#' @param POI Character vector: ID of the protein or proteins to be highlighted on the plot.
#' @param plot.settings List of graphical settings for plot. Defaults are: 
#' list(labels=TRUE, label.text.size=2.5, label.force=1.3,
#' xlims=c(-max(c(abs(hit_data$dAUC)),2), +max(c(2,abs(hit_data$dAUC)))),
#' ylims=c(min(-log10(hit_data$CI)),max(-log10(hit_data$CI))),
#'  point.sizes=c(1,2), point.colors=c("gray","red","green"),
#'  axis.title.size=18, axis.text.size=16,
#'  legend.position="bottom", legend.text.size=8,legend.title=element_blank(),
#'  POI.color="blue", POI.size=3.5)
#' 
#' @import tidyverse
#' @import magrittr
#' @importFrom gtools mixedsort
#' @importFrom minpack.lm nlsLM
#' @importFrom ggrepel geom_label_repel
#' 
#' @return A list with three elements. $data is a data frame with dAUC and p-values for all proteins, $plot is the volcano plot based on the hits and 
#' $hitlist is a list with two vectors: one for Stabilized and one for Destabilized hits.
#' @examples 
#' data_fitted <- ITDRMS.hit(data_fitted, nMAD=0.2)
#' @export

ITDRMS.hit <- function(
    data=NULL,
    weighted=TRUE,
    R2w=10,
    nMAD=3,
    mindAUC=0.1,
    minresponse=NA,
    R2line=0.6,
    POI=c(),
    plot.settings=list()) 
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
  
  fit_01sigmoid <- function(data) {
    if(length(data$x)!=length(data$y)) {
      stop("Vectors x and y do not have the same length.")
    }
    # guess the initial value of inflection point for fitting. 
    # The guess is made when difference between subsequent scaled abundances is >0.25
    
    
    scaledata=data
    scaledata$x<-range.scale(scaledata$x)
    for (point in 1:(length(scaledata$y)-1)) {
      y_dif <- abs(scaledata$y[point+1] - scaledata$y[point])
      if(abs(y_dif)>0.25) {
        e_guess <- (data$x[point+1]+data$x[point])/2
        e_vec <- point+1
        break
      } else {
        y_dif_saved <- y_dif
        e_guess <- median(data$x)
        e_vec <- 5
      }
    }
    
    #guess slope
    b_guess <- ((data$y[e_vec-1]-data$y[e_vec])*10)^3/7
    b_guess=b_guess*(-1)
    
    # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
    assign("my.fit.dat",
           try(minpack.lm::nlsLM(formula=y~(((1)/(1+exp(b*(log(x)-log(e)))))),
                                 data=data,
                                 start=list(b=b_guess, e=e_guess),
                                 lower=c(-200,0.0005),
                                 upper=c(200,10),
                                 control=list(maxiter=100)),
               silent=TRUE)
    )
    
    return(my.fit.dat)
  }
  
  fitted <- fit_01sigmoid(data.frame(x=c(1,1.5,2),y=c(0.05,0.5,0.95)))  # equation for calculating p-value
  
  #log4 axis transformation
  
  require(scales) # trans_new() is in the scales library
  log4_trans = function() trans_new("log4", function(x) log(x,4), function(x) 4^x)
  
  
  controlcond <- grep("^[[:digit:]]*C$", unique(data$condition), value=TRUE)
  conditions <- unique(data$condition)[-grep(controlcond,unique(data$condition))]
  
  ratio_columns <- data %>% dplyr::select(starts_with("conf.int")) %>% names() %>% str_split(pattern="_") %>%
    as.data.frame() %>% slice(2) %>% unlist() %>% unname()
  
  top.conc <- as.character(max(as.numeric(ratio_columns)))

  hit_data <- data %>%
    mutate(response=!!sym(paste0("fit_",top.conc))) %>%
    mutate(across(all_of(ratio_columns), ~ifelse(is.na(.x), NaN, .x))) %>% # to distinguish missing values (NA) and outliers (NaN)
    dplyr::select(id,condition,R2orig,response,matches(ratio_columns)) %>%
    mutate(R2orig=ifelse(is.na(R2orig), 0, R2orig)) %>%
    mutate(R2=ifelse(condition==controlcond,0,R2orig)) %>%
    rename_with(~ paste0("ratio_",.x),all_of(ratio_columns)) %>%
    # mutate(response=!!sym(paste0("ratio_",top.conc))) %>%
    pivot_longer(cols=!c(id,condition,R2,response), names_sep="_", names_to=c(".value","Dose")) %>%
    filter(Dose!=0) %>%
    
    mutate(conf.int=ifelse(is.na(conf.int)&(!is.nan(ratio)), ratio/2,conf.int)) %>%
    
    filter(if_any(everything(), ~ !is.nan(.x))) %>% # removal of rows with removed outliers (NaN)
    pivot_wider(id_cols=c(id,Dose), names_from=condition, values_from=c(ratio,conf.int,R2,response), names_sep="_") %>%
    mutate(across(starts_with("ratio_"), ~ .x - !!sym(paste0("ratio_",controlcond)), .names = "sub.{col}" )) %>%
    mutate(across(starts_with("conf.int_"), ~ .x + !!sym(paste0("conf.int_",controlcond)), .names="sum.{col}")) %>%
    dplyr::select(!contains(controlcond)) %>%
    pivot_longer(cols=!c(id,Dose), names_to = c(".value", "condition"), names_sep = "_") %>%
    na.omit() %>% #to remove NAs from conditions in which the protein was not detected
    mutate(Dose=as.numeric(Dose)) %>%
    mutate(gap=sum.conf.int/abs(sub.ratio)) %>%
    mutate(ci=predict(fitted,data.frame(x=.$gap))) %>%
    
    
    group_by(condition) %>%
    mutate(ci.adj=p.adjust(ci,method="BH")) %>%
    ungroup() %>%
    
    group_by(id,condition,R2,response) %>%
    summarise(sub.ratio=sum(sub.ratio,na.rm=TRUE), ci.mean=prod(ci.adj,na.rm=TRUE)^(1/n())) %>%
    ungroup() %>%
    
    mutate(weighted=weighted) %>%
    mutate(R2weight=ifelse(weighted,(2-R2)^R2w,1)) %>%
    mutate(ci.wt = ci.mean*R2weight) %>%
    
    group_by(id) %>%
    summarise(dAUC=sum(sub.ratio)/n(), 
              CI= prod(p=ci.wt,na.rm=TRUE), 
              R2mean=1-prod(1-R2,na.rm=TRUE),
              max.response=max(response,na.rm=TRUE)) %>% 
    rename(R2=R2mean) %>%
    na.omit() %>%
    ungroup()
  
  if(is.na(nMAD) & is.na(mindAUC)) {
    stop("Please give either nMAD or mindAUC, or both.")
  } else if(is.na(nMAD)) {
    limits <- c(mindAUC,-mindAUC)
  } else if(is.na(mindAUC)) {
    limits <- nMAD*mad(hit_data$dAUC, na.rm=TRUE) * c(-1,1) + median(hit_data$dAUC)
  } else {
    limits1 <- c(mindAUC,-mindAUC)
    limits2 <- nMAD*mad(hit_data$dAUC, na.rm=TRUE) * c(-1,1) + median(hit_data$dAUC)
    
    limits <- ifelse(limits1[1]>limits1[2],limits1,limits2)
    cat("Both nMAD and mindAUC given. Using the higher limit of the two.\n")
  }
  
  if(!is.na(minresponse)) {
    minresponse = minresponse + 1
  } else {
    minresponse=1
  }

  hit_data <- hit_data %>%
    mutate(hit=ifelse(CI<=0.05&(dAUC<=limits[1]|dAUC>=limits[2])&R2>=R2line&abs(max.response)>=minresponse,"hit","")) %>%
    mutate(Stabilization=ifelse(dAUC<0,"Destabilized","Stabilized"))
  
  labels <- interaction(unique(hit_data$Stabilization),unique(hit_data$hit), sep=" ")
  
  
  plot.settings.defaults <- list(labels=TRUE,
                                 label.text.size=2.5,
                                 label.force=1.3,
                                 xlims=c(-max(c(abs(hit_data$dAUC)),2), +max(c(2,abs(hit_data$dAUC)))),
                                 ylims=c(min(-log10(hit_data$CI)),max(-log10(hit_data$CI))),
                                 point.sizes=c(1,2),
                                 point.colors=c("gray","red","green"),
                                 axis.title.size=18,
                                 axis.text.size=16,
                                 legend.position="bottom",
                                 legend.text.size=8,
                                 legend.title=element_blank(),
                                 POI.color="blue",
                                 POI.size=3.5
  )
  
  for(psp in names(plot.settings.defaults)) {
    if(!psp %in% names(plot.settings)) {
      plot.settings[[psp]] <- plot.settings.defaults[[psp]]
    }
  }
  
  if(length(POI)==0) {
    addPOI=FALSE 
  } else {
    POIdata <- hit_data %>%
      filter(id %in% POI)
    
    if(nrow(POIdata)==0) {
      message("POI not present in the data.")
      addPOI=FALSE
    } else {
      addPOI=TRUE
    }
  }
  
  
  
  hit_plot <- hit_data %>% ggplot(aes(dAUC,-log(CI,10))) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray80") +
    geom_vline(xintercept=limits[1], linetype="dashed", color="gray80") + geom_vline(xintercept=limits[2], linetype="dashed", color="gray80") +
    geom_point(aes(color=interaction(Stabilization,hit, sep=" "),alpha=hit,size=interaction(Stabilization,hit, sep=" "))) +
    scale_color_manual(values=c(plot.settings$point.colors[1],plot.settings$point.colors[1],plot.settings$point.colors[2],plot.settings$point.colors[3]), name="Proteins", drop=FALSE)  +
    scale_alpha_manual(values=c(0.3,1), name="Stabilization") +
    scale_size_manual(values=c(plot.settings$point.sizes[1],plot.settings$point.sizes[1],plot.settings$point.sizes[2],plot.settings$point.sizes[2]), name="Proteins",drop=FALSE) +
    customPlot +
    scale_x_continuous(limits=plot.settings$xlims, name="dAUC", expand = c(0,0)) +
    scale_y_continuous(name="Confidence index", limits=plot.settings$ylims) +
    guides(alpha="none") + 
    theme(legend.position=plot.settings$legend.position,
          axis.title=element_text(size=plot.settings$axis.title.size),
          axis.text=element_text(size=plot.settings$axis.text.size),
          legend.text=element_text(size=plot.settings$legend.text.size),
          legend.title=plot.settings$legend.title)
  
  if(addPOI) {
    hit_plot <- hit_plot + 
      geom_point(data=POIdata, aes(dAUC,-log(CI,10)), color=plot.settings$POI.color,size=plot.settings$POI.size)
  }
  
  if(plot.settings$labels) {
    hit_plot <- hit_plot + geom_label_repel(data=hit_data%>%filter(hit=="hit"), aes(label=id), size=plot.settings$label.text.size, force=plot.settings$label.force,
                                            fill=alpha(c("white"),0), max.time=2, box.padding=0, label.size=NA)
  }
  

  hit_list=list(
    "Stabilized"=hit_data %>% filter(hit=="hit"&Stabilization=="Stabilized") %>% pull(id),
    "Destabilized"=hit_data %>% filter(hit=="hit"&Stabilization=="Destabilized") %>% pull(id)
  )
  
  output <- list("data"=hit_data, "plot"=hit_plot, "hitlist"=hit_list)
  return(output)
  
}