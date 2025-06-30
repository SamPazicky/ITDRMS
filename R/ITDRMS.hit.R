#' ITDRMS.hit
#'
#' Identifies hits from fitted mass spec data.
#' @param data Data frame: Scaled data with removed outliers and fitting statistics, ideally $data element from ITDRMS.fit output.
#' @param CIthreshold Numeric vector: Confidence index cutoff. Default is c(0.05,0.1) for hits and candidates, respectively.
#' @param minresponse Numeric vector: Minimal change in the solubility to be considered hit or candidate. Default is c(1,0.5) for hits and candidates, respectively.
#' @param R2line Numeric vector: R2  cutoff: Proteins with R-squared value below this value will not be considered as hits or candidates despite favorable response and confidence index.
#' Default is c(0.6,0.7) for hits and candidates, respectively.
#' @param POI Character vector: ID of the protein or proteins to be highlighted on the plot.
#' @param plot.settings List of graphical settings for plot. The defaults are:
#' \preformatted{
#' list(
#'   labels = TRUE,
#'   label.text.size = 2.5,
#'   label.force = 1.3,
#'   xlims = c(-max(1, abs(hit_data$total.response)) - 1, 
#'             max(1, abs(hit_data$total.response)) + 1),
#'   ylims = c(min(-log10(hit_data$CI)), max(-log10(hit_data$CI))),
#'   point.sizes = c(1, 2),
#'   point.colors = c("gray", "red3", "green3", "coral2", "darkolivegreen2"),
#'   axis.title.size = 18,
#'   axis.text.size = 16,
#'   legend.position = "bottom",
#'   legend.text.size = 8,
#'   legend.title = element_blank(),
#'   POI.color = "blue",
#'   POI.size = 3.5
#' )
#' }
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
#' data_fitted <- ITDRMS.hit(data_fitted)
#' @export

ITDRMS.hit <- function(
    data=NULL,
    CIthreshold=c(0.05,0.1),
    minresponse=c(1,0.5),
    R2line=c(0.6,0.7),
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
  
  data$maxresp <- apply(data[,grepl("fit_",names(data))],1,function(x) max(x,na.rm=TRUE))
  data$minresp <- apply(data[,grepl("fit_",names(data))],1,function(x) min(x,na.rm=TRUE))
  
  responses <- data %>%
    dplyr::select(id,condition,starts_with("fit_")) %>%
    pivot_longer(cols=starts_with("fit"),names_to="concentration",values_to="fraction_soluble") %>%
    pivot_wider(id_cols=c(id,concentration),names_from=condition,values_from=fraction_soluble) %>%
    mutate(across(!c(id,concentration), ~ .x - !!sym(controlcond), .names = "response.{col}")) %>%
    dplyr::select(id,concentration,starts_with("response")) %>%
    pivot_longer(cols=starts_with("response"), names_to="condition",values_to="response", names_prefix="response.") %>%
    filter(condition!=controlcond) %>%
    na.omit() %>%
    group_by(id,condition) %>%
    dplyr::summarise(maxresp=max(response,na.rm=TRUE), minresp=min(response,na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(response=ifelse(abs(maxresp)>abs(minresp),maxresp,minresp)) %>%
    dplyr::select(!ends_with("resp"))
    
  
  hit_data <- data %>%
    dplyr::select(!ends_with("resp")) %>%
    mutate(across(all_of(ratio_columns), ~ifelse(is.na(.x), NaN, .x))) %>% # to distinguish missing values (NA) and outliers (NaN)
    dplyr::select(id,condition,R2orig,matches(ratio_columns)) %>%
    mutate(R2orig=ifelse(is.na(R2orig), 0, R2orig)) %>%
    mutate(R2=ifelse(condition==controlcond,0,R2orig)) %>%
    rename_with(~ paste0("ratio_",.x),all_of(ratio_columns)) %>% # organizing the data to longer format
    pivot_longer(cols=!c(id,condition,R2), names_sep="_", names_to=c(".value","Dose")) %>%
    filter(Dose!=0) %>% # removing zero-dose values
    
    mutate(conf.int=ifelse(is.na(conf.int)&(!is.nan(ratio)), ratio/2,conf.int)) %>%
    
    filter(if_any(everything(), ~ !is.nan(.x))) %>% # removal of rows with removed outliers (NaN)
    pivot_wider(id_cols=c(id,Dose), names_from=condition, values_from=c(ratio,conf.int,fit,R2), names_sep="_") %>% # back to wide
    mutate(across(starts_with("fit_"), ~ .x - !!sym(paste0("ratio_",controlcond)), .names = "sub.{col}" )) %>% # subtract 37 ratio from higher-temp ratios
    mutate(across(starts_with("conf.int_"), ~ .x + !!sym(paste0("conf.int_",controlcond)), .names="sum.{col}")) %>% # sum 37 conf. interval and each higher-temp conf. interval
    dplyr::select(!contains(controlcond)) %>% # remove columns with control temperature condition
    pivot_longer(cols=!c(id,Dose), names_to = c(".value", "condition"), names_sep = "_") %>%# longer format
    na.omit() %>% #to remove NAs from conditions in which the protein was not detected
    mutate(Dose=as.numeric(Dose)) %>%
    mutate(gap=sum.conf.int/abs(sub.fit)) %>%
    mutate(ci=predict(fitted,data.frame(x=.$gap))) %>%
    
    
    group_by(condition) %>%
    mutate(ci.adj=p.adjust(ci,method="BH")) %>%
    ungroup() %>%
    
    group_by(id,condition,R2) %>%
    summarise(sub.fit=sum(sub.fit,na.rm=TRUE), ci.mean=prod(ci.adj,na.rm=TRUE)^(1/n())) %>%
    ungroup() %>%
    
    left_join(responses,by=c("id","condition")) %>%
    
    group_by(id) %>%
    summarise(dAUC=sum(sub.fit)/n(),
              CI=prod(p=ci.mean,na.rm=TRUE),
              R2prod=1-prod(1-R2,na.rm=TRUE),
              R2max=max(R2,na.rm=TRUE),
              sum.response=sum(response,na.rm=TRUE),
              max.response=ifelse(sum.response>0,max(response,na.rm=TRUE),min(response,na.rm=TRUE))) %>%
    na.omit() %>%
    ungroup() %>%
    mutate(total.response=sum.response+max.response) %>%
    mutate(Stabilization=ifelse(sum.response<0,"Destabilized","Stabilized"))
    
  
  if(!is.na(minresponse[2] & !is.na(CIthreshold[2] & !is.na(R2line[2])))) {
    hit_data <- hit_data %>%
      mutate(hit=ifelse(CI<=CIthreshold[2]&R2max>=R2line[2]&abs(total.response)>=minresponse[2],"candidate",""))
    secondline=TRUE
  } else {
    hit_data <- hit_data %>% mutate(hit="")
    secondline=FALSE
  }
  hit_data <- hit_data %>%
    mutate(hit=ifelse(CI<=CIthreshold[1]&R2max>=R2line[1]&abs(total.response)>=minresponse[1],"hit",hit))
  
  labels <- interaction(unique(hit_data$Stabilization),unique(hit_data$hit), sep=" ")
  
  
  plot.settings.defaults <- list(labels=TRUE,
                                 label.text.size=2.5,
                                 label.force=1.3,
                                 xlims=c(-max(1,abs(hit_data$total.response))-1, +max(1,abs(hit_data$total.response))+1),
                                 ylims=c(min(-log10(hit_data$CI)),max(-log10(hit_data$CI))),
                                 point.sizes=c(1,2),
                                 point.colors=c("gray","red3","green3","coral2","darkolivegreen2"),
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
  
  hit_plot <- hit_data %>% arrange(hit) %>%
    ggplot(aes(total.response, -log10(CI))) +
    geom_hline(yintercept = -log10(ifelse(secondline,CIthreshold[2],CIthreshold[1])), linetype = "dashed", color = "gray80") +
    geom_vline(xintercept = min(minresponse), linetype = "dashed", color = "gray80") +
    geom_vline(xintercept = (-1) * min(minresponse), linetype = "dashed", color = "gray80") +
    geom_point(aes(
      color = interaction(Stabilization, hit, sep = " "),
      alpha = interaction(Stabilization, hit, sep = " "),
      size = interaction(Stabilization, hit, sep = " ")
    )) +
    scale_color_manual(
      values = c(
        `Destabilized ` = plot.settings$point.colors[1],
        `Stabilized ` = plot.settings$point.colors[1],
        `Destabilized hit` = plot.settings$point.colors[2],
        `Stabilized hit` = plot.settings$point.colors[3],
        `Destabilized candidate` = plot.settings$point.colors[4],
        `Stabilized candidate` = plot.settings$point.colors[5]
      ),
      name = "Proteins", drop = TRUE
    ) +
    scale_alpha_manual(
      values = c(
        `Destabilized ` = 0.3,
        `Stabilized ` = 0.3,
        `Destabilized hit` = 1,
        `Stabilized hit` = 1,
        `Destabilized candidate` = 1,
        `Stabilized candidate` = 1
      ),
      name = "Proteins", drop = TRUE
    ) +
    scale_size_manual(
      values = c(
        `Destabilized ` = plot.settings$point.sizes[1],
        `Stabilized ` = plot.settings$point.sizes[1],
        `Destabilized hit` = plot.settings$point.sizes[2],
        `Stabilized hit` = plot.settings$point.sizes[2],
        `Destabilized candidate` = plot.settings$point.sizes[1],
        `Stabilized candidate` = plot.settings$point.sizes[1]
      ),
      name = "Proteins", drop = TRUE
    ) +
    customPlot +
    scale_x_continuous(limits=plot.settings$xlims, name="Total response", expand = c(0,0)) +
    scale_y_continuous(name="Confidence index", limits=plot.settings$ylims) +
    guides(alpha="none") + 
    theme(legend.position=plot.settings$legend.position,
          axis.title=element_text(size=plot.settings$axis.title.size),
          axis.text=element_text(size=plot.settings$axis.text.size),
          legend.text=element_text(size=plot.settings$legend.text.size),
          legend.title=plot.settings$legend.title)
  
  if(addPOI) {
    hit_plot <- hit_plot + 
      geom_point(data=POIdata, aes(total.response,-log(CI,10)), color=plot.settings$POI.color,size=plot.settings$POI.size)
  }
  
  if(plot.settings$labels) {
    hit_plot <- hit_plot + geom_label_repel(data=hit_data%>%filter(hit=="hit"), aes(label=id), size=plot.settings$label.text.size, force=plot.settings$label.force,
                                            fill=alpha(c("white"),0), max.time=2, box.padding=0, label.size=NA)
  }
  
  
  hit_list=list(
    "Stabilized_hits"=hit_data %>% filter(hit=="hit"&Stabilization=="Stabilized") %>% pull(id),
    "Destabilized_hits"=hit_data %>% filter(hit=="hit"&Stabilization=="Destabilized") %>% pull(id),
    "Stabilized_candidates"=hit_data %>% filter(hit=="candidate"&Stabilization=="Stabilized") %>% pull(id),
    "Destabilized_candidates"=hit_data %>% filter(hit=="candidate"&Stabilization=="Destabilized") %>% pull(id)
  )
  
  output <- list("data"=hit_data, "plot"=hit_plot, "hitlist"=hit_list)
  return(output)
  
}
