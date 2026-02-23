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
    CIthreshold=c(1,0.75),
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
    dplyr::select(id,condition,starts_with("fit_")) %>% # extracts fit values
    pivot_longer(cols=starts_with("fit"),names_to="concentration",values_to="fraction_soluble") %>% # everything into long format
    pivot_wider(id_cols=c(id,concentration),names_from=condition,values_from=fraction_soluble) %>% # one column per each temperature
    mutate(across(!c(id,concentration), ~ .x - !!sym(controlcond), .names = "response.{col}")) %>% # subtract control condition, call it "response"
    dplyr::select(id,concentration,starts_with("response")) %>%  # keep only id, conc and response columns
    pivot_longer(cols=starts_with("response"), names_to="condition",values_to="response", names_prefix="response.") %>% # back to long format
    filter(condition!=controlcond) %>% # remove control condition
    na.omit() %>% # remove NA values
    group_by(id,condition) %>% # for each id and temperature, pick max and min response
    dplyr::summarise(maxresp=max(response,na.rm=TRUE), minresp=min(response,na.rm=TRUE),.groups="drop") %>%
    ungroup() %>%
    mutate(response=ifelse(abs(maxresp)>abs(minresp),maxresp,minresp)) %>% # keep only max or min response, depends on which is larger
    dplyr::select(!ends_with("resp"))
    
  
  hit_data <- data %>%
    dplyr::select(!ends_with("resp")) %>% # remove columns maxresp and minresp
    mutate(across(all_of(ratio_columns), ~ifelse(is.na(.x), NaN, .x))) %>% # to distinguish missing values (NA) and outliers (NaN)
    dplyr::select(id,condition,R2orig,matches(ratio_columns)) %>% # columns: id, temperatures, R2 and ratio columns, conf.ints and fits
    # filter(!(condition == controlcond & R2orig < 0.5)) %>% # remove low-quality control conditions
    mutate(R2orig=ifelse(is.na(R2orig), 0, R2orig)) %>% # adjust R2 to 0 if NA
    mutate(R2=ifelse(condition==controlcond,0,R2orig)) %>% # adjust R2 in control condition to 0
    rename_with(~ paste0("ratio_",.x),all_of(ratio_columns)) %>% # organizing the data to longer format
    pivot_longer(cols=!c(id,condition,R2), names_sep="_", names_to=c(".value","Dose")) %>%
    filter(Dose!=0) %>% # removing zero-dose values. They were only used for alignmenet
    # mutate(conf.int=ifelse(is.nan(conf.int), 0,conf.int)) %>% # I am not convinced that replacing with 0 is best. removing is better. Anyways it is mutliplied.
    # filter(!is.na(conf.int)) %>%
    filter(!is.na(conf.int)&!is.nan(conf.int)&!conf.int==0) %>%
    # mutate(conf.int=ifelse(is.na(conf.int)&(!is.nan(ratio)), ratio/2,conf.int)) %>%
     
    # filter(if_any(everything(), ~ !is.nan(.x))) %>% # removal of rows with removed outliers (NaN)
    pivot_wider(id_cols=c(id,Dose), names_from=condition, values_from=c(ratio,conf.int,fit,R2), names_sep="_") %>% # back to wide
    mutate(across(starts_with("fit_"), ~ .x - !!sym(paste0("ratio_",controlcond)), .names = "sub.{col}" )) %>% # subtract 37 ratio from higher-temp ratios
    mutate(across(starts_with("conf.int_"), ~ .x + !!sym(paste0("conf.int_",controlcond)), .names="sum.{col}")) %>% # sum 37 conf. interval and each higher-temp conf. interval
    dplyr::select(!contains(controlcond)) %>% # remove columns with control temperature condition
    pivot_longer(cols=!c(id,Dose), names_to = c(".value", "condition"), names_sep = "_") %>%# longer format
    na.omit() %>% #to remove NAs from conditions in which the protein was not detected
    mutate(Dose=as.numeric(Dose)) %>%
    mutate(pointgap=abs(sub.fit)-sum.conf.int/2) %>%
    mutate(adsign=sign(sub.fit)) %>%
    mutate(pgsign=sign(pointgap)) %>%
    # filter(pgsign==pgsign2) %>%
    # mutate(pointgap=ifelse(pgsign>0,1+pointgap,1-pointgap)) %>%
    left_join(responses,by=c("id","condition"),relationship = "many-to-one") %>%
    
    group_by(id, condition) %>%
    dplyr::mutate(
      R2   = mean(R2, na.rm = TRUE),
      dAUC = sum(sub.fit) / n(),
      response=mean(response,na.rm=TRUE),
      .groups = "drop"
    ) %>%
    
    group_by(id) %>%
    summarise(dAUC=sum(dAUC)/n(),
              CI= sum(pointgap),
              R2prod=1-prod(1-R2,na.rm=TRUE),
              R2max=max(R2,na.rm=TRUE),
              mean.response=mean(response,na.rm=TRUE),
              max.response=ifelse(mean.response>0,max(response,na.rm=TRUE),min(response,na.rm=TRUE))) %>%
    ungroup() %>%
    mutate(total.response=mean.response+max.response) %>%
    mutate(Stabilization=ifelse(mean.response<0,"Destabilized","Stabilized"))
    
  
  if(!is.na(minresponse[2] & !is.na(CIthreshold[2] & !is.na(R2line[2])))) {
    hit_data <- hit_data %>%
      mutate(hit=ifelse(CI>=CIthreshold[2]&R2max>=R2line[2]&abs(total.response)>=minresponse[2],"candidate",""))
    secondline=TRUE
  } else {
    hit_data <- hit_data %>% mutate(hit="")
    secondline=FALSE
  }
  hit_data <- hit_data %>%
    mutate(hit=ifelse(CI>=CIthreshold[1]&R2max>=R2line[1]&abs(total.response)>=minresponse[1],"hit",hit))
  
  labels <- interaction(unique(hit_data$Stabilization),unique(hit_data$hit), sep=" ")
  
  
  plot.settings.defaults <- list(labels=TRUE,
                                 label.text.size=2.5,
                                 label.force=1.3,
                                 xlims=c(-max(1,abs(hit_data$total.response))-1, +max(1,abs(hit_data$total.response))+1),
                                 ylims=c(min(hit_data$CI),max(hit_data$CI)),
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
    ggplot(aes(total.response, CI)) +
    geom_hline(yintercept = ifelse(secondline,CIthreshold[2],CIthreshold[1]), linetype = "dashed", color = "gray80") +
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
      geom_point(data=POIdata, aes(total.response,CI), color=plot.settings$POI.color,size=plot.settings$POI.size)
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
