#` @keywords internal
plot.ITDRcurve <- function(mergeddata, fakedata,print.stats,hits,ratio_columns,conditions,scale,pallete) {
  
  suppressPackageStartupMessages({
    library(tidyverse)
  })
  
  # extract variables
  ratio_data <- mergeddata$data
  fits <- mergeddata$fit
  p <- unique(ratio_data$id)
  

  # scale settings
  ticks <- 10^(-10:10)
  ticks <- ticks[c(min(which(ticks>min(as.numeric(ratio_columns[-1])))-1),which(ticks>min(as.numeric(ratio_columns[-1]))))]
  ticks <- ticks[c(which(ticks<max(as.numeric(ratio_columns))),max(which(ticks<max(as.numeric(ratio_columns))))+1)]
  ticks <- log10(ticks)

  
  cur_data <- ratio_data %>% 
    pivot_longer(cols=all_of(ratio_columns), names_to="Dose", values_to="Fraction_soluble",names_transform=list(Dose=as.double))
  label = str_trunc(unique(cur_data$label),35)
  
  cur_fakedata <- fakedata %>% setNames("Dose")
  for(T in unique(cur_data$condition)) {
    Tfit <- fits[[paste0(p,";",T)]]
    if(!is.null(Tfit)) {
      cur_fakedata <- bind_cols(cur_fakedata,
                                predict_clean(fits[[paste0(p,";",T)]], fakedata)%>%as.data.frame()%>%setNames(T)
      )
    }
  }
  
  if(length(names(cur_fakedata))==1) {
    return(NULL)
  }
  cur_fakedata <- cur_fakedata %>% 
    pivot_longer(cols=!Dose, names_to="condition", values_to="Fraction_soluble") %>%
    mutate(condition=factor(condition, levels=names(pallete)))
  
  cur_conditions <- sapply(cur_data$condition%>%unique()%>%gtools::mixedsort(), function(x) grep(x,conditions)) %>% unname()
  wishmin=0
  wishmax=2
  if(scale) {
    
    scmin <- min(cur_fakedata$Fraction_soluble,na.rm=TRUE)
    scmax <- max(cur_fakedata$Fraction_soluble,na.rm=TRUE)
    
    cur_fakedata <- cur_fakedata %>%
      mutate(Fraction_soluble=rescale(Fraction_soluble,from=c(scmin,scmax), to=c(0,1)))
    
    cur_data <- cur_data %>%
      mutate(Fraction_soluble=rescale(Fraction_soluble,from=c(scmin,scmax), to=c(0,1)))
  }
  
  miny <- min(wishmin,min(cur_data$Fraction_soluble, na.rm=TRUE))
  maxy <- max(wishmax,max(cur_data$Fraction_soluble, na.rm=TRUE))
  
  title <- paste0(p, "\n", label)
  
  if (!print.stats) {
    subtitle <- NULL
  } else {
    hit <- hits %>% dplyr::filter(id == p)
    hitR2 <- round(hit$R2, 2)
    hitp.val <- round(hit$p.adj, 2)
    subtitle <- paste0("R2=", hitR2, " p.adj=", hitp.val)
  }
  
  plot <- ggplot(data = cur_data) +
    geom_point(
               aes(x = log10(Dose), y = Fraction_soluble, color = condition),
               show.legend = TRUE) +
    geom_line(data = cur_fakedata,
              aes(x = log10(Dose), y = Fraction_soluble, color = condition),
              show.legend = TRUE) +
    scale_color_manual(values = pallete, drop = FALSE, limits = names(pallete)) +
    scale_x_continuous(breaks = ticks, name = "log(Dose)") +
    scale_y_continuous(limits = c(miny, maxy), name = "Fraction soluble") +
    ggtitle(title, subtitle = subtitle) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = grid::unit(c(5,5,5,5), "mm"),
      plot.title = element_text(size = 8, hjust = 0.5)
    )
  
  plot$plot_env <- NULL
  plot$mapping <- NULL
  plot$theme <- NULL
  plot$scales <- NULL
  plot$guides <- NULL
  plot$coordinates <- NULL
  plot$facet <- NULL
  
  # plotelems <- plot[c("theme","scales","guides","coordinates","facet")] #
  # saveRDS(plotelems,"inst/extdata/plotelems.RDS")
  return(plot)
  
}
