#' ITDRMS.plot
#'
#' Identifies hits from fitted mass spec data.
#' @param data Data frame: Scaled data with removed outliers and fitting statistics, ideally $data element from ITDRMS.fit output.
#' @param fits List containing all fit objects, ideally $fits element from ITDRMS.fit output.
#' @param hits Data frame with hits coming from ITDRMS.hit. Only needed when print.stats=TRUE
#' @param calc.POS vector: for which temperatures should POS be calculated? Default is NULL, in which case POS will not be calculated. Use with care, ideally only with a small subselection of proteins, such as identified hits.
#' @param POS.source Data frame with melting curve data to calculate POS from. 
#' @param scale Logical: Should the curves for each protein be scaled from 0 to 1? Default is FALSE. If TRUE, the amplitude of response will not
#' be comparable across proteins!
#' @param print.stats Logical: Should the R2 and p-value be printed in the plots? Default is FALSE.
#' @param color.scheme Character string: "rainbow" will plot curves from red to violet from lowest to highest temperature, "distinct" will use distinct colors
#' @param label.col Character string: Name of the column to use for plot subtitles.
#' @param fit.length Integer: How many points should be used for fitting curves. Default is 100 which is sufficient for plotting.
#' @param pdf.export Logical: If TRUE (default), a pdf with all plots will be exported.
#' @param pdf.folder Character string: Name of the directory for pdf export. Default is the working directory.
#' @param pdf.name Character string: Name of tthe exported pdf file. Default is 'ITDR_curves'
#'
#' @import tidyverse
#' @import magrittr
#' @importFrom gtools mixedsort
#' @import patchwork
#' 
#' @return A list of all plots.
#' @examples 
#' data_plotted <- ITDRMS.plot(data_fitted$data, data_fitted$fits)
#' @export


ITDRMS.plot <- function(
    data=NULL,
    fits=NULL,
    hits=NULL,
    calc.POS=NULL,
    POS.source=Pf_meltcurve_lysate,
    scale=FALSE,
    print.stats=FALSE,
    color.scheme="rainbow",
    label.col=NA,
    fit.length=100,
    pdf.export=FALSE,
    pdf.folder=".",
    pdf.name="ITDR_curves"
    ) 
{
  
  customPlot <- list(
    theme_bw(base_size = 12),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=ggplot2::margin(5,5,5,5, "pt")
    )
  )
  
  if(print.stats&is.null(hits)) {
    stop("If you want to print stats, include data.frame in argument 'hits'.")
  }
  
  if(is.null(data)) {
    stop("Please include data")
  } else {
    data=as.data.frame(data)
  }
  
  if(!is.null(calc.POS)) {
    
    cat("Extracting POS curves...\n")
    
    all_hits <- data %>%
      as.data.frame() %>%
      dplyr::select(id,Slope)
    
    char_table <- data.frame()
    for(i in 1:nrow(data)) {
      plaus_table <- ITDRMS.POS(protein=data[i,"id"], source.data=POS.source, temperatures=calc.POS) %>%
        mutate(across(everything(), as.numeric))

      type=all_hits[i,"Slope"]
      if(is.na(type)) {
        seletion="Stabilization.plausibility"
      } else {
        if(type>0) {
          selection="Stabilization.plausability"
        } else {
          selection="Destabilization.plausability"
        }
      }
     
      
      char_table <- bind_rows(char_table,
                              plaus_table %>% dplyr::select(Temperature, all_of(selection)) %>%
                                column_to_rownames("Temperature") %>%
                                t() %>% as.data.frame() %>%
                                remove_rownames()
      )
    }
    

    POS_plots <- list()
    if(length(char_table)>0) {
      all_hits <- bind_cols(all_hits,char_table)
      for (i in all_hits$id) {
        POS_plots[[i]] <- all_hits %>%
          filter(id==i) %>%
          setNames(make.names(names(.))) %>%
          pivot_longer(cols=starts_with("X"), names_to="Temperature",values_to="Value") %>%
          mutate(Temperature=str_remove(Temperature,"X")) %>%
          ggplot() +
          geom_tile(aes(x=Temperature,y=id,fill=Value)) +
          geom_text(aes(x=Temperature,y=id,label=as.numeric(Temperature)), size=1.5) +
          scale_fill_gradient2(
            low = "red2",
            mid = "white",
            high = "green2",
            midpoint = 0.5,
            limits=c(0,1)
          ) +
          customPlot +
          scale_x_discrete(expand=c(0,0)) +
          scale_y_discrete(expand=c(0,0)) +
          # coord_fixed(5) +
          theme(legend.position="none",
                axis.text=element_blank(),
                axis.title=element_blank(),
                plot.margin=margin(0,0,0,0,"cm"),
                axis.ticks=element_blank(),
                panel.background=element_rect(fill="transparent"),
                plot.background=element_rect(fill="transparent")
          ) +
          coord_equal()
      }
    }
  }
  
  if(is.null(data$Replicate)) {
    data$Replicate=1
  }
  ratio_columns <- names(data)[4:13]
  if(is.na(label.col)) {
    ratio_data <- data %>%
      mutate(label="")
  } else {
    ratio_data <- data %>% 
      rename(label=!!sym(label.col)) %>%
      dplyr::select(id,condition,label,all_of(ratio_columns))
  }
 
  conditions <- data$condition %>% unique() %>% gtools::mixedsort()
  
  if(is.na(pdf.name)) {
    pdf.name="ITDRcurves"
  }
  
  ### DESIGN SETTINGS ###
  
  # color settings
  
  if(color.scheme=="rainbow") {
    pallete <- rainbow(length(unique(data$condition)))
  } else if(color.scheme=="distinct") {
    pallete <- Pal25 <- c(
      "dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
      "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
      "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown"
    )
  } else if (length(color.scheme)>=length(unique(data$condition))) {
    pallete <- color.scheme
    color.scheme <- "distinct"
  } else {
    stop("color.scheme must be 'rainbow', 'distinct' or a vector with color names or color codes. The vector must be at least as long as number of temperatures used in the assay")
  }
  
  # scale settings
  ticks <- 10^(-10:10)
  ticks <- ticks[c(min(which(ticks>min(as.numeric(ratio_columns[-1])))-1),which(ticks>min(as.numeric(ratio_columns[-1]))))]
  ticks <- ticks[c(which(ticks<max(as.numeric(ratio_columns))),max(which(ticks<max(as.numeric(ratio_columns))))+1)]
  ticks <- log10(ticks)
  
  # points and lines
  pointshapes<-c(16,15,17,18,8,4,3,9,10,11,12,13,14)
  linetypes <- c("solid","dashed","dotted")
  

  
  topconc <- max(as.numeric(ratio_columns))
  divfactor <- mynthroot(topconc/min(as.numeric(ratio_columns)[-which(ratio_columns=="0")]),(fit.length-1))
  fakedata <- data.frame("x"=topconc/divfactor^(0:(fit.length-1)) )
  
  proteins <- unique(data$id)
  plots <- list()
  
  
  message("Curve plotting in progress...")
  pb <- txtProgressBar(min=0, max=length(proteins), style=3, initial="")
  for(p in proteins) {
    
    cur_data <- ratio_data %>% 
      filter(id==p) %>%
      mutate(condition=factor(condition,levels=gtools::mixedsort(unique(ratio_data$condition)))) %>%
      pivot_longer(cols=all_of(ratio_columns), names_to="Dose", values_to="Fraction_soluble",names_transform=list(Dose=as.double))
    label = str_trunc(unique(cur_data$label),35)
    
    cur_fakedata <- fakedata %>% setNames("Dose")
    for(T in unique(cur_data$condition)) {
      Tfit <- fits[[paste0(p,";",T)]]
      if(class(Tfit) %in% c("lm","nls","drc")) {
        cur_fakedata <- bind_cols(cur_fakedata,
                                  predict(fits[[paste0(p,";",T)]], fakedata)%>%as.data.frame()%>%setNames(T)
        )
      }
    }
    
    if(length(names(cur_fakedata))==1) {
      next
    }
    cur_fakedata <- cur_fakedata %>% 
      pivot_longer(cols=!Dose, names_to="condition", values_to="Fraction_soluble")
    
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
    
    plots[[p]] <- ggplot() +
      geom_point(data=cur_data, aes(x=log10(Dose),y=Fraction_soluble,color=condition),show.legend=TRUE) +
      geom_line(data=cur_fakedata, aes(x=log10(Dose), y=Fraction_soluble,color=condition),show.legend=TRUE) +
      scale_color_manual(values=pallete, drop=FALSE) +
      scale_x_continuous(breaks=ticks,name="log(Dose)") +
      scale_y_continuous(limits=c(miny,maxy), name="Fraction soluble") +
      theme_bw(base_size = 12) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # legend.position=legend.position,
            plot.margin=grid::unit(c(5,5,5,5), "mm"),
            plot.title=element_text(size=8, hjust=0.5)
      )
    
    if(!print.stats) {
      plots[[p]] <- plots[[p]] +
        ggtitle(paste0(p,"\n",label))
    } else {
      hitR2=round(hits%>%filter(id==p)%>%pull(R2),2)
      hitp.val=round(hits%>%filter(id==p)%>%pull(p.adj),2)
      subtitle=paste0("R2=",hitR2," p.adj=",hitp.val)
      plots[[p]] <- plots[[p]] +
        ggtitle(paste0(p,"\n",label),
                subtitle=subtitle)
    }
    
    setTxtProgressBar(pb, which(proteins==p))
  }
  close(pb)
  
  output <- plots
  
  if(pdf.export) {
    
    # prepare legend
    
    realcondNo <- length(unique(data$condition))
    realrepNo <- length(unique(data$Replicate))
    legenddata <- data.frame("Dose"=rep(ratio_columns,realcondNo*realrepNo),
                           "condition"=rep(conditions,realrepNo,each=length(ratio_columns)),
                           "Replicate"=rep(realrepNo,each=length(ratio_columns)*realcondNo)
    )
    
    legenddata$value <- sample(0:1000000,size=nrow(legenddata))/1000000
    legenddata$condition <- factor(legenddata$condition,levels=gtools::mixedsort(unique(legenddata$condition)))
    
    legendplot <- ggplot(legenddata, aes(x=Dose,y=value)) +
      geom_point(aes(color=condition)) +
      geom_line(aes(color=condition)) +
      customPlot +
      scale_linetype_manual(name="Replicate",values=linetypes) +
      guides(color=guide_legend(title.position="top",title.hjust=0.5,nrow=1, ncol=length(conditions))) +
      scale_color_manual(name="Dose",values=pallete)
    
    cond_legend <- get_legend(legendplot)
    
    
    legendplot <- ggplot(legenddata, aes(x=Dose,y=value)) +
      geom_point(aes(shape=as.character(Replicate))) +
      geom_line(aes(linetype=as.character(Replicate))) +
      customPlot +
      scale_shape_manual(name="Replicate",values=pointshapes) +
      scale_linetype_manual(name="Replicate",values=linetypes) +
      guides(shape=guide_legend(nrow=realrepNo, ncol=1)) +
      guides(linetype=guide_legend(nrow=realrepNo, ncol=1)) +
      scale_color_manual(name="Dose",values=pallete)
  
    rep_legend <- get_legend(legendplot)
    rep_legend[["heights"]][1] <- unit(20.0000,"cm")
  }
  
  plots <- lapply(plots, function(x) x + theme(axis.text=element_text(size=6),
                                               plot.title=element_text(hjust=0.5, size=8)
                                               )
  )
  
  cat("Saving pdf file. This may take several minutes, depending on the number of detected proteins...\n")

  plots <- lapply(plots, function(x) x + 
                    theme(axis.title=element_blank(),
                          legend.title=element_blank()) + 
                    guides(colour = guide_legend(nrow = 1)) 
                  )
  
  if(!is.null(calc.POS)) {
    for(p in proteins) {
      plots[[p]] <- plots[[p]] + 
        patchwork::inset_element(POS_plots[[p]],
                                 right=(0.05+0.1*length(calc.POS)), left=0.05,
                                 bottom=0.85,top=0.95, ignore_tag=TRUE)
    }
  }
  
  output <- plots
  

  layout <- "
  ABBBBBBBBB
  #CCCCCCCCC
  #DDDDDDDDD"
  glob_lab <- "Fraction soluble"
  
  emptyplot <- ggplot() + customPlot + theme(panel.border=element_blank())
  eplist=list()
  for (ep in 1:20) {
    eplist[[ep]]=emptyplot
  }
  y_lab <- 
    ggplot() + 
    annotate(geom = "text", x = 1, y = 1, label = glob_lab, angle = 90) +
    coord_cartesian(clip = "off")+
    theme_void()
  
  glob_lab <- expression("log(Dose, " ~ mu~"M)")
  x_lab <- 
    ggplot() + 
    annotate(geom = "text", x = 1, y = 1, label = glob_lab, angle = 0) +
    coord_cartesian(clip = "off")+
    theme_void()
 
  pdf(paste0(pdf.folder,"/",pdf.name,".pdf"), width = 10 , height = 14.5)
  for (page in 1:ceiling(length(plots)/20)) {
    if(page==ceiling(length(plots)/20)) {
      print(
        patchwork::wrap_elements(y_lab) +
          patchwork::wrap_plots(c(plots[((20*page)-19):(length(plots))],eplist[0:(20-length(((20*page)-19):(length(plots))))]), ncol = 4, nrow=5 ) + 
          patchwork::wrap_elements(x_lab) +
          patchwork::guide_area() + patchwork::plot_layout(guides = 'collect', design=layout, heights=c(1,0.1,0.1))
  
      )
    } else {
      print(
        patchwork::wrap_elements(y_lab) + 
          patchwork::wrap_plots(plots[((20*page)-19):(20*page)], ncol = 4, nrow=5) + 
          patchwork::wrap_elements(x_lab) +
          patchwork::guide_area() + patchwork::plot_layout(guides = 'collect', design=layout, heights=c(1,0.1,0.1))
        
      )
    }
  }
  dev.off()
    
  return(output)
}