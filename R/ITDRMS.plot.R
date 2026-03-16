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
#' @param plots_per_page Integer: How many plots on one pdf page? Default is 20.
#' @param ncores Integer: How many cores to use for fitting. Default is 1.
#' @param ram Integer: Size allowed for parallel computing in GB.
#' 
#' @import tidyverse
#' @import magrittr
#' @importFrom gtools mixedsort
#' @import patchwork
#' @import furrr
#' @importFrom qpdf pdf_combine
#' 
#' @return A list of all plots.
#' @examples 
#' data_plotted <- ITDRMS.plot(data_fitted$data, data_fitted$fits, label.col="label",ncores=4)
#' @export

ITDRMS.plot <- function(
    data=NULL,
    fits=NULL,
    hits=NULL,
    calc.POS=NULL,
    POS.source=ITDRMS::Pf_meltcurve_lysate,
    scale=FALSE,
    print.stats=FALSE,
    color.scheme="rainbow",
    label.col=NA,
    fit.length=100,
    pdf.export=FALSE,
    pdf.folder=".",
    pdf.name="ITDR_curves",
    plots_per_page=20,
    ncores=1,
    ram=8
) 
{
  
  on.exit({
    dev.off()
    closeAllConnections()
    if(exists("default_warning_handler")) {
      options(warn = default_warning_handler)
    }
  })
  
  # forcing evaluation of parameters directly used in the future loop.
  invisible(force(pdf.folder))
  invisible(force(pdf.name))
  invisible(force(plots_per_page))

  customPlot <- list(
    theme_bw(base_size = 12),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=ggplot2::margin(5,5,5,5, "pt")
    )
  )
  
  data <- data %>%
    mutate(condition=factor(condition,levels=gtools::mixedsort(unique(condition))))
  
  if(color.scheme=="rainbow") {
    pallete <- rainbow(length(levels(data$condition)))
  } else if(color.scheme=="distinct") {
    pallete <- Pal25 <- c(
      "dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
      "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
      "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown"
    ) %>% .[1:length(levels(data$condition))]
  } else if (length(color.scheme)>=length(levels(data$condition))) {
    pallete <- color.scheme
    color.scheme <- "distinct"
  } else {
    stop("color.scheme must be 'rainbow', 'distinct' or a vector with color names or color codes. The vector must be at least as long as number of temperatures used in the assay")
  }
  names(pallete) <- levels(data$condition)
  
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
        selection="Stabilization.plausibility"
      } else {
        if(type>0) {
          selection="Stabilization.plausibility"
        } else {
          selection="Destabilization.plausibility"
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
  
  conditions <- levels(data$condition)
  
  if(is.na(pdf.name)) {
    pdf.name="ITDRcurves"
  }
  
  ### DESIGN SETTINGS ###

  topconc <- max(as.numeric(ratio_columns))
  divfactor <- mynthroot(topconc/min(as.numeric(ratio_columns)[-which(ratio_columns=="0")]),(fit.length-1))
  fakedata <- data.frame("x"=topconc/divfactor^(0:(fit.length-1)) )
  
  proteins <- unique(data$id)

  # Split fits by accession
  accessions <- sub(";.*$", "", names(fits))
  fits_split <- split(fits, accessions)
  datafits <- ratio_data %>%
    as.data.table() %>%
    split(by="id")
  all_accessions <- names(datafits)
  merged <- setNames(
    Map(function(acc) list(
      data = datafits[[acc]],    # NULL if missing
      fit  = fits_split[[acc]]   # NULL if missing
    ),
    all_accessions
    ),
    all_accessions
  )
  merged <- merged %>%
    # Remove top-level elements where $fit becomes empty
    keep(~ {
      if (!is.null(.x$fit)) {
        # Remove NULLs inside $fit
        .x$fit <- discard(.x$fit, is.null)
      }
      # Keep only if $fit exists and is not empty
      !is.null(.x$fit) && length(.x$fit) > 0
    })
  
  plots <- list()
  # plotelems <- readRDS("inst/extdata/plotelems.RDS")
  plotelems <- system.file("extdata", "plotelems.RDS", package = "ITDRMS")
  plotelems <- readRDS(plotelems)
  
  cat("Curve plotting in progress...\n")
  if(ncores==1) {
    pb <- txtProgressBar(min=0, max=length(merged), style=3, initial="")
    plots <- progress_lapply(merged, function(mergeddata) plot.ITDRcurve(mergeddata,fakedata,print.stats,hits,ratio_columns,conditions,scale,pallete),pb)
    close(pb)
    names(plots) <- names(merged)
  } else {
    
    plan(multisession, workers = ncores)
    handlers("txtprogressbar")
    
    with_progress({
      pr <- progressor(steps = length(merged))  
      
      plots <- future_map(
        merged,
        function(mergeddata) {
          pr()  
          plot.ITDRcurve(
            mergeddata, fakedata, print.stats, hits,
            ratio_columns, conditions, scale, pallete
          )
        },
        .options = furrr_options(
          packages = c("ggplot2"),
          globals = c("fakedata", "print.stats", "hits", "ratio_columns", "conditions", "scale", "pallete","pr")
        )
      )
    })
    plan(sequential)
  }
  
  # add back plot elements and set y limits
  plots <- lapply(plots, function(x) {
    for (nm in names(plotelems)) {
      x[[nm]] <- plotelems[[nm]]
    }
    x
  })
  
  for(i in names(plots)) {
    plots[[i]] <- ggplot2:::plot_clone(plots[[i]])
    dmin <- min(plots[[i]]$data$Fraction_soluble, na.rm = TRUE) - 0.5
    dmax <- max(plots[[i]]$data$Fraction_soluble, na.rm = TRUE) + 0.5
    plots[[i]]$scales$scales[[3]]$limits <- c(min(0,dmin),max(2,dmax))
  }

  output <- plots
  
  if(!pdf.export) {
    return(output)
  }
    
  cat("Saving pdf file. This may take several minutes, depending on the number of detected proteins...\n")

  # points and lines
  pointshapes<-c(16,15,17,18,8,4,3,9,10,11,12,13,14)
  linetypes <- c("solid","dashed","dotted")
  
  # prepare legend
  
  realcondNo <- length(unique(data$condition))
  realrepNo <- length(unique(data$Replicate))
  legenddata <- data.frame("Dose"=rep(ratio_columns,realcondNo*realrepNo),
                           "condition"=rep(conditions,realrepNo,each=length(ratio_columns)),
                           "Replicate"=rep(realrepNo,each=length(ratio_columns)*realcondNo)
  )
  
  legenddata$value <- sample(0:1000000,size=nrow(legenddata))/1000000
  legenddata$condition <- factor(legenddata$condition,levels=gtools::mixedsort(unique(legenddata$condition)))
  
  suppressMessages({
    legendplot <- ggplot(legenddata, aes(x=Dose,y=value)) +
      geom_point(aes(color=condition)) +
      geom_line(aes(color=condition)) +
      customPlot +
      scale_linetype_manual(name="Replicate",values=linetypes) +
      guides(color=guide_legend(title.position="top",title.hjust=0.5,nrow=1, ncol=length(conditions))) +
      scale_color_manual(name="Dose (μM)",values=pallete)
    
    invisible(cond_legend <- get_legend(legendplot))
  
    plots <- lapply(plots, function(x) x + theme(axis.text=element_text(size=6),
                                                 plot.title=element_text(hjust=0.5, size=8),
                                                 axis.title=element_blank(),
                                                 legend.position="none"
    ))
  })

  if(!is.null(calc.POS)) {
    for(p in proteins) {
      plots[[p]] <- plots[[p]] + 
        patchwork::inset_element(POS_plots[[p]],
                                 right=(0.05+0.1*length(calc.POS)), left=0.05,
                                 bottom=0.85,top=0.95, ignore_tag=TRUE)
    }
  }
  
  layout <- "
  AB
  #C"
  
  glob_lab <- "Fraction soluble"
  # y_lab <- 
  #   ggplot() + 
  #   annotate(geom = "text", x = 1, y = 1, label = glob_lab, angle = 90) +
  #   coord_cartesian(clip = "off")+
  #   theme_void()
  y_lab <- textGrob(
    "Fraction soluble",
    rot = 90,
    gp = gpar(fontsize = 12)
  )
  
  render_page <- function(page_plots, temp.folder, pdf.name, y_lab, cond_legend, 
                          layout, plots_per_page, pr = NULL) {
    
    page_index <- names(page_plots)[1]
    temp_file <- file.path(temp.folder, paste0(pdf.name, "_page_", page_index, ".pdf"))
    
    if(length(page_plots) < plots_per_page) {
      emptyplot <- ggplot() + theme(panel.border = element_blank(), 
                                    panel.background = element_rect(fill = "white"))
      epno <- plots_per_page - length(page_plots)
      eplist <- list()
      for (ep in 1:epno) eplist[[ep]] <- emptyplot
      page_plots <- c(page_plots, eplist)
    }
    
    pdf(temp_file, width = 10, height = 14.5)
    suppressWarnings(
      print(
        patchwork::wrap_elements(y_lab) +
          patchwork::wrap_plots(page_plots, ncol = 4, nrow = 5) +
          patchwork::wrap_elements(cond_legend) +
          patchwork::plot_layout(design = layout, heights = c(1, 0.05), widths = c(0.1, 1))
      )
    )
    dev.off()
    
    if (!is.null(pr)) tryCatch(pr(), warning = function(w) NULL)
    
    return(temp_file)  # full path to temp location
  }

  environment(render_page) <- asNamespace("ITDRMS")
  
  # Split plots into pages
  plot_chunks <- split(plots, (seq_along(plots) - 1) %/% 20)
  
  # prepare temp folder
  temp.folder <- tempfile()  # generates a unique path
  dir.create(temp.folder)    # creates it in the main session's tempdir
  
  if(ncores<=1) {
    pb <- txtProgressBar(min=0, max=length(plot_chunks), style=3, initial="")
    pdf_pages <- progress_lapply(plot_chunks, function(pc) render_page(pc, temp.folder, pdf.name, y_lab, cond_legend, layout,plots_per_page),pb)
    close(pb)
  } else {

    plan(multisession, workers = ncores)
    handlers("txtprogressbar")  # or "progress" for RStudio
    options(future.globals.maxSize = ram * 1024^3)
    handlers(
      handler_txtprogressbar(
        reporters = list(
          on_message = function(...) NULL  # silence
        )
      )
    )
    
    
    default_warning_handler <- getOption("warn")
    options(warn = -1)  # suppress all warnings
    
    
    with_progress({
      pr <- progressor(steps = length(plot_chunks))  
      
      pdf_pages <- future_map(
        plot_chunks, render_page, temp.folder, pdf.name = pdf.name, y_lab = y_lab, cond_legend = cond_legend, 
        layout = layout, plots_per_page = plots_per_page, pr = pr,
        .options = furrr_options(
          scheduling = Inf,
          packages = c("patchwork", "ggplot2", "grDevices","progressr"),
          globals = list(
            render_page=render_page, temp.folder=temp.folder, pdf.name=pdf.name,y_lab=y_lab,
            cond_legend=cond_legend,layout=layout,plots_per_page=plots_per_page
          )
        )
      )
    })
    plan(sequential)
    options(warn = default_warning_handler)
  }

  # Remove any failed pages
  pdf_pages <- pdf_pages[!sapply(pdf_pages, is.null)]

  # combine pdf pages
  output_pdf <- file.path(pdf.folder, paste0(pdf.name, ".pdf"))
  qpdf::pdf_combine(input = pdf_pages, output = output_pdf)

  invisible(suppressMessages(file.remove(unlist(pdf_pages))))
  cat("Final PDF saved to: ", output_pdf,"\n")
  
  return(output)
}

