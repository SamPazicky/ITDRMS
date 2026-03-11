#' range.scale

range.scale <- function(x) {
  z <- (x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
  return(z)
}



#' get_legend

get_legend<-function(myggplot) { # to extract legend from a single ggplot
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' fit_sigmoid


fit_sigmoid <- function(data, tryfixedslope=TRUE) {
  
  require(drc)
  
  if(length(data$x)!=length(data$y)) {
    stop("Vectors x and y do not have the same length.")
  }
  
  my.fit.dat <-
    tryCatch(
      drc::drm(formula = y ~ x, data=data, fct = LL.4(fixed=c(NA,NA,1,NA)), na.action=na.omit,
               # start=c(b=b_guess,c=c_guess, e=e_guess),
               lower=c(1,-30,data$x[2]),
               upper=c(20,30,data$x[length(data$x)-1])
      ),
  
      error=function(cond){
        return(NA)
      },
      finally={
      }
    )
  
  if(class(my.fit.dat)=="logical") {
    if(tryfixedslope) {
      my.fit.dat <-
        tryCatch(
          drc::drm(formula = y ~ x, data=data, fct = LL.4(fixed=c(NA,NA,1,1)), na.action=na.omit,
                   # start=c(b=b_guess,c=c_guess, e=e_guess),
                   lower=c(1,-30),
                   upper=c(20,30)
          ),
          
          error=function(cond){
            return(NA)
          },
          finally={  
          }
          
        )
    }
    if(class(my.fit.dat)=="try-error") {
      my.fit.dat=NA
    }
  }
  return(my.fit.dat)
}

#' fit_any_sigmoid

fit_any_sigmoid <- function(data, e_guess_lim=0.2) {
  
  require(drc)
  
  if(length(data$x)!=length(data$y)) {
    stop("Vectors x and y do not have the same length.")
  }
  
  my.fit.dat <-
    tryCatch(
      drc::drm(formula = y ~ x, data=data, fct = LL.4(fixed=c(NA,NA,NA,NA)), na.action=na.omit,
               # start=c(b=b_guess,c=c_guess, e=e_guess),
               lower=c(1,-30,-40,data$x[2]),
               upper=c(20,30,20,data$x[length(data$x)-1])
      ),
      
      error=function(cond){
        return(NA)
      },
      finally={  
      }
      
    )
  
  if(class(my.fit.dat)=="try-error") {
    my.fit.dat=NA
  }
  return(my.fit.dat)
}

#' mynthroot


mynthroot = function(x,n) {
  (abs(x)^(1/n))*sign(x)
}

#' rescale

rescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE), ...) {
  
  (x - from[1]) / diff(from) * diff(to) + to[1]
}


#' fit_MC
fit_MC <- function(data) {
  if(length(data$x)!=length(data$y)) {
    stop("Vectors x and y do not have the same length.")
  }
  # guess the initial value of inflection point for fitting.
  # The guess is made when difference between subsequent scaled abundances is >0.25
  for (point in 1:(length(data$y)-1)) {
    y_dif <- abs(data$y[point+1] - data$y[point])
    if(y_dif>0.25) {
      e_guess <- data$x[point+1]
      e_vec <- point+1
      break
    } else {
      y_dif_saved <- y_dif
      e_guess <- median(data$x)
      e_vec <- 5
    }
  }
  #guess slope
  b_guess <- ((data$y[e_vec-1]-data$y[e_vec+1])*10)^3/7
  d_guess <- max(data$y)
  c_guess <- min(data$y)
  # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
  assign("my.fit.dat",
         try(minpack.lm::nlsLM(formula=y~(c+a*x+((d-c)/(1+exp(b*(log(x)-log(e)))))),
                               data=data,
                               start=list(a=0,b=b_guess,c=c_guess,d=d_guess, e=e_guess),
                               lower=c(-0.004,0,-0.35,0.85,37),
                               upper=c(0.004,150,0.26,1.1,73),
                               control=list(maxiter=100)),
             silent=TRUE)
  )

  return(my.fit.dat)
}

#` progress_lapply
progress_lapply <- function(X, FUN, pb, silent=TRUE, ...) {
  on.exit({
    sink(NULL, type = "message")
    if(exists("zz")) {
      close(zz)
    }
  })
  result <- vector("list", length(X))
  
  if(silent) {
    zz <- file(tempfile(pattern = "fit_", fileext = ".txt"), open = "wt")
    sink(zz, type="message")
  }
  
  for (i in seq_along(X)) {
    result[[i]] <- FUN(X[[i]], ...)  # Apply function
    setTxtProgressBar(pb, i)  # Update progress
  }
  close(pb)  # Close progress bar when done
  return(result)
}


#` remove_outliers
remove_outliers <- function(
    x=c(),
    y=c(),
    noise=4,
    z.cutoff=0.4,
    max.out=8,
    if.max="remove" # 'remove' to remove all, 'outstanding' to remove outliers standing out, 'none' to keep
) {
  
  if(length(x)!=length(y)) {
    stop("Length x and y must be equal.")
  }
  if(max.out>length(y)-2) {
    max.out<-length(y)-2
    cat("max.out too high for the data size. Adjusted to ", max.out,".\n")
  }
  
  
  removing <- TRUE
  j=1
  outliers <- out.positions <- zs <- c()
  out.positions[1] <- 100 # arbitrarily high number to not remove any outlier in the first round of while loop.
  while(removing & j<=(max.out)) {
    p2pdiffs <-  (c(0,diff(y[-out.positions]),0) %>% diff())*c(2,rep(1,length(y[-out.positions])-2),2) %>% scale() %>% .[,1]*c(1,rep(1,length(y[-out.positions])-2),0.5)
    # if median is too high it means it is way too noisy to perform any removal:
    if(median(abs(p2pdiffs))/max(abs(p2pdiffs))>=noise) {
      removing <- FALSE
      next
    }
    p2ptable <- data.frame(x[-out.positions],p2pdiffs) %>% setNames(c("x","zdiff"))
    outlier <- p2ptable %>% filter(abs(zdiff)>=z.cutoff) %>% slice_max(abs(zdiff)) %>% pull(x) %>% as.character()
    if(length(outlier)==0) {
      removing <- FALSE
      next
    }
    outliers[j] <- outlier
    
    zs[j] <- p2ptable %>% filter(abs(zdiff)>=z.cutoff) %>% slice_max(abs(zdiff)) %>% pull(zdiff) %>% abs()
    out.positions[j] <- which(x==as.numeric(outliers[j]))
    j=j+1
  }
  
  
  if(length(outliers)==0) {
    outliers <- NULL
    # in case it just peeled of outliers from top to bottom conc
  } else if(length(outliers)>1 & outliers[1]==x[length(x)] & sum(diff(out.positions))==(-1)*length(out.positions)+1) {
    outliers <- NULL
  } else if(j>max.out) {
    zs <- sort(zs,TRUE)
    
    
    if(if.max=="remove") {
      outliers <- x
    } else if(if.max=="outstanding") {
      outliers <- x[out.positions[which(scale(zs,0)[,1]>=0.75)]]
    } else {
      outliers <- NULL
    }
  }
  
  return(outliers)
  
}

#` predict_clean
predict_clean <- function(clean_fit, x, dil.factor = exp(1)) {
  
  # =========================
  # LL.4 MODEL
  # =========================
  if (clean_fit$model_type == "LL4") {
    
    # Default all parameters to 1
    pars <- c(b = 1, c = 1, d = 1, e = 1)
    
    coefs <- clean_fit$coef
    names(coefs) <- sub(":.*", "", names(coefs))
    
    pars[names(coefs)] <- coefs
    
    return(
      pars["c"] + (pars["d"] - pars["c"]) /
        (1 + exp(pars["b"] * (log(x) - log(pars["e"]))))
    )
  }
  
  # =========================
  # LINEAR MODEL
  # =========================
  if (clean_fit$model_type == "lm") {
    
    coefs <- clean_fit$coef
    
    # Case 1: intercept-only model (y ~ 1)
    if (length(coefs) == 1) {
      return(rep(coefs[1], length(x)))
    }
    
    # Case 2: y ~ log(x, dil.factor)
    intercept <- coefs[1]
    slope <- coefs[2]
    
    return(intercept + slope * log(x, base = dil.factor))
  }
}

#` clean.fit
clean.fit <- function(fit)
{
  if (inherits(fit, "drc")) {
    
    cfit <- list(
      model_type = "LL4",
      coef = coef(fit)
    )
    
    # ---- linear fits ----
  } else if (inherits(fit, "lm")) {
    
    cfit <- list(
      model_type = "lm",
      coef = coef(fit)
    )
    
  } else {
    cfit <- NULL
  }
  return(cfit)
}

#` render_page
render_page <- function(page_plots, pdf.folder, pdf.name, y_lab, cond_legend, layout, plots_per_page) {
  # require(patchwork)
  # require(ggplot2)
  
  page_index <- names(page_plots)[1]
  temp_file <- file.path(pdf.folder, paste0(pdf.name, "_page_", page_index, ".pdf"))
  if(length(page_plots)<plots_per_page) {
    emptyplot <- ggplot() +  theme(panel.border=element_blank(), panel.background = element_rect(fill="white"))
    epno <- plots_per_page-length(page_plots)
    eplist=list()
    for (ep in 1:epno) {
      eplist[[ep]]=emptyplot
    }
    page_plots <- c(page_plots,eplist)
  }
  pdf(temp_file, width = 10, height = 14.5)
  print(
    wrap_elements(y_lab) +
      wrap_plots(page_plots, ncol = 4, nrow = 5) +
      wrap_elements(cond_legend)  +
      patchwork::plot_layout(
        design = layout,
        heights = c(1, 0.05),
        widths = c(0.1,1)
      )
  )
  dev.off()
  return(temp_file)
}
