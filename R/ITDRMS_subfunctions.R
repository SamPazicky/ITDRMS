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


fit_sigmoid <- function(data, e_guess_lim=0.2) {
  
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

#' 
#' #' fit_MC
#' fit_MC <- function(data) {
#'   if(length(data$x)!=length(data$y)) {
#'     stop("Vectors x and y do not have the same length.")
#'   }
#'   # guess the initial value of inflection point for fitting. 
#'   # The guess is made when difference between subsequent scaled abundances is >0.25
#'   for (point in 1:(length(data$y)-1)) {
#'     y_dif <- abs(data$y[point+1] - data$y[point])
#'     if(y_dif>0.25) {
#'       e_guess <- data$x[point+1]
#'       e_vec <- point+1
#'       break
#'     } else {
#'       y_dif_saved <- y_dif
#'       e_guess <- median(data$x)
#'       e_vec <- 5
#'     }
#'   }
#'   #guess slope
#'   b_guess <- ((data$y[e_vec-1]-data$y[e_vec+1])*10)^3/7
#'   d_guess <- max(data$y)
#'   c_guess <- min(data$y)
#'   # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
#'   assign("my.fit.dat",
#'          try(minpack.lm::nlsLM(formula=y~(c+a*x+((d-c)/(1+exp(b*(log(x)-log(e)))))),
#'                                data=data,
#'                                start=list(a=0,b=b_guess,c=c_guess,d=d_guess, e=e_guess),
#'                                lower=c(-0.004,0,-0.35,0.85,37),
#'                                upper=c(0.004,150,0.26,1.1,73),
#'                                control=list(maxiter=100)),
#'              silent=TRUE)
#'   )
#'   
#'   return(my.fit.dat)
#' }
