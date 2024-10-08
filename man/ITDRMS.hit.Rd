% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITDRMS.hit.R
\name{ITDRMS.hit}
\alias{ITDRMS.hit}
\title{ITDRMS.hit}
\usage{
ITDRMS.hit(
  data = NULL,
  weighted = TRUE,
  R2w = 10,
  nMAD = 3,
  R2line = 0.6,
  plot.settings = list()
)
}
\arguments{
\item{data}{Data frame: Scaled data with removed outliers and fitting statistics, ideally $data element from ITDRMS.fit output.}

\item{weighted}{Logical: If TRUE (default), the resulting adjusted p values will be weighted by R-squared measure of goodness of fit.}

\item{R2w}{Integer: How much should the R-squared weight in? Default is 10.}

\item{nMAD}{Integer: How much higher or lower does the area under dose-response curve (dAUC) be than mean absolute deviation (MAD)? Defualt is 3.}

\item{R2line}{Double: Soft R2 cut-off. Proteins with R-squared value below this value will not be considered as hits despite favourable dAUC and p-value.}

\item{plot.settings}{List of graphical settings for plot. Defaults are: 
list(labels=TRUE, label.text.size=2.5, label.force=1.3,
xlims=c(-max(c(abs(hit_data$dAUC)),2), +max(c(2,abs(hit_data$dAUC)))),
ylims=c(min(-log10(hit_data$CI)),max(-log10(hit_data$CI))),
 point.sizes=c(1,2), point.colors=c("gray","red","green"),
 axis.title.size=18, axis.text.size=16,
 legend.position="bottom", legend.text.size=8,legend.title=element_blank())}
}
\value{
A list with three elements. $data is a data frame with dAUC and p-values for all proteins, $plot is the volcano plot based on the hits and 
$hitlist is a list with two vectors: one for Stabilized and one for Destabilized hits.
}
\description{
Identifies hits from fitted mass spec data.
}
\examples{
data_fitted <- ITDRMS.hit(data_fitted, nMAD=0.2)
}
