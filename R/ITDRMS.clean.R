#' ITDRMS.clean
#'
#' Cleans the loaded mass spec data
#' @param data Data frame with the loaded mass spec data
#' @param keratin Logical: Should the keratins be removed?
#' @param serum Logical: Should the serum proteins be removed?
#' @param trypsin Logical: Should trypsin be removed?
#' 
#' @import tidyverse
#' @import magrittr
#' 
#' @return List of two elements. $data is the data frame with contaminants removed and $removed is the data frame listing the removed contaminants.
#' @examples 
#' data_cleaned <- ITDRMS.clean(data_raw)
#' @export


ITDRMS.clean <- function(
    data=NULL,
    keratins=TRUE,
    serum=TRUE,
    trypsin=TRUE
) {
  
  require(tidyverse)
  
  if(is.null(data)) {
    stop("Please include data argument")
  }
  
  ratio_columns <- suppressWarnings(try(names(data) %>% as.numeric() %>% .[!is.na(.)] %>% as.character(),silent=TRUE))
  
  # remove NA columns
  NAcolumns <- apply( (data%>%dplyr::select(all_of(ratio_columns))),1,function(x) any(is.na(x))) %>% which()
  removed <- data %>% slice(NAcolumns) %>% dplyr::select(id,condition) %>% mutate(removal="NAvalues")
  data <- data %>% slice(-NAcolumns)
  
  # remove keratins 
  if(keratins) {
    removed <- bind_rows(removed,
                         data %>% filter(grepl("Keratin",description)) %>% dplyr::select(id,condition) %>% mutate(removal="Keratin")
                         )
    data <- data %>% filter(!grepl("Keratin", description))
  }
  
  # remove serum albumin
  if(serum) {
    removed <- bind_rows(removed,
                         data %>% filter(grepl("Serum albumin",description)) %>% dplyr::select(id,condition) %>% mutate(removal="Serum_albumin")
    )
    data <- data %>% filter(!grepl("Serum albumin", description))
  }
  
  # remove trypsin
  if(trypsin) {
    removed <- bind_rows(removed,
                         data %>% filter(grepl("Trypsin",description)) %>% dplyr::select(id,condition) %>% mutate(removal="Trypsin")
    )
    data <- data %>% filter(!grepl("Trypsin", description))
  }
  
  output <- list("data"=data, "removed"=removed)
  return(output)
  
  cat("Contaminants removed.")
  
}