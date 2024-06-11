#' ITDRMS.plot
#'
#' Identifies hits from fitted mass spec data.
#' @param protein Character string: Identifier of a protein that can be found in 'id' column of source.data.
#' @param source.data Data frame with measured melting curve, with one column for protein 'id' and other columns for scaled(zero to one) values
#' @param temperatures vector: for which temperatures should POS be calculated? Default is c(50,55,60)
#' @return Data frame with calculated plausibility values.
#' @examples 
#' plausibilities <- ITDRMS.plot("PF3D7_1017500", ITDRMS::source.data.lysate)
#' @export

ITDRMS.POS <- function(
    protein=NA,
    source.data=Pf_meltcurve_lysate,
    temperatures=c(50,55,60) # temperatures to be assessed
) {
  
  if(is.na(protein)) {
    stop("Input protein name.")
  }
  
  protein_reference <- source.data %>%
    as.data.frame() %>%
    filter(id==protein) %>%
    pivot_longer(cols=!id, names_to="temperature",values_to="value") %>%
    dplyr::select(!id) %>%
    mutate(temperature=as.numeric(str_remove(temperature,"X")))
  
  if(nrow(protein_reference)==0) {
    plausabilities <- data.frame("Temperature"=temperatures,
                                 "Stabilization plausability"=rep(NA,length(temperatures)),
                                 "Destabilization plausability"=rep(NA,length(temperatures)),
                                 "SD"=rep(NA,length(temperatures))
    )
  } else {
    
    plausabilities <- data.frame()
    data <- protein_reference %>% setNames(c("x","y"))%>%na.omit()
    fit <- fit_MC(data)
    
    if(class(fit)=="try-error") {
     
      plausabilities <- bind_rows(plausabilities,
                                  data.frame("Temperature"=temperatures,
                                             "Prediction"=rep(NA,length(temperatures))
                                  )
      )
    } else {
      
      prediction <- predict(fit, data.frame(x=temperatures))
      
      plausabilities <- bind_rows(plausabilities,
                                  data.frame("Temperature"=temperatures,
                                             "Prediction"=prediction
                                  )
      )
    }
      
    plausabilities <- plausabilities %>%
      group_by(Temperature) %>%
      dplyr::summarise(mean_pred=mean(Prediction,na.rm=TRUE)) %>%
      mutate(op_mean_pred=1-mean_pred) %>%
      dplyr::select(Temperature, op_mean_pred,mean_pred) %>%
      setNames(c("Temperature","Stabilization.plausability", "Destabilization.plausability"))
  }
  
  return(plausabilities)

}
