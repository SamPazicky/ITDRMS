#' ITDRMS.load
#'
#' Loads ITDR-CETSA datafiles from Proteome Discover (Versions 2.0 and above) to R. Protein abundances should be exported as .txt files and their filenames must
#' indicate organism and temperature, for example 'CETSAITDR_cipargamin_37_Pf.txt'.
#' @param datafolder Directory with datafiles. The name of the files needs to the pattern (see below), temperature (e.g. _37_) and organism abbreviation (e.g. _Hs). Default is "datafiles"
#' @param pattern Character string that is contained in all filenames of files that are to be loaded. Default is ".*".
#' @param orgs Character vector listing files with data of what organisms there are in the datafolder. Default is c("Pf","Hs").
#' @param temperatures Integer vector with temperatures at which the assay was conducted. Default is c(37,53,59,65).
#' @param top.conc Integer. Top concentration of the drug dilution series.
#' @param dilutions Integer stating the dilution factor of the drug dilution series.
#' @param direction Character string. "decreasing" when the loaded data contains highest concentration first, "increasing" otherwise.
#' @param control.tube Character string. "last" if the vehicle control is last in the loaded data, otherwise "first".
#' @param tubes Integer. How many tubes (including vehicle control)? Default is 10 assuming TMT10 labeling.
#' @param concentrations Integer vector with the concentrations of drug used in the assay. Default is NA. If specified, it overwrites
#' dilution, tp.conc, direction, control.tube and tubes.
#' 
#' @import tidyverse
#' @import magrittr
#' @importFrom gtools mixedsort
#' 
#' @return List with two data.frames: $curves can be used in the next step of ITDRMS pipeline, $data can be saved as a csv.
#' @examples 
#' ITDRMS.load(temperatures=c(37,53,59,65))
#' @export


ITDRMS.load=function(datafolder="datafiles",
                   pattern=".*", 
                   orgs=c("Hs","Pf"),
                   temperatures=c(37,53,59,65),
                   top.conc=10,
                   dilution=4,
                   direction="decreasing",
                   control.tube="last",
                   tubes=10,
                   concentrations=NA) {
  

  if(any(is.na(concentrations))) {
    if(control.tube=="none") {
      concentrations=round(top.conc/dilution^(seq(0, tubes-1, by = 1)),digits=5)
    } else {
      concentrations=round(top.conc/dilution^(seq(0, tubes-2, by = 1)),digits=5)
    }
    if(direction=="increasing") {
      concentrations=rev(concentrations)
    }
    if(control.tube=="last") {
      concentrations=c(concentrations,0)
    } else if (control.tube=="first") {
      concentrations=c(0,concentrations)
    }
  }
  

  # adjust folder and datafolder
  # if(!grepl("\\/$",destfolder)) {
  #   destfolder=paste0(destfolder,"/")
  # }

  if(!grepl("\\/$",datafolder)) {
    datafolder=paste0(datafolder,"/")
  }
  
  raw_data <- data.frame()
  new_names <- paste0("Abundance",1:tubes)
  pfiles <- list.files(path=datafolder,pattern=pattern,full.names=TRUE)
  for(t in temperatures) {
    for(i in seq_along(orgs)) {
      readfile <- grep(orgs[i],pfiles, value=TRUE)
      readfile <- grep(paste0("_",as.character(t),"_"),readfile,value=TRUE)
      newfile <- read.csv(readfile,sep="\t")
      
      # find columns with abundances
      ab_cols <- newfile %>%
        dplyr::select(!contains("CV")) %>%
        dplyr::select(!contains("Ratio")) %>% 
        names()
      
      # find last column before all abundance and ratio columns
      last_notabcol <- names(newfile) %>%
        str_detect("bundanc") %>%
        which() %>% min()
      last_notabcol <- last_notabcol -1
      
      raw_data <- bind_rows(raw_data,
                            newfile %>%
                              mutate(condition=t) %>%
                              mutate(Organism=orgs[i]) %>%
                              dplyr::select(1:last_notabcol, condition, Organism, all_of(ab_cols)) %>%
                              rename_with(~new_names, all_of(ab_cols)) %>%
                              mutate(Ratio1=1) %>%
                              mutate(across(all_of(paste0("Abundance",2:tubes)), ~ .x/Abundance1, 
                                            .names="{str_replace(.col,'Abundance','Ratio')}"))
      )
    }
  }
  
  #rename raw data columns
  CoverageColumn <- which(grepl("Coverage",names(raw_data)))
  raw_data <- raw_data  %>%
    rename(Coverage....=all_of(CoverageColumn))
  
  nrow1 <- nrow(raw_data)
  
  raw_data <- raw_data %>%
    mutate(Accession=str_remove(Accession,"\\.[[:digit:]]-p[[:digit:]]+")) %>%
    group_by(condition) %>%
    distinct(Accession, .keep_all=TRUE) %>%
    ungroup()
  
  nrow2 <- nrow(raw_data)
  cat("Removed", nrow1-nrow2,"data points upon alternative transcript reduction.\n")
  
  # constructing meltcurve table
  meltcurve <- raw_data %>%
    # mutate("1"=1) %>%
    mutate(countNum=2)%>%
    dplyr::select(Accession,Description,condition,starts_with("Ratio"), starts_with("Abundance"), X..Unique.Peptides, X..PSMs,countNum, Organism) %>%
    setNames(c("id","description","condition",as.character(concentrations), paste0("Abundance_",concentrations),"sumUniPeps","sumPSMs","countNum","Organism")) %>%
    mutate(condition=as.character(condition)) %>%
    mutate(condition=ifelse(condition==as.character(min(temperatures)), paste0(condition,"C"),condition)) %>% # paste0(condition,".r1"))) %>%
    mutate(description=gsub(".*organism","organism", description)) %>%
    mutate(description=gsub("*gene_","", description))
  
  output=list(curves=meltcurve,data=raw_data)
  return(output)
}