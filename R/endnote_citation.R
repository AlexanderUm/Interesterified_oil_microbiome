
#-------------------------------------------------------------------------------
# Save the R citation as EndNote
#-------------------------------------------------------------------------------
endnote_citation <- function(package_name, file_to_save){
  
  require(tidyverse)
  
  bibtex.cit <- capture.output(utils:::print.bibentry(citation(package_name), 
                                        style = "Bibtex"))
  
  trim_ent <- function(x) {grep(x, bibtex.cit, value = TRUE) %>% 
                           gsub(".* = \\{", "", .) %>% 
                           gsub("\\},.*", "", .)}
  
  authors <- trim_ent("author") %>%  
                str_split(., " and ") %>% 
                .[[1]] %>% 
                str_split(., " ") %>% 
                lapply(., function(x) {c("%A", 
                                         x[length(x)], 
                                         ",",
                                          x[1:(length(x)-1)])}) %>% 
                lapply(., function(x){paste(x, collapse = " ")}) %>% 
                unlist() %>% 
                gsub(" ,", ",", .)
  
  other <- c()
  
  
  other.filds <- grep(" = ", bibtex.cit, value = TRUE) %>% 
                      gsub(" = .*", "", .) %>% 
                      trimws() %>% 
                      .[. != "author"]
  
  for(i in other.filds) {
    
    prefix <- paste0("%", toupper(substr(i, 1, 1)), " ")
    
    if(i == "year") {prefix <- "%D "}
    if(i == "doi") {prefix <- "%R "}
  
    line.inst <- paste0(prefix, trim_ent(i)) 
    
    other <- c(other, line.inst)
    
  }
  
  other <- other[nchar(other) > 3]
  
  out.string <- c("%0 Journal Article", authors, other)
  
  write_lines(out.string, paste0(file_to_save, ".enw"))
  
}
