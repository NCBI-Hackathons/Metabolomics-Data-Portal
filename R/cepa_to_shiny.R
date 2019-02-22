
#' Cepa to shiny
#'
#' @param z_vec 
#' @param pathway.name 
#' @param pmap.path 
#' @param type 
#' @param thresh 
#' @param cen 
#' @param cen.name 
#' @param iter 
#'
#' @return CEPA
#' @export
#'
#' @examples
shiny.get.cepa <- function(z_vec, pathway.name = NULL, pmap.path="extdata", type = "ora",  
  thresh = 1.96, #required arguments
  cen = "betweenness", 
  cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
  iter = 1000){
  
  if(type == "ora"){
    cepa.ora.metab.all(names(z_vec[abs(z_vec) > thresh]), names(z_vec), pmap.path, pathway.name, cen, cen.name, iter)
  }else{
    cepa.univariate.metab.all(z_vec, pmap.path, pathway.name, cen, cen.name, iter)
  }
}
