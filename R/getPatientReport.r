#' Get Patient Report
#' @param PatientID - The patient identifier string associated with the patient's profile. 
#' @param all_data - The entire data matrix loaded based on the diagnosis selected in the dropdown menu input$diagClass.
#' @return patientReport - a data table with the metabolites and z-scores associated with the selected patient ID.
#' @export getPatientReport
#'
#' @examples
#' data(Miller2015_Heparin)
#' # Input is supplied by R shiny app, but you can hard code parameters as a list object, too, to test functionality.
#' input = list()
#' input$ptIDs = colnames(Miller2015_Heparin)[4]
#' input$diagClass = "paa"
#' rpt = getPatientReport(input, Miller2015_Heparin)
#' head(rpt$patientReport)
getPatientReport = function(input, all_data) {
  print(input$diagClass)
  print(input$ptIDs)
  
  # MetaboliteName Zscore
  all_data = data.matrix(all_data)
  tmp.zscore = rownames(all_data)
  if (length(input$ptIDs)>1) {
    zscore.data = apply(all_data[ , which(colnames(all_data) %in% input$ptIDs)], 1, function(i) mean(na.omit(i)))
  } else {
    zscore.data = all_data[ , which(colnames(all_data)==input$ptIDs)]
  }
  names(zscore.data) = tmp.zscore
  print(head(zscore.data))
  
  ind = which(is.na(zscore.data))
  if (length(ind)>0) {
    zscore.data = zscore.data[-ind]
  }
  
  data = data.frame(Metabolite=character(), Zscore=numeric(), stringsAsFactors = FALSE)
  for (row in 1:length(zscore.data)) {
    data[row, "Metabolite"] = names(zscore.data)[row]
    data[row, "Zscore"] = round(zscore.data[names(zscore.data)[row]], 2)
  }
  
  # Remove mets that were NA in zscore 
  ind0 = which(is.na(data[,"Zscore"]))
  if (length(ind0)>0) {
    data = data[-ind0,]
  }
  print(dim(data))
  
  # Order by abs(Zscore)
  class(data[,"Zscore"]) = "numeric"
  data = data[order(abs(data[,"Zscore"]), decreasing = TRUE), ]
  names(data) = c("Metabolite", "Z-score")
  
  return(list(patientReport=data))
}