#' @title Converting atomic proportions into atomic proportions with different normalization methods and/or element choices
#' @description \code{Cation2Cation} converts atomic proportions into atomic proportions with different normalization methods and/or element choices.
#'
#' @param table (dataframe) Atomic proportion data
#' @param ElementList (dataframe) Information of considered elements, cation/oxygen numbers per unit formula.
#' Use MFRC::ElementList_all or its subset (e.g., MFRC::ElementList_MnNCKFMASTCr, MFRC::ElementList_NCFMAS).
#'
#' @param NormMode (chr) Set "Oxygen" or "Cation" (default = "Oxygen")
#' @param NormValue (num) Normalization value (default = 24)
#' @param padding (logical) TRUE for padding concentrations for element absent in the input data (default = FALSE)
#'
#'
#' @export
#' @examples
#' Mass2Cation(MFRC::testwtpc, MFRC::ElementList_NCFMAS, NormMode = "Oxygen", NormValue = 24, padding = FALSE) |>
#'           Cation2Cation(MFRC::ElementList_MnNCKFMASTCr, NormMode = "Cation", NormValue = 4, padding = TRUE)
#'
#'
Cation2Cation <- function(table, ElementList, NormMode = "Oxygen", NormValue = 24, padding = FALSE){

  if(NormMode != "Oxygen" && NormMode != "Cation"){
    stop("Set NormMode either as 'Oxygen' or 'Cation'")
  }

  CationList <- data.frame(matrix(ncol=0, nrow = nrow(table)))
  OxygenList <- data.frame(matrix(ncol=0, nrow = nrow(table)))

  ElementFilter <- ElementList[2,] %in% colnames(table)
  ElementFilterIndex <- which(ElementFilter == TRUE)

  if (padding == TRUE){

    jvector <- 1:ncol(ElementList)

  }else{

    jvector <- ElementFilterIndex

  }

  for (j in jvector){

    if (ElementFilter[j] == FALSE){

      numerator <- numeric(nrow(table))

    }else{

      numerator <- table[,ElementList[2,j]]

    }

    CationList[,ElementList[2,j]] <- numerator
    OxygenList[,ElementList[2,j]] <- numerator / as.numeric(ElementList[3,j]) * as.numeric(ElementList[4,j])

  }

  TemporaryNcol <- ncol(CationList)

  for (i in 1:nrow(CationList)){

    CationList$Total[i] <- sum(CationList[i,1:TemporaryNcol])
    CationList$Oxygen[i] <- sum(OxygenList[i,1:TemporaryNcol])

    if (NormMode == "Cation"){
      CationList[i,] <- CationList[i,] / CationList$Total[i] * NormValue
    }

    if (NormMode == "Oxygen"){
      CationList[i,] <- CationList[i,] / CationList$Oxygen[i] * NormValue
    }


  }

  return(CationList)
}
