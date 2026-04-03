#' @title Converting atomic proportions into wt%
#' @description \code{Cation2Cation} converts atomic proportions into wt%.
#'
#' @param table (dataframe) Atomic proportion data
#' @param ElementList (dataframe) Information of considered elements, cation/oxygen numbers per unit formula.
#' Use MFRC::ElementList_all or its subset (e.g., MFRC::ElementList_MnNCKFMASTCr, MFRC::ElementList_NCFMAS).
#'
#' @param NormValue (num) Normalization value for total (default = 100)
#' @param padding (logical) TRUE for padding concentrations for element absent in the input data (default = FALSE)
#'
#'
#' @export
#' @examples
#' Mass2Cation(MFRC::testwtpc, MFRC::ElementList_NCFMAS, NormMode = "Oxygen", NormValue = 24, padding = FALSE) |>
#'           Cation2MassPc(MFRC::ElementList_MnNCKFMASTCr, NormValue = 100, padding = TRUE)
#'
#'
Cation2MassPc <- function(table, ElementList, NormValue = 100, padding = FALSE){

  MassPcList <- data.frame(matrix(ncol=0, nrow = nrow(table)))

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


    MassPcList[,ElementList[1,j]] <- numerator * MFRC::MassList[,ElementList[1,j]] / as.numeric(ElementList[3,j])

  }

  TemporaryNcol <- ncol(MassPcList)

  for (i in 1:nrow(MassPcList)){

    MassPcList$Total[i] <- sum(MassPcList[i,1:TemporaryNcol])

    MassPcList[i,] <- MassPcList[i,] / MassPcList$Total[i] * NormValue

  }

  return(MassPcList)


}
