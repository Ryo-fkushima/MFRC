#' @title Recalculation for omphacite (NCFMAS, Al_T is omitted, Ca + Mg + Fe in M2 = Mg + Fe in M1, Mg/Fe equipartition)
#' @description \code{Omp_AlT_Omit}: For omphacite (NCFMAS, Al_T is omitted)
#'
#' @param CationTable (dataframe) Atomic proportion data
#'
#' @export
#' @examples
#' Mass2Cation(MFRC::testwtpc, MFRC::ElementList_MnNCKFMASTCr) |> Omp_AlT_Omit()
#'

Omp_AlT_Omit <- function(CationTable){

  CT <- Cation2Cation(CationTable, MFRC::ElementList_NCFMAS, NormMode = "Oxygen", NormValue = 6, padding = TRUE)

  Result <- data.frame(matrix(ncol=0, nrow = nrow(CT)))

  for (i in 1:nrow(CT)){

    Fe3 <- max(CT[i, "Na"] - CT[i, "Al"], 0)

    ccc <- (CT[i, "Fe2"] - Fe3) / CT[i, "Mg"] # Fe/Mg in M1, M2 site

    Result$Si_T[i] <- CT[i, "Si"]
    Result$Al_M1[i] <- CT[i, "Al"]

    Result$Fe2_M1[i] <- 0.5* ccc * (CT[i, "Mg"] + CT[i, "Ca"]/(ccc + 1))
    Result$Fe3_M1[i] <- Fe3
    Result$Mg_M1[i] <- Result$Fe2_M1[i] / ccc

    Result$Fe2_M2[i] <- max((0.5* ccc * (CT[i, "Mg"] - CT[i, "Ca"]/(ccc + 1))), 0)
    Result$Mg_M2[i] <- Result$Fe2_M2[i] / ccc

    Result$Ca_M2[i] <- CT[i, "Ca"]
    Result$Na_M2[i] <- CT[i, "Na"]

    Result$Total[i] <- sum(Result[i,1:9])

    Result$Oxygen[i] = 2 * (Result$Si_T[i]) + 1.5 * (Result$Al_M1[i] + Result$Fe3_M1[i]) +
      1 * (Result$Fe2_M1[i] + Result$Mg_M1[i] + Result$Fe2_M2[i] + Result$Mg_M2[i] +
             Result$Ca_M2[i])+ 0.5 * (Result$Na_M2[i])

    Result[i,] <- Result[i,] / Result$Oxygen[i] * 6


  }


  return(Result)

}
