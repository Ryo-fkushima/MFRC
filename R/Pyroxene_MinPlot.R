#' @title Recalculation for pyroxene (method from MinPlot)
#' @description \code{Pyroxene_MinPlot} calculates site fractions for pyroxene with the method used in MinPlot (Walters, 2022).
#'
#' @param CationTable (dataframe) Atomic proportion data
#'
#' @export
#' @examples
#' Mass2Cation(MFRC::testwtpc, MFRC::ElementList_MnNCKFMASTCr) |> Pyroxene_MinPlot()
#' @references Walters, J. B. (2022). A mineral formula recalculation and plotting program for electron probe microanalysis.
#' Mineralogia, 53(1), 51-66.
#'
Pyroxene_MinPlot <- function(CationTable){

  CT <- Cation2Cation(CationTable, MFRC::ElementList_MnNCKFMASTCr, NormMode = "Cation", NormValue = 4, padding = TRUE)

  Result <- data.frame(matrix(ncol=0, nrow = nrow(CT)))

  for (i in 1:nrow(CT)){

    Oxygen_def <- max(0, (6 - CT[i, "Oxygen"]))
    Fe3 <- min(2 * Oxygen_def, CT[i, "Fe2"])
    XMg <- CT[i, "Mg"] / (CT[i, "Mg"] + CT[i, "Fe2"] - Fe3)

    # T site

    Result$Si_T[i] <- min(2, CT[i, "Si"])
    Result$Al_T[i] <- min((2 - Result$Si_T[i]), CT[i, "Al"])
    Result$Fe3_T[i] <- min((2 - Result$Si_T[i] - Result$Al_T[i]), Fe3)
    Result$Sum_T[i] <- Result$Si_T[i] + Result$Al_T[i] + Result$Fe3_T[i]

    # M1 site

    Result$Al_M1[i] <-  CT[i, "Al"] - Result$Al_T[i]
    Result$Ti_M1[i] <-  CT[i, "Ti"]
    Result$Cr_M1[i] <-  CT[i, "Cr"]
    Result$Fe3_M1[i] <-  Fe3 - Result$Fe3_T[i]
    Result$Mn_M1[i] <-  CT[i, "Mn2"]

    MgFe_M1_by_subtraction <- 1 - (Result$Al_M1[i] + Result$Ti_M1[i] + Result$Cr_M1[i] +
                                     Result$Fe3_M1[i] + Result$Mn_M1[i])

    Result$Mg_M1[i] <- min((XMg * MgFe_M1_by_subtraction), CT[i, "Mg"])

    Fe_M1_by_subtraction <- MgFe_M1_by_subtraction - Result$Mg_M1[i]

    if (Fe_M1_by_subtraction < 0){

      Result$Fe2_M1[i] <- 0

    }else{

      Result$Fe2_M1[i] <- min(Fe_M1_by_subtraction, (CT[i, "Fe2"] - Fe3))
    }

    Result$Sum_M1[i] <- Result$Al_M1[i] + Result$Ti_M1[i] + Result$Cr_M1[i] +
      Result$Fe3_M1[i] + Result$Mn_M1[i] + Result$Mg_M1[i] + Result$Fe2_M1[i]


    # M2 site

    Result$Mg_M2[i] <- max(0, (CT[i, "Mg"] - Result$Mg_M1[i]))
    Result$Fe2_M2[i] <- max(0, (CT[i, "Fe2"] - Fe3 - Result$Fe2_M1[i]))
    Result$Ca_M2[i] <-  CT[i, "Ca"]
    Result$Na_M2[i] <-  CT[i, "Na"]
    Result$K_M2[i] <-  CT[i, "K"]
    Result$Sum_M2[i] <- Result$Mg_M2[i] + Result$Fe2_M2[i] + Result$Ca_M2[i] + Result$Na_M2[i] + Result$K_M2[i]



  }


  return(Result)

}
