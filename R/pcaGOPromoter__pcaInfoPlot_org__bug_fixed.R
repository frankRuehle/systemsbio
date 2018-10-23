pcaInfoPlot <- function (eData, inputType = "hgu133plus2", org = "Hs", groups, 
          printNames = TRUE, plotCI = TRUE, noProbes = 1365, GOtermsAnnotation = TRUE, 
          primoAnnotation = TRUE) 
{
  pcaObj <- pca(eData)
  if (!is.na(pcaObj$expressionData[1])) {
    chipType <- pcaObj$expressionData[[1]]
    if (chipType == "hgu133plus2") {
      inputType <- chipType
      org <- "Hs"
    }
    if (chipType == "hugene10st" | chipType == "hugene10stv1") {
      inputType <- "hugene10st"
      org <- "Hs"
    }
    if (chipType == "mouse4302") {
      inputType <- chipType
      org <- "Mm"
    }
    if (chipType == "mogene10st" | chipType == "mogene10stv1") {
      inputType <- "mogene10st"
      org <- "Mm"
    }
    if (chipType == "mogene20st") {
      inputType <- "mogene20st"
      org <- "Mm"
    }
    if (chipType == "rat2302") {
      inputType <- chipType
      org <- "Rn"
    }
  }
  probesPC1pos <- getRankedProbeIds(pcaObj, pc = 1, decreasing = TRUE)[1:noProbes]
  probesPC1neg <- getRankedProbeIds(pcaObj, pc = 1, decreasing = FALSE)[1:noProbes]
  probesPC2pos <- getRankedProbeIds(pcaObj, pc = 2, decreasing = TRUE)[1:noProbes]
  probesPC2neg <- getRankedProbeIds(pcaObj, pc = 2, decreasing = FALSE)[1:noProbes]
  if (GOtermsAnnotation) {
    GOtreePC1pos <- GOtree(probesPC1pos, inputType = inputType, org = org)
    GOtreePC1neg <- GOtree(probesPC1neg, inputType = inputType, org = org)
    GOtreePC2pos <- GOtree(probesPC2pos, inputType = inputType, org = org)
    GOtreePC2neg <- GOtree(probesPC2neg, inputType = inputType, org = org)
    GOtreeObjs <- list(GOtreePC1pos, GOtreePC1neg, GOtreePC2pos, GOtreePC2neg)
  }
  else {
    GOtreeObjs <- NA
  }
  if (primoAnnotation) {
    primoPC1pos <- primo(probesPC1pos, inputType = inputType, org = org)
    primoPC1neg <- primo(probesPC1neg, inputType = inputType, org = org)
    primoPC2pos <- primo(probesPC2pos, inputType = inputType, org = org)
    primoPC2neg <- primo(probesPC2neg, inputType = inputType, org = org)
    primoObjs <- list(primoPC1pos, primoPC1neg, primoPC2pos, primoPC2neg)
  }
  else {
    primoObjs <- NA
  }
  plot(pcaObj, groups = groups, printNames = printNames, plotCI = plotCI, 
       GOtreeObjs = GOtreeObjs, primoObjs = primoObjs)
}
