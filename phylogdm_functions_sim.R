formatsitepair_phylo <- function (bioData, bioFormat, dist = "bray", abundance = FALSE, 
          siteColumn = NULL, XColumn, YColumn, sppColumn = NULL, abundColumn = NULL, 
          sppFilter = 0, predData, distPreds = NULL, weightType = "equal", 
          custWeights = NULL, sampleSites = 1, verbose = FALSE, community_data, thetree) 
{
  if (!(is(bioData, "data.frame") | is(bioData, "matrix") | 
        is(bioData, "gdmData"))) {
    "bioData object needs to be either of class data.frame or matrix in one of\n    the acceptable formats listed in the help document."
  }
  if (is(bioData, "data.frame") | is(bioData, "matrix")) {
    bioData <- as.data.frame(bioData, stringsAsFactors = F)
  }
  if (!(is(predData, "data.frame") | is(predData, "matrix") | 
        is(predData, "RasterStack") | is(predData, "RasterLayer") | 
        is(predData, "RasterBrick"))) {
    "predData object needs to either of class data.frame, matrix, or raster."
  }
  if (is(predData, "data.frame") | is(predData, "matrix")) {
    predData <- as.data.frame(predData, stringsAsFactors = F)
  }
  if (bioFormat %in% c(1:4)) {
  }
  else {
    stop("Acceptable values for the bioFormat argument are: 1, 2, 3, or 4")
  }
  if (!(abundance == TRUE | abundance == FALSE)) {
    stop("abundance argument must be either TRUE or FALSE")
  }
  if (is.numeric(sampleSites) == FALSE | sampleSites <= 0 | 
      sampleSites > 1) {
    stop("sampleSites argument must be a number 0 < x <= 1")
  }
  if ((is.numeric(sppFilter) == FALSE & is.null(sppFilter) == 
       FALSE) | sppFilter < 0) {
    stop("sppFilter argument must be a positive integer")
  }
  if (weightType %in% c("equal", "richness", "custom")) {
  }
  else {
    stop("Acceptable values for the weightType argument are: equal, richness, or custom")
  }
  if (weightType == "custom" & is.null(custWeights) == T) {
    stop("weightType argument = 'custom', but no custWeights vector provided.")
  }
  if (bioFormat == 2 & is.null(sppColumn) == TRUE) {
    stop("Need to define sppColumn argument when bioFormat==2")
  }
  if (bioFormat == 2 & is.null(siteColumn) == TRUE) {
    if (!(is(predData, "RasterStack") | is(predData, "RasterLayer") | 
          is(predData, "RasterBrick"))) {
      stop("A siteColumn needs to be provided in either the bioData or predData objects.")
    }
  }
  if (is.null(siteColumn) == FALSE) {
    if (!is(siteColumn, "character")) {
      stop("siteColumn argument needs to be of class = 'character'.")
    }
    else if (!(siteColumn %in% colnames(bioData)) & (bioFormat == 
                                                     1 | bioFormat == 2)) {
      stop("Cannot find a match for the name of the siteColumn in the columns\n           of the bioData object.")
    }
  }
  if (bioFormat != 4) {
    if (!is(XColumn, "character")) {
      stop("XColumn argument needs to be of class 'character'.")
    }
    else if (!is(YColumn, "character")) {
      stop("YColumn argument needs to be of class 'character'.")
    }
    else if (!(XColumn %in% colnames(bioData) | XColumn %in% 
               colnames(predData))) {
      stop("XColumn not found in either the bioData or predData arguments")
    }
    else if (!(YColumn %in% colnames(bioData) | YColumn %in% 
               colnames(predData))) {
      stop("YColumn not found in either the bioData or predData arguments")
    }
  }
  if (bioFormat == 3) {
    if (weightType == "richness") {
      stop("Cannot weight by site richness when supplying the biological data\n           as a distance matrix.")
    }
    else if (nrow(bioData) != (ncol(bioData) - 1)) {
      stop("Biological dissimilarity matrix must have the same number of rows\n           and columns. Did you forget to add a column for site ID's?")
    }
  }
  for (mat in distPreds) {
    if (!is(mat, "matrix") & !is(mat, "data.frame")) {
      warning("One or more of the provided predictor distance matrices are not\n              of class 'matrix'.")
    }
  }
  if (is.null(custWeights) == FALSE & (!is(custWeights, "data.frame") & 
                                       !is(custWeights, "matrix"))) {
    stop("The argument custWeights needs to be of class 'data.frame' or 'matrix'.")
  }
  toRemove <- NULL
  removeRand <- NULL
  distData <- NULL
  if (bioFormat == 1 | bioFormat == 2) {
    if (bioFormat == 2) {
      if ((sppColumn %in% colnames(bioData))) {
      }
      else {
        stop("Cannot find sppColumn in bioData - check name?")
      }
      if (is.null(siteColumn)) {
        colnames(bioData)[which(colnames(bioData) == 
                                  XColumn)] <- "myXness"
        colnames(bioData)[which(colnames(bioData) == 
                                  YColumn)] <- "myYness"
        bioData <- transform(bioData, siteUltimateCoolness = as.numeric(interaction(bioData$myXness, 
                                                                                    bioData$myYness, drop = TRUE)))
        siteColumn <- "siteUltimateCoolness"
        colnames(bioData)[which(colnames(bioData) == 
                                  "myXness")] <- XColumn
        colnames(bioData)[which(colnames(bioData) == 
                                  "myYness")] <- YColumn
      }
      if (is.null(abundColumn)) {
        warning("No abundance column was specified, so the species data are\n                assumed to be presences.")
        bioData["reallysupercoolawesomedata"] <- 1
        abundColumn <- "reallysupercoolawesomedata"
      }
      preCastBio <- bioData
      colnames(preCastBio)[which(colnames(preCastBio) == 
                                   siteColumn)] <- "siteUltimateCoolness"
      colnames(preCastBio)[which(colnames(preCastBio) == 
                                   sppColumn)] <- "spcodeUltimateCoolness"
      castData <- dcast(preCastBio, fill = 0, siteUltimateCoolness ~ 
                          spcodeUltimateCoolness, value.var = abundColumn)
      uniqueCoords <- unique(preCastBio[which(colnames(preCastBio) %in% 
                                                c("siteUltimateCoolness", XColumn, YColumn))])
      bioData <- merge(castData, uniqueCoords, by = "siteUltimateCoolness")
      colnames(bioData)[which(colnames(bioData) == "siteUltimateCoolness")] <- siteColumn
    }
    if ((XColumn %in% colnames(bioData)) == FALSE | (YColumn %in% 
                                                     colnames(bioData) == FALSE)) {
      xCol <- which(colnames(predData) == XColumn)
      yCol <- which(colnames(predData) == YColumn)
      locs <- predData[c(xCol, yCol)]
    }
    else {
      xCol <- which(colnames(bioData) == XColumn)
      yCol <- which(colnames(bioData) == YColumn)
      locs <- bioData[c(xCol, yCol)]
    }
    if (is(predData, "RasterStack") | is(predData, "RasterLayer") | 
        is(predData, "RasterBrick")) {
      warning("When using rasters for prediction data, sites are assigned to the\n              cells in which they are located and then aggreagted as necessary (e.g.,\n              if more than one site falls in the same raster cell - common for rasters\n              with large cells).")
      cellID <- as.data.frame(raster::cellFromXY(predData, locs))
      colnames(cellID)[which(colnames(cellID) == "raster::cellFromXY(predData, locs)")] <- "cellName"
      if (nrow(cellID) == sum(is.na(cellID$cellName))) {
        stop("None of the data points provided intersect with the rasters. Double check spatial data.")
      }
      cellLocs <- as.data.frame(raster::xyFromCell(predData, cellID$cellName))
      rastBioData <- cbind(cellID, cellLocs, bioData[-c(which(colnames(bioData) %in% 
                                                                c(XColumn, YColumn)))])
      if (weightType == "custom" & !is.null(custWeights)) {
        nameTab <- unique(rastBioData[c("cellName", siteColumn)])
        tempWeightTab <- merge(x = nameTab, y = custWeights, 
                               by = siteColumn)
        siteNum <- which(colnames(tempWeightTab) == "cellName")
        custWeights <- tempWeightTab[-siteNum]
        colnames(custWeights)[1] <- siteColumn
      }
      siteNum <- which(colnames(rastBioData) == siteColumn)
      rastBioData <- rastBioData[-siteNum]
      cellNum <- which(colnames(rastBioData) == "cellName")
      bioData <- aggregate(rastBioData, rastBioData[cellNum], 
                           FUN = mean)
      bioData <- bioData[-cellNum]
      rastEx <- as.data.frame(extract(predData, bioData$cellName))
      colnames(bioData)[which(colnames(bioData) == "cellName")] <- siteColumn
      colnames(bioData)[which(colnames(bioData) == "x")] <- XColumn
      colnames(bioData)[which(colnames(bioData) == "y")] <- YColumn
      xCol <- which(colnames(bioData) == XColumn)
      yCol <- which(colnames(bioData) == YColumn)
      locs <- bioData[c(xCol, yCol)]
      siteCol <- which(colnames(bioData) == siteColumn)
      predData <- cbind(bioData[siteCol], locs, rastEx)
    }
    ##filters out sites with low species counts
    ##first isolates the species data
    siteCol <- which(colnames(bioData) == siteColumn)
    xCol <- which(colnames(bioData) == XColumn)
    yCol <- which(colnames(bioData) == YColumn)
    sppDat <- bioData[-c(siteCol, xCol, yCol)]
    sppDat[sppDat >= 1] <- 1
    sppDat[sppDat == 0] <- 0
    sppDat[is.na(sppDat)] <- 0
    sppDat[] <- apply(sppDat, 2, function(x) ifelse(x < 1 & x > 0, 1, x))
    sppTotals <- cbind.data.frame(bioData[, siteCol], apply(sppDat, 
                                                            1, function(m) {
                                                              sum(as.numeric(m))
                                                            }))
    colnames(sppTotals) <- c(siteColumn, "totals")
    filterBioDat <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] >= 
                             sppFilter)
    toRemove <- bioData[, siteCol][which(!(bioData[, siteCol] %in% 
                                             filterBioDat[, 1]))]
    names(filterBioDat)[1] <- siteColumn
    spSiteCol <- filterBioDat[1]
    bioData <- unique(merge(spSiteCol, bioData, by = siteColumn))
    if (sampleSites < 1) {
      fullSites <- bioData[, siteCol]
      randRows <- sort(sample(1:nrow(bioData), round(nrow(bioData) * 
                                                       sampleSites, 0)))
      bioData <- bioData[c(randRows), ]
      removeRand <- fullSites[which(!(fullSites %in% bioData[, 
                                                             siteCol]))]
    }
    colnames(bioData)[colnames(bioData) == siteColumn] <- "theSiteColumn"
    colnames(predData)[colnames(predData) == siteColumn] <- "theSiteColumn"
    predData <- unique(predData)
    predData <- predData[which(predData$theSiteColumn %in% 
                                 as.character(bioData$theSiteColumn)), ]
    if (weightType == "custom" & !is.null(custWeights)) {
      colnames(custWeights)[colnames(custWeights) == siteColumn] <- "theSiteColumn"
      custWeights <- custWeights[which(predData$theSiteColumn %in% 
                                         custWeights[, "theSiteColumn"]), ]
      hwap <- custWeights[, "theSiteColumn"]
      hwap <- order(hwap)
      custWeights <- custWeights[hwap, ]
      colnames(custWeights)[colnames(custWeights) == "theSiteColumn"] <- siteColumn
    }
    colnames(bioData)[colnames(bioData) == "theSiteColumn"] <- siteColumn
    colnames(predData)[colnames(predData) == "theSiteColumn"] <- siteColumn
    predSite <- which(names(predData) == siteColumn)
    bioSite <- which(names(bioData) == siteColumn)
    hwap <- predData[, predSite]
    hwap <- order(hwap)
    predData <- predData[hwap, ]
    rosetta <- bioData[, bioSite]
    rosetta <- order(rosetta)
    bioData <- bioData[rosetta, ]
    bx <- which(names(bioData) == XColumn)
    by <- which(names(bioData) == YColumn)
    sppData <- bioData[,-c(bioSite, bx, by)]
    if (abundance == F) {
      sppData[sppData >= 1] <- 1
      sppData[sppData == 0] <- 0
      sppData[is.na(sppData)] <- 0
      sppData[] <- apply(sppData, 2, function(x) ifelse(x < 1 & x > 0, 1, x))
      species_assemblages_comp_sparse <- dense2sparse(sppData)
      thetree <- drop.tip(thetree, which(is.na(match(thetree$tip.label, names(sppData)))))
      distData <- phylobeta(species_assemblages_comp_sparse, thetree)$phylo.beta.sim
      }
    else {
      sppData[is.na(sppData)] <- 0
      species_assemblages_comp_sparse <- dense2sparse(sppData)
      thetree <- drop.tip(thetree, which(is.na(match(thetree$tip.label, names(sppData)))))
      distData <- phylobeta(species_assemblages_comp_sparse, thetree)$phylo.beta.sim   }
  }
  else if (bioFormat == 3) {
    holdSite <- bioData[, which(siteColumn %in% colnames(bioData))]
    bioData <- bioData[, -which(siteColumn %in% colnames(bioData))]
    orderedData <- as.matrix(as.dist(bioData[order(holdSite), 
                                             order(holdSite)]))
    distData <- lower.tri(as.matrix(orderedData), diag = FALSE)
    distData <- as.vector(orderedData[distData])
    predData <- unique(predData)
    hwap <- predData[siteColumn][, 1]
    hwap <- order(hwap)
    predData <- predData[hwap, ]
  }
  else if (bioFormat == 4) {
    outTable <- bioData
  }
  else {
    stop(paste("bioFormat argument of '", as.character(bioFormat), 
               "' is not an accepted input value", sep = ""))
  }
  if (bioFormat != 4) {
    outTable <- as.data.frame(createsitepair(dist = distData, 
                                             spdata = bioData, envInfo = predData, dXCol = XColumn, 
                                             dYCol = YColumn, siteCol = siteColumn, weightsType = weightType, 
                                             custWeights = custWeights))
  }
  else {
    outTable <- bioData
  }
  if (length(distPreds) > 0) {
    baseMat <- distPreds[[1]]
    lapply(distPreds, function(mat, mat1) {
      if ((dim(mat1)[1] != dim(mat)[1]) & (dim(mat1)[2] != 
                                           dim(mat)[2])) {
        stop("The dimensions of your predictor matrices are not the same.")
      }
    }, mat1 = baseMat)
    holdSiteCols <- lapply(distPreds, function(dP) {
      dP[, which(siteColumn %in% colnames(dP))]
    })
    distPreds <- lapply(distPreds, function(dP) {
      dP[, -which(siteColumn %in% colnames(dP))]
    })
    distPreds <- mapply(function(dP, hSC) {
      as.matrix(as.dist(dP[order(hSC), order(hSC)]))
    }, dP = distPreds, hSC = holdSiteCols, SIMPLIFY = FALSE)
    orderSiteCols <- lapply(holdSiteCols, function(hSC) {
      hSC[order(hSC)]
    })
    rmSites <- c(toRemove, removeRand)
    if (length(rmSites) > 0) {
      rmIndex <- lapply(orderSiteCols, function(hSC, tR) {
        which((hSC %in% tR))
      }, tR = rmSites)
      distPreds <- mapply(function(mat, tR) {
        mat[-c(tR), -c(tR)]
      }, mat = distPreds, tR = rmIndex, SIMPLIFY = FALSE)
    }
    baseMat <- distPreds[[1]]
    baseMatDat <- lower.tri(as.matrix(baseMat), diag = FALSE)
    baseMatDat <- as.vector(baseMat[baseMatDat])
    if (nrow(outTable) != length(baseMatDat)) {
      stop("The dimensions of the distance predictor matrices do not match the biological data.")
    }
    for (num in 1:length(distPreds)) {
      matrixDat <- lower.tri(as.matrix(distPreds[[num]], 
                                       diag = FALSE))
      sweetFreakinPredMatrix <- as.vector(distPreds[[num]][matrixDat])
      if (ncol(outTable) > 6) {
        baseSitePair <- outTable[, 1:6]
        otherSitePair <- outTable[, 7:ncol(outTable)]
        otherNames <- colnames(otherSitePair)
        s1SitePair <- as.data.frame(otherSitePair[, 1:(ncol(otherSitePair)/2)])
        colnames(s1SitePair) <- otherNames[1:(ncol(otherSitePair)/2)]
        s2SitePair <- as.data.frame(otherSitePair[, (ncol(otherSitePair)/2 + 
                                                       1):ncol(otherSitePair)])
        colnames(s2SitePair) <- otherNames[(ncol(otherSitePair)/2 + 
                                              1):ncol(otherSitePair)]
        s1Zeros <- as.data.frame(rep(0, length(sweetFreakinPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, 
                                   sep = "")
        s2Mat <- as.data.frame(sweetFreakinPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep = "")
        outTable <- cbind(baseSitePair, s1SitePair, s1Zeros, 
                          s2SitePair, s2Mat)
      }
      else {
        s1Zeros <- as.data.frame(rep(0, length(sweetFreakinPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, 
                                   sep = "")
        s2Mat <- as.data.frame(sweetFreakinPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep = "")
        outTable <- cbind(outTable, s1Zeros, s2Mat)
      }
    }
  }
  if (verbose) {
    if (weightType[1] == "equal") {
      print("Site weighting type: Equal")
    }
    else if (weightType[1] == "custom") {
      print("Site weighting type: Custom")
    }
    else {
      print("Site weighting type: Richness")
    }
    print(paste0("Site-pair table created with ", nrow(outTable), 
                 " rows ", "(", nrow(unique(outTable[, 3:4])) + 1, 
                 " unique sites)", " and ", ncol(outTable), " columns (", 
                 (ncol(outTable) - 6)/2, " environmental variables)."))
  }
  class(outTable) <- c("gdmData", "data.frame")
  return(outTable)
}





createsitepair <- function(dist, spdata, envInfo, dXCol, dYCol, siteCol,
                           weightsType, custWeights){
  ###########################
  ##lines used to quickly test function
  #dist = distData
  #spdata = sppTab
  #envInfo = predData
  #dXCol = XColumn
  #dYCol = YColumn
  #siteCol = "site"
  #siteColumn <- siteCol
  #weightsType = weightType
  #custWeights = custWeights
  ###########################
  
  locs <- c("Long","Lat")
  XColumn <- "Long"
  YColumn <- "Lat"
  cellID <- as.data.frame(raster::cellFromXY(rasts, spdata[locs]))
  colnames(cellID)[which(colnames(cellID) == "raster::cellFromXY(rasts, spdata[locs])")] <- "cellName"
  if (nrow(cellID) == sum(is.na(cellID$cellName))) {
    stop("None of the data points provided intersect with the rasters. Double check spatial data.")
  }
  cellLocs <- as.data.frame(raster::xyFromCell(rasts, cellID$cellName))
  rastBioData <- cbind(cellID, cellLocs, spdata[-c(which(colnames(spdata) %in% 
                                                            c(XColumn, YColumn)))])
  
  weightsType <- as.character(weightsType)
  distance <- as.vector(dist)
  
  names(rastBioData)[2:3] <- c(XColumn, YColumn)
  
  #rastBioData <- rastBioData[-siteNum]
  cellNum <- which(colnames(rastBioData) == "cellName")
  bioData <- aggregate(rastBioData, rastBioData[cellNum], 
                       FUN = mean)
  bioData <- bioData[-cellNum]
  rastEx <- as.data.frame(extract(rasts, bioData$cellName))
  colnames(bioData)[which(colnames(bioData) == "cellName")] <- "site"
  colnames(bioData)[which(colnames(bioData) == "x")] <- XColumn
  colnames(bioData)[which(colnames(bioData) == "y")] <- YColumn
  xCol <- which(colnames(bioData) == XColumn)
  yCol <- which(colnames(bioData) == YColumn)
  locs <- bioData[c(xCol, yCol)]
  siteCol <- which(colnames(bioData) == siteColumn)[1]
  predData <- cbind(bioData[siteCol][1], locs, rastEx)
  ##calculates total richness = the sum of the two most diverse sites
  if(weightsType[1]=="richness"){
    sppOnly <- spdata[, -c(1,2,3)]
    sppSums <- rowSums(sppOnly)
    sppSiteSums <- cbind(spdata[1], sppSums)
    orderedSums <- sppSiteSums[order(-sppSiteSums[,2]),]
    richTotal <- orderedSums[1,2]+orderedSums[2,2]
  }
  #envInfo <- predData
  
  ##Builds index needed for site-pair table format
  s1.xCoord <- s1.yCoord <- s2.xCoord <- s2.yCoord <- NULL
  s1 <- s2 <- NULL
  

    count <- seq(length(unique(envInfo[,siteCol]))-1,1,-1)
  
  s1 <- unlist(sapply(seq(length(count),1), function(y){c(s1, rep((max(count)-y)+1, times=y))}))
  s2 <- unlist(sapply(seq(length(count),1), function(y){c(s2, (max(count)-y+2):(max(count)+1))}))
  
  
  
  if(length(s1)!=length(distance)){
    stop("The length of distance values is not the same as the expected number of rows in the site-pair table, unable to proceed.")
  }
  
  if(weightsType[1]=="equal"){
    #print("Site weighting type: Equal")
    weights <- rep(1, times=length(distance))
  }else if(weightsType[1]=="custom"){
    #print("Site weighting type: Custom")
    weights <- (custWeights[s1, "weights"] + custWeights[s2, "weights"]) / 2
  }else{
    #print("Site weighting type: Richness")
    weights <- (sppSiteSums[s1, "sppSums"] + sppSiteSums[s2, "sppSums"]) / richTotal
  }
  
  
  gdmTable <- cbind(distance, weights)
  
  ##from environmental or species table, copy coordinates for site-pair table

    s1.xCoord <- envInfo[s1, dXCol]
    s2.xCoord <- envInfo[s2, dXCol]
    s1.yCoord <- envInfo[s1, dYCol]
    s2.yCoord <- envInfo[s2, dYCol]

  
  ##sets up output table
  gdmForm <- cbind(gdmTable, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord)
  xhold <- which(names(envInfo)==dXCol)
  yhold <- which(names(envInfo)==dYCol)
  sitehold <- which(names(envInfo)=="site")
  #sitehold2 <- which(names(envInfo)=="siteUltimateCoolness")
  envInfo <- envInfo[-c(xhold, yhold, sitehold)]
  
  ##fills output table
  if(ncol(envInfo)>0){
    gdmTableFill <- cbind(gdmForm, envInfo[s1,1:ncol(envInfo)], envInfo[s2,1:ncol(envInfo)])
    names.s1 <- paste("s1.",names(envInfo[1:ncol(envInfo)]), sep="")
    names.s2 <- paste("s2.",names(envInfo[1:ncol(envInfo)]), sep="")
    colnames(gdmTableFill) <- c(colnames(gdmTableFill)[1:6], names.s1, names.s2)
  }else{
    gdmTableFill <- gdmForm
  }
  
  ##returns results
  return(gdmTableFill)
}

