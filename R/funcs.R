#' Calculates the proportion of a species that could migrate from a cell to its
#' neighbours
#'
#' Calculates the proportion of a species that could migrate from a cell to its
#' neighbours given a raster, the maximum dispersal distance, and the number of
#' times the time-lapse will have
#' @param Raster a raster with the space where the species will be inhabiting
#' @param Distance the maximum dispersal distance of the species
#' @param Time the number of time-slices to be used
#' @return a dataframe with the cell id from, the cell id to, and beta, that is
#' the proportion of individuals that will go from cell from to cell to
#' @examples
#' data("r")
#' DistConect(Raster = r, Distance = 1000000, Time = 7)
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom gdistance accCost
#' @importFrom gdistance geoCorrection
#' @importFrom gdistance transition
#' @importFrom magrittr "%>%"
#' @importFrom raster ncell
#' @importFrom raster xyFromCell
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Javier Fajardo <javierfajnolla@gmail.com >
#' @export

DistConect<- function(Raster, Distance, Time = 7){
  #First we make a transition layer with the function transition from gdistance
  h16  <- transition(Raster, transitionFunction=function(x){1},16,symm=FALSE)
  #Then geocorrect for projection
  h16   <- geoCorrection(h16, scl=FALSE)
  #Since transition layers work with XY rather than IDs, get a matrix of XY coordinates
  B <- xyFromCell(Raster, cell = 1:ncell(Raster))
  #This nested loop is where the Bottle neck is
  #Start a list
  connections <- list()
  #For each pair of cells in B
  accCost2 <- function(x, fromCoords) {

    fromCells <- cellFromXY(x, fromCoords)
    tr <- transitionMatrix(x)
    tr <- rBind(tr, rep(0, nrow(tr)))
    tr <- cBind(tr, rep(0, nrow(tr)))
    startNode <- nrow(tr)
    adjP <- cbind(rep(startNode, times = length(fromCells)), fromCells)
    tr[adjP] <- Inf
    adjacencyGraph <- graph.adjacency(tr, mode = "directed", weighted = TRUE)
    E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
    return(shortest.paths(adjacencyGraph, v = startNode, mode = "out")[-startNode])
  }
  for (i in 1:nrow(B)){
    #Create a temporal raster for each row with the distance from cell xy to all other cells
    temp <- accCost2(h16,B[i,])
    index <- which(temp < Distance)
    connections[[i]] <- cbind(i, index, temp[index])
  }
  #Get everything together as a large data frame
  connections <- do.call("rbind", connections)
  connections <- as.data.frame(connections)
  colnames(connections) <- c("from", "to", "dist")
  connections$Beta <- exp(-(connections$dist/max(connections$dist)))
  b <- connections %>% group_by(from) %>% summarize(TotalBeta = sum(Beta))
  connections <-merge(connections, b)
  connections$beta <-connections$Beta /connections$TotalBeta
  connections<- dplyr::filter(connections, beta > quantile(beta, 0.05))
  connections <- connections[,c(1,2,6)]
  n <- nrow(connections)
  connections <- do.call("rbind", replicate((Time), connections, simplify = FALSE))
  connections$Time <- rep(c(0:(Time-1)), each =n)
  return(connections)
}

#' Generates an AMPL dat file from a Stack
#'
#' Generates an AMPL dat file from a Stack in which each file is the projection
#' of a species distribution model into a time slice
#' @param Stack a raster with the space where the species will be inhabiting
#' @param maxalpha the maximum rate of change of the population in optimal
#' coditions, defaults in 10
#' @param maxbiomass the maximum initial biomass of the population in optimal
#' coditions, defaults in 2
#' @param maxcapacidad the maximum biomass of the population in optimal coditions,
#' defaults in 10
#' @param name the name of the .dat file that will be exported, defaults in
#' stack
#' @param Dist the maximum dispersal distance of the species modeled in the
#' stack
#' @return exports a .dat file to feed the AMPL model
#' @examples
#' data("univariate")
#' RasterToAmplDat(Stack = univariate, Threshold = 0.5)
#' @importFrom raster nlayers
#' @importFrom raster values
#' @importFrom tidyr spread
#' @importFrom tidyr unite_
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Javier Fajardo <javierfajnolla@gmail.com >
#' @export

RasterToAmplDat <- function(Stack, maxalpha = 10, maxbiomass = 2, maxcapacidad = 10, name = "Stack", Dist = 1000000, Threshold){

  TempStack <- Stack
  values(TempStack)[values(TempStack) < Threshold] = 0
  values(TempStack)[values(TempStack) >= Threshold] = 1

  TempRaster <- sum(TempStack)

  TempRaster[values(TempRaster) > 0] = 1
  TempRaster[values(TempRaster) == 0] = NA

  Stack <- Stack*TempRaster

  Alpha <- list()
  for (i in 1:nlayers(Stack)){
    temp <- data.frame(Alpha = values(Stack[[i]]), ID = 1:length(values(Stack[[i]])), Time = i)
    Alpha[[i]] <- temp[complete.cases(temp),]
  }

  Alpha <-  do.call(rbind, Alpha)
  Alpha$Alpha <- Alpha$Alpha/max(Alpha$Alpha)*maxalpha


  Biomasa <- data.frame(Biomasa = (values(Stack[[1]])/max(values(Stack[[1]]), na.rm = T))*maxbiomass, ID = 1:length(values(Stack[[1]])))
  Biomasa <-  Biomasa[complete.cases(Biomasa),]
  Biomasa$Biomasa <- round(Biomasa$Biomasa, 4)

  Capacidad <- list()
  for (i in 1:nlayers(Stack)){
    temp <- data.frame(Capacidad = values(Stack[[i]]), ID = 1:length(values(Stack[[i]])), Time = i)
    Capacidad[[i]] <- temp[complete.cases(temp),]
  }

  Capacidad <-  do.call(rbind, Capacidad)
  Capacidad$Capacidad <- Capacidad$Capacidad/max(Capacidad$Capacidad)*maxcapacidad

  Nodos <- merge(Alpha, Capacidad)

  Alphas <- spread(Alpha, key = Time, value = Alpha)
  Alphas$ID <- paste("[", Alphas$ID, ",*]", sep = "")
  Alphas$T1 <- 0
  Alphas$T2 <- 1
  Alphas$T3 <- 2
  Alphas$T4 <- 3
  Alphas$T5 <- 4
  Alphas$T6 <- 5
  Alphas$T7 <- 6
  Alphas$T8 <- 7
  Alphas <- Alphas[,c(1,10,2,11,3,12,4,13,5,14,6,15,7,16,8,17,9)]
  Alphas$line <- "\n"

  Biomasas <- Biomasa[,c(2,1)]
  Biomasas$line <- "\n"

  Capacidades <- spread(Capacidad, key = Time, value = Capacidad)
  Capacidades$ID <- paste("[", Capacidades$ID, ",*]", sep = "")
  Capacidades$T1 <- 0
  Capacidades$T2 <- 1
  Capacidades$T3 <- 2
  Capacidades$T4 <- 3
  Capacidades$T5 <- 4
  Capacidades$T6 <- 5
  Capacidades$T7 <- 6
  Capacidades$T8 <- 7
  Capacidades <- Capacidades[,c(1,10,2,11,3,12,4,13,5,14,6,15,7,16,8,17,9)]
  Capacidades$line <- "\n"

  Cost <- data.frame(ID = Biomasas$ID, Cost = 1, line = "\n")

  Beta <- DistConect(Stack[[1]], Distance = Dist, Time = nlayers(Stack))
  temp <-  split(Beta, Beta$Time)
  Betas <- do.call(cbind, lapply(1:length(temp), function(i){
    if (i == 1){
      setNames(data.frame(paste("[",temp[[i]][["from"]], ",", temp[[i]][["to"]], ",*","]", sep = ""), temp[[i]]["Time"], temp[[i]]["beta"]),
               c("V", paste("T", i, sep = ""), i-1))
    } else if (i == length(temp)){
      setNames(data.frame(temp[[i]]["Time"], temp[[i]]["beta"], rep("\n", NROW(temp[[i]]))),
               c(paste("T", i, sep = ""), i-1, "line"))
    } else {
      setNames(data.frame(temp[[i]]["Time"], temp[[i]]["beta"]),
               c(paste("T", i, sep = ""), i-1))
    }
  }))

  sink(paste0(name, ".dat"))
  cat(c("set V :=", unique(Alpha$ID), ";"))
  cat("\n")
  cat("\n")
  cat(c("set E :=", paste0("(",unique(unite_(Beta, col = "V", sep = ",", from = c("from", "to"))$V), ")"), ";"))
  cat("\n")
  cat("param T:= 7;")
  cat("\n")
  cat("param alpha :=")
  cat("\n")
  cat(do.call(paste, Alphas))
  cat(";")
  cat("\n")
  cat("param beta :=")
  cat("\n")
  cat(do.call(paste, Betas))
  cat(";")
  cat("\n")
  cat("param b0 :=")
  cat("\n")
  cat(do.call(paste, Biomasas))
  cat(";")
  cat("\n")
  cat("param u :=")
  cat("\n")
  cat(do.call(paste, Capacidades))
  cat(";")
  cat("\n")
  cat("param c :=")
  cat("\n")
  cat(do.call(paste, Cost))
  cat(";")
  cat("\n")
  sink()
  return(list(Nodos = Nodos, Biomasa = Biomasa, Alphas = Alphas, Alpha = Alpha))
}

#' Generates an AMPL dat file from a Stack
#'
#' Generates an AMPL dat file from a Stack in which each file is the projection
#' of a species distribution model into a time slice
#' @param Stack a raster with the space where the species will be inhabiting
#' @param Distance the maximum dispersal distance of the species modeled in the
#' stack
#' @param Time number of time slices present in the Stack
#' @param name the name of the .dat file that will be exported, defaults in
#' stack
#' @param Threshold minimum value in the model to allow the species to exist
#' @return exports a .dat file to feed the AMPL model
#' @examples
#' data("univariate")
#' RtoQuadAmplDat(Stack = univariate, Distance = 1000000, Threshold = 0.5,
#' name = "TeSt")
#' @importFrom gdistance accCost
#' @importFrom gdistance geoCorrection
#' @importFrom gdistance transition
#' @importFrom raster ncell
#' @importFrom raster nlayers
#' @importFrom raster values
#' @importFrom raster xyFromCell
#' @importFrom tidyr unite_
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Javier Fajardo <javierfajnolla@gmail.com >
#' @export


RtoQuadAmplDat <- function(Stack, Distance, Threshold, name){
  accCost2 <- function(x, fromCoords) {

    fromCells <- cellFromXY(x, fromCoords)
    tr <- transitionMatrix(x)
    tr <- rBind(tr, rep(0, nrow(tr)))
    tr <- cBind(tr, rep(0, nrow(tr)))
    startNode <- nrow(tr)
    adjP <- cbind(rep(startNode, times = length(fromCells)), fromCells)
    tr[adjP] <- Inf
    adjacencyGraph <- graph.adjacency(tr, mode = "directed", weighted = TRUE)
    E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
    return(shortest.paths(adjacencyGraph, v = startNode, mode = "out")[-startNode])
  }

  values(Stack)[values(Stack) < Threshold] = 0
  values(Stack)[values(Stack) >= Threshold] = 1

  Suitability <- list()
  for (i in 1:nlayers(Stack)){
    temp <- data.frame(Suitability = values(Stack[[i]]), ID = 1:length(values(Stack[[i]])), Time = i-1)
    Suitability[[i]] <- temp[complete.cases(temp),]
  }

  Suitability <- do.call("rbind", Suitability)

  Raster <- sum(Stack)

  Raster[values(Raster) > 0] = 1
  Raster[values(Raster) == 0] = NA

  h16  <- transition(Raster, transitionFunction=function(x){1},16,symm=FALSE)

  h16   <- geoCorrection(h16, scl=FALSE)

  ID <-c(1:ncell(Raster))[!is.na(values(Raster))]

  B <- xyFromCell(Raster, cell = ID)

  connections <- list()
  #For each pair of cells in B
  for (i in 1:nrow(B)){
    #Create a temporal raster for each row with the distance from cell xy to all other cells
    temp <- accCost2(h16,B[i,])
    index <- which(temp < Distance)
    connections[[i]] <- cbind(ID[i], index, temp[index])
  }
  #Get everything together as a large data frame
  connections <- do.call("rbind", connections)
  connections <- as.data.frame(connections)
  colnames(connections) <- c("from", "to", "dist")
  temp <-  split(Suitability, Suitability$Time)
  Suitability <- do.call(cbind, lapply(1:length(temp), function(i){
    if (i == 1){
      setNames(data.frame(paste("[",temp[[i]][["ID"]], ",*","]", sep = ""), temp[[i]]["Time"], temp[[i]]["Suitability"]),
               c("V", paste("T", i, sep = ""), i-1))
    } else if (i == length(temp)){
      setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Suitability"], rep("\n", NROW(temp[[i]]))),
               c(paste("T", i, sep = ""), i-1, "line"))
    } else {
      setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Suitability"]),
               c(paste("T", i, sep = ""), i-1))
    }
  }))

  sink(paste0(name, ".dat"))
  cat(c("set V :=", unique(unique(connections$to), unique(connections$to)), ";"))
  cat("\n")
  cat("\n")
  cat(c("set E :=", paste0("(",unique(unite_(connections, col = "V", sep = ",", from = c("from", "to"))$V), ")"), ";"))
  cat("\n")
  cat("\n")
  cat(paste0("param T:= ", nlayers(Stack),";"))
  cat("\n")
  cat("\n")
  cat("param u :=")
  cat("\n")
  cat(do.call(paste, Suitability))
  cat(";")
  cat("\n")
  sink()
  return(list(connections = connections, Suitability = Suitability))
}
