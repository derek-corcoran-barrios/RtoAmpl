library(gdistance)
library(RtoAmpl)

DistMax <- function(m){3.31*(m^0.65)*10 * 1000}
Density <- function(m) {m*4.23*(m^-0.75)}

DistConect2 <- function(Raster, m, Time = 7){
  #First we make a transition layer with the function transition from gdistance
  h16  <- transition(Raster, transitionFunction=function(x){1},16,symm=FALSE)
  #Then geocorrect for projection
  h16   <- geoCorrection(h16, scl=FALSE)
  #Since transition layers work with XY rather than IDs, get a matrix of XY coordinates
  B <- xyFromCell(Raster, cell = 1:ncell(Raster))
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
    index <- which(temp < DistMax(m))
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
  #connections<- dplyr::filter(connections, beta > quantile(beta, 0.05))
  connections <- connections[,c(1,2,6)]
  n <- nrow(connections)
  connections <- do.call("rbind", replicate((Time), connections, simplify = FALSE))
  connections$Time <- rep(c(0:(Time-1)), each =n)
  return(connections)
}

library(RtoAmpl)

data(Cost)
costlayer <- Cost
data("univariate")
data("bivariate")

extent(Cost) <- c(-2.5,2.5,-2.5,2.5)
extent(univariate) <- c(-2.5,2.5,-2.5,2.5)
extent(bivariate) <- c(-2.5,2.5,-2.5,2.5)


Stacklist <- list(univariate, bivariate)
names(Stacklist) <- c("SPA", "SPB")
Threshold = 0.0
mass <- c(25,50)
alphamax <- function(m){(3.5 * 7.5  * (10^-4) * (m^-0.25) )* 3650}

multisppBiomass <- function(Stacklist, name = "Stack", Threshold, costlayer, mass = c(50,50), N = 4169){

  Masklayer <- costlayer
  values(Masklayer) <- ifelse(is.na(values(Masklayer)), NA, 1)
  TempRaster <- list()
  for (i in 1:length(Stacklist)){
    Stacklist[[i]] <- Stacklist[[i]] * Masklayer
    TempStack <- Stacklist[[i]]

    values(TempStack)[values(TempStack) < Threshold] = 0
    values(TempStack)[values(TempStack) >= Threshold] = 1

    TempRaster[i] <- sum(TempStack)
  }

  TempRaster <- do.call("sum", TempRaster)
  TempRaster[values(TempRaster) > 0] = 1
  TempRaster[values(TempRaster) == 0] = NA
  for (i in 1:length(Stacklist)){
    Stacklist[[i]] <- Stacklist[[i]]*TempRaster
  }


  Alphas <- list()
  for(j in 1:length(Stacklist)){
    Alpha <- list()
    for (i in 1:nlayers(Stacklist[[j]])){
      temp <- data.frame(Alpha = values(Stacklist[[j]][[i]]), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i-1)
      Alpha[[i]] <- temp[complete.cases(temp),]
    }
    Alphas[[j]]<- do.call("rbind", Alpha)
    Alphas[[j]]$Spp <- names(Stacklist)[j]
    Alphas[[j]]$Alpha <- Alphas[[j]]$Alpha*alphamax(mass[j])
  }

  Alpha <-  do.call(rbind, Alphas)
  s <- Alpha %>% group_by(ID) %>% summarise(SUMA = sum(Alpha)) %>% filter(SUMA > 0)
  Alpha <- Alpha[Alpha$ID %in% s$ID,]
  Spps <- unique(Alpha$Spp)
  Nodos <- unique(Alpha$ID)

  Alphas <- list()
  for (i in Spps){
    Alphas[[i]] <- dplyr::filter(Alpha, Spp == i)

    temp <-  split(Alphas[[i]], Alphas[[i]]$Time)
    Alphas[[i]] <- do.call(cbind, lapply(1:length(temp), function(i){
      if (i == 1){
        setNames(data.frame(paste("[",temp[[i]][["Spp"]],",",temp[[i]][["ID"]], ",*","]", sep = ""), temp[[i]]["Time"], temp[[i]]["Alpha"]),
                 c("V", paste("T", i, sep = ""), i-1))
      } else if (i == length(temp)){
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Alpha"], rep("\n", NROW(temp[[i]]))),
                 c(paste("T", i, sep = ""), i-1, "line"))
      } else {
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Alpha"]),
                 c(paste("T", i, sep = ""), i-1))
      }
    }))
  }

  Alpha <-do.call("rbind", Alphas)

  Biomasa <- list()

  for (i in 1:length(Stacklist)){
    Biomasa[[i]] <- data.frame(Sp = rep(Spps[i], times = length(values(Stacklist[[i]][[1]]))), ID = 1:length(values(Stacklist[[i]][[1]])), Biomasa = values(Stacklist[[i]][[1]])*Density(mass[i])*cellStats(area(costlayer), mean))
    Biomasa[[i]] <-  Biomasa[[i]][complete.cases(Biomasa[[i]]),]
    Biomasa[[i]]$Biomasa <- round(Biomasa[[i]]$Biomasa, 4)
    Biomasa[[i]]$Sp <- Spps[i]
  }

  Biomasa <- do.call(rbind, Biomasa)
  Biomasa$line <- "\n"
  Biomasa$ID  <-  paste0("[",Biomasa$Sp, ",", Biomasa$ID, "]")
  Biomasa <- Biomasa[,-1]

  Capacidades <- list()
  for(j in 1:length(Stacklist)){
    Capacidad <- list()
    for (i in 1:nlayers(Stacklist[[j]])){
      temp <- data.frame(Capacidad = values(Stacklist[[j]][[i]])*Density(mass[j])*cellStats(area(costlayer), mean), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i-1)
      Capacidad[[i]] <- temp[complete.cases(temp),]
    }
    Capacidades[[j]]<- do.call("rbind", Capacidad)
    Capacidades[[j]]$Spp <- names(Stacklist)[j]
  }

  Capacidad <-  do.call(rbind, Capacidades)
  Capacidad <- Capacidad[Capacidad$ID %in% s$ID,]

  Capacidades <- list()
  for (i in Spps){
    Capacidades[[i]] <- dplyr::filter(Capacidad, Spp == i)

    temp <-  split(Capacidades[[i]], Capacidades[[i]]$Time)
    Capacidades[[i]] <- do.call(cbind, lapply(1:length(temp), function(i){
      if (i == 1){
        setNames(data.frame(paste("[",temp[[i]][["Spp"]],",",temp[[i]][["ID"]], ",*","]", sep = ""), temp[[i]]["Time"], temp[[i]]["Capacidad"]),
                 c("V", paste("T", i, sep = ""), i-1))
      } else if (i == length(temp)){
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Capacidad"], rep("\n", NROW(temp[[i]]))),
                 c(paste("T", i, sep = ""), i-1, "line"))
      } else {
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["Capacidad"]),
                 c(paste("T", i, sep = ""), i-1))
      }
    }))
  }

  Capacidad <-do.call("rbind", Capacidades)


  cost <-data.frame(ID = Nodos,cost = values(costlayer)[Nodos])
  cost$ID <- paste0("[", cost$ID, "]")
  cost$line <- "\n"

  BF <- data.frame(Spp = Spps, BF = mass* N, Space = "\n")


  Beta <- list()
  for(i in 1:length(Stacklist)){
    Beta[[i]] <- DistConect2(Stacklist[[i]][[1]], m = mass[i], Time = nlayers(Stacklist[[i]]))
    Beta[[i]] <- dplyr::filter(Beta[[i]], from %in% Nodos, to %in% Nodos)
    Beta[[i]]$Sp <- Spps[i]
  }
  Beta[[1]] <-right_join(Beta[[1]], Beta[[2]][,c(1,2,4)])
  Beta[[1]]$beta <- ifelse(is.na(Beta[[1]]$beta), 0, Beta[[1]]$beta)
  Beta[[1]]$Sp <- ifelse(is.na(Beta[[1]]$Sp), names(Stacklist)[1], Beta[[1]]$Sp)

  Beta <- do.call(rbind, Beta)

  Beta <- dplyr::filter(Beta, Time != 7)

  temp <-  split(Beta, Beta$Time)
  Betas <- do.call(cbind, lapply(1:length(temp), function(i){
    if (i == 1){
      setNames(data.frame(paste("[",temp[[i]][["Sp"]], ",",temp[[i]][["from"]], ",", temp[[i]][["to"]], ",*","]", sep = ""), temp[[i]]["Time"], temp[[i]]["beta"]),
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
  cat(c("set V :=", Nodos, ";"))
  cat("\n")
  cat(c("set Sp :=", names(Stacklist), ";"))
  cat("\n")
  cat(c("set E :=", paste0("(",unique(unite_(Beta, col = "V", sep = ",", from = c("from", "to"))$V), ")"), ";"))
  cat("\n")
  cat("param T:= 7;")
  cat("\n")
  cat("param alpha :=")
  cat("\n")
  cat(do.call(paste, Alpha))
  cat(";")
  cat("\n")
  cat("param beta :=")
  cat("\n")
  cat(do.call(paste, Betas))
  cat(";")
  cat("\n")
  cat("param b0 :=")
  cat("\n")
  cat(do.call(paste, Biomasa))
  cat(";")
  cat("\n")
  cat("param u :=")
  cat("\n")
  cat(do.call(paste, Capacidad))
  cat(";")
  cat("\n")
  cat("param bf := ")
  cat("\n")
  cat(do.call(paste, BF))
  cat(";")
  cat("\n")
  cat("param c :=")
  cat("\n")
  cat(do.call(paste, cost))
  cat(";")
  cat("\n")
  sink()
  return(list(Nodos = Nodos, Biomasa = Biomasa, Alphas = Alphas, Alpha = Alpha))
}


setwd("/home/derek/Documents/PostdocPablo/AMPL")

library(beepr)

Ns <- seq(from = 100, to = 7300, length.out =  73)
for (i in 1:length(Ns)){
  multisppBiomass(Stacklist = Stacklist, costlayer = Cost, Threshold = 0.5, name = paste0("Multispp","_", i), N = Ns[i], mass = mass)
  print(i)
  beep(1)
}

SOLS1 <- list.files(pattern = "Multispp_[0-9].txt")
SOLS2 <-  list.files(pattern = "Multispp_[0-9][0-9].txt")
SOLS <- c(SOLS1, SOLS2)
library(readr)

bla <- list()
for (i in 1:length(SOLS)){
  DF <- read_delim(SOLS[i], " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

  colnames(DF) <- c("ID", "z")

  temp <- univariate[[1]]

  values(temp) <- NA
  values(temp)[DF$ID] <- DF$z
  bla[[i]] <- temp
}

Costs1 <-  list.files(pattern = "Multispp_[0-9]cost.txt")
Costs2 <-  list.files(pattern = "Multispp_[0-9][0-9]cost.txt")
Costs <- c(Costs1, Costs2)
library(readr)
DF <- list()
for(i in 1:length(Costs)){
  DF[[i]] <- read_delim(Costs[i],  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(DF[[i]]) <- c("bf", "cost")
  DF[[i]]$SP <- c("A", "B")
}

DF <- do.call(rbind, DF)

DF <- DF %>% filter(SP == "A") %>% mutate(N = bf/mass[1])

library(ggplot2)
library(dplyr)
DF <- filter(DF, cost > 0)
DF$Buy <- "Continuous"
ggplot(DF, aes(x = cost, y = N)) + geom_line() + theme_classic() + geom_hline(yintercept = 4169, lty = 2, col = "red")

animation::saveGIF(for(i in 1:length(bla)){plot(bla[[i]], colNA = "black", main = paste("Target N =", round(DF$N[i],0)))}, movie.name = "Continuous.gif")


#####################################################################################
################################Binary#####################################
#####################################################################################

SOLS1 <- list.files(pattern = "MultisppBIN_[0-9].txt")
SOLS2 <-  list.files(pattern = "MultisppBIN_[0-9][0-9].txt")
SOLS <- c(SOLS1, SOLS2)
library(readr)

blaBIN <- list()
for (i in 1:length(SOLS)){
  DFBIN <- read_delim(SOLS[i], " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

  colnames(DFBIN) <- c("ID", "z")

  temp <- univariate[[1]]

  values(temp) <- NA
  values(temp)[DFBIN$ID] <- DFBIN$z
  blaBIN[[i]] <- temp
}

Costs1 <-  list.files(pattern = "MultisppBIN_[0-9]cost.txt")
Costs2 <-  list.files(pattern = "MultisppBIN_[0-9][0-9]cost.txt")
Costs <- c(Costs1, Costs2)

DFBIN <- list()
for(i in 1:length(Costs)){
  DFBIN[[i]] <- read_delim(Costs[i],  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(DFBIN[[i]]) <- c("bf", "cost")
  DFBIN[[i]]$SP <- c("A", "B")
}

DFBIN <- do.call(rbind, DFBIN)

DFBIN <- DFBIN %>% filter(SP == "A") %>% mutate(N = bf/mass[1])

library(ggplot2)
library(dplyr)
DFBIN <- filter(DFBIN, cost > 0)
DFBIN$Buy <- "Binary"

ggplot(DFBIN, aes(x = cost, y = N)) + geom_line() + theme_classic() + geom_hline(yintercept = 4169, lty = 2, col = "red")

DF <- rbind(DF, DFBIN)
ggplot(DF, aes(x = cost, y = N)) + geom_line(aes(color = Buy)) + theme_classic()+ geom_hline(yintercept = 4169, lty = 2, col = "red")

animation::saveGIF(for(i in 1:length(blaBIN)){plot(blaBIN[[i]], colNA = "black", main = paste("Target N =", round(DF$N[i])))}, movie.name = "Binary.gif")


###############################################################################
##################Accuracy##############################################
###############################################################################
library(RtoAmpl)

setwd("/home/derek/Documents/PostdocPablo/AMPL")

SOLS1 <- list.files(pattern = "Multispp_[0-9].txt")
SOLS2 <-  list.files(pattern = "Multispp_[0-9][0-9].txt")
SOLS <- c(SOLS1, SOLS2)
library(readr)
library(dplyr)

DF <- list()
blaRoundC <- list()
for (i in 1:length(SOLS)){
  DF[[i]] <- read_delim(SOLS[i], " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)
  colnames(DF[[i]]) <- c("ID", "Index")
  DF[[i]]$BF <- i
}

DF <- do.call(rbind, DF)


SOLS1 <- list.files(pattern = "MultisppBIN_[0-9].txt")
SOLS2 <-  list.files(pattern = "MultisppBIN_[0-9][0-9].txt")
SOLS <- c(SOLS1, SOLS2)
library(readr)

DFBIN <- list()

for (i in 1:length(SOLS)){
  DFBIN[[i]] <- read_delim(SOLS[i], " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)
  colnames(DFBIN[[i]]) <- c("ID", "z")
  DFBIN[[i]]$BF <- i
}

DFBIN <- do.call(rbind, DFBIN)

Accuracy <- full_join(DF, DFBIN)

Values <- seq(from = 0, to = 1, by = 0.01)

for(i in 1:length(Values)){
  Accuracy[,4+i] <- ifelse(Accuracy$Index > Values[i], 1, 0)
  colnames(Accuracy)[4+i] <- paste0("Threshold_", Values[i])
}


#Accuracy <- Accuracy %>% mutate(eval(parse(paste0("Threshold_","0.1"))) = ifelse(Index >= 0.1, 1, 0))

library(caret)

ForGraph <- data.frame(Threshold = Values, Accuracy = NA)
for(i in 1:length(Values)){
  ForGraph$Accuracy[i] <- confusionMatrix(as.numeric(data.matrix(Accuracy[,(4+i)])),reference =  Accuracy$z, positive = "1")$overall[1]
  ForGraph$Kappa[i] <- confusionMatrix(as.numeric(data.matrix(Accuracy[,(4+i)])),reference =  Accuracy$z, positive = "1")$overall[2]
  ForGraph$TSS[i] <- ((confusionMatrix(as.numeric(data.matrix(Accuracy[,(4 + i)])),reference =  Accuracy$z, positive = "1")$byClass[1] + confusionMatrix(as.numeric(data.matrix(Accuracy[,(4 + i)])),reference =  Accuracy$z, positive = "1")$byClass[2]) -1)
}
library(tidyr)

Selected <- ForGraph %>% summarize(Accuracy = sample(Threshold[Accuracy == max(Accuracy)],1), Kappa = Threshold[Kappa == max(Kappa)], TSS = Threshold[TSS == max(TSS)])
Selected <- gather(Selected)

ForGraph <- dplyr::arrange(ForGraph, desc(Accuracy)) %>% gather(key = Parameter, value = value, -Threshold)

library(ggplot2)

ggplot(ForGraph, aes(x = Threshold, y = value)) + geom_line(aes(color = Parameter)) + theme_classic() + geom_vline(aes(xintercept = value, color = key), lty = 2, data = Selected)


Reoptim <- select(Accuracy, ID, BF, z, contains(as.character(Selected$value[2])))

CostFun <- data.frame(ID = unique(Reoptim$ID), Costo = values(Cost)[unique(Reoptim$ID)])

Reoptim <- left_join(Reoptim, CostFun)

library(tidyr)

CurvaCostos <- Reoptim %>% group_by(BF) %>% summarise(Binary = sum(Costo*z), Optimo = sum(Costo*Threshold_0.11)) %>% gather(key = Modelo, value = Costo, -BF)

ggplot(CurvaCostos, aes(x = Costo, y = BF)) + geom_line(aes(color = Modelo)) + theme_classic() #+ geom_hline(yintercept = 4169, lty = 2, col = "red")

library(purrr)

a <- Reoptim %>%split(.$BF) %>% map_chr(~confusionMatrix(.x$Threshold_0.16, reference = .x$z, positive = "1")$overall[1])

AccuracyByBf <- data.frame(Accuracy = a, BF = unique(Reoptim$BF))
AccuracyByBf$Accuracy <- as.numeric(as.character(AccuracyByBf$Accuracy))

ggplot(AccuracyByBf, aes(x= BF, y = Accuracy)) + geom_line() + theme_classic()
###gif


SOLS1 <- list.files(pattern = "Multispp_[0-9].txt")
SOLS2 <-  list.files(pattern = "Multispp_[0-9][0-9].txt")
SOLS <- c(SOLS1, SOLS2)
library(readr)

blaAcc <- list()
for (i in 1:length(SOLS)){
  DF <- read_delim(SOLS[i], " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

  colnames(DF) <- c("ID", "z")

  temp <- univariate[[1]]

  values(temp) <- NA
  values(temp)[DF$ID] <- ifelse(DF$z > ForGraph$Threshold[1], 1, 0)
  blaAcc[[i]] <- temp
}

animation::saveGIF(for(i in 1:length(blaAcc)){plot(blaAcc[[i]], colNA = "black", main = paste("Sp B,", "Final biomass =", i))}, movie.name = "RedondeoAccurate.gif")




##################OptimoContinuo
SOLS <- list.files(pattern = "Optimo_[0-9].txt")

library(readr)


DF <- read_delim("Optimo_1.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

colnames(DF) <- c("ID", "z")

temp <- univariate[[1]]

values(temp) <- NA
values(temp)[DF$ID] <- DF$z
OptimContRast <- temp


SOLS <- list.files(pattern = "Optimo_[0-9].txt")

library(readr)


DF <- read_delim("OptimoBIN_1.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

colnames(DF) <- c("ID", "z")

temp <- univariate[[1]]

values(temp) <- NA
values(temp)[DF$ID] <- DF$z
OptimBinRast <- temp


library(readr)

DF <- read_delim("OptimoBIN_1cost.txt",  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(DF) <- c("bf", "cost")
DF$SP <- c("A", "B")

DF <- DF %>% filter(SP == "A") %>% mutate(N = bf/mass[1])

DF <- filter(DF, cost > 0)
DF$Buy <- "Continuous"

DFOptimCostBin <- DF


library(RtoAmpl)

setwd("/home/derek/Documents/PostdocPablo/AMPL")


library(readr)
library(dplyr)


DF <- read_delim("Optimo_1.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)
colnames(DF) <- c("ID", "Index")
DF$N <- 4169




DFBIN <- read_delim("OptimoBIN_1.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)
colnames(DFBIN) <- c("ID", "z")
DFBIN$N <- 4169



Accuracy <- full_join(DF, DFBIN)

Values <- seq(from = 0, to = 1, by = 0.01)

for(i in 1:length(Values)){
  Accuracy[,4+i] <- ifelse(Accuracy$Index > Values[i], 1, 0)
  colnames(Accuracy)[4+i] <- paste0("Threshold_", Values[i])
}


library(caret)

ForGraph <- data.frame(Threshold = Values, Accuracy = NA)
for(i in 1:length(Values)){
  ForGraph$Accuracy[i] <- confusionMatrix(as.numeric(data.matrix(Accuracy[,(4+i)])),reference =  Accuracy$z, positive = "1")$overall[1]
  ForGraph$Kappa[i] <- confusionMatrix(as.numeric(data.matrix(Accuracy[,(4+i)])),reference =  Accuracy$z, positive = "1")$overall[2]
  ForGraph$TSS[i] <- ((confusionMatrix(as.numeric(data.matrix(Accuracy[,(4 + i)])),reference =  Accuracy$z, positive = "1")$byClass[1] + confusionMatrix(as.numeric(data.matrix(Accuracy[,(4 + i)])),reference =  Accuracy$z, positive = "1")$byClass[2]) -1)
}
library(tidyr)

Selected <- ForGraph %>% summarize(Accuracy = sample(Threshold[Accuracy == max(Accuracy)],1), Kappa = sample(Threshold[Kappa == max(Kappa)],1), TSS = sample(Threshold[TSS == max(TSS)],1))
Selected <- gather(Selected)

ForGraph <- dplyr::arrange(ForGraph, desc(Accuracy)) %>% gather(key = Parameter, value = value, -Threshold)

library(ggplot2)

ggplot(ForGraph, aes(x = Threshold, y = value)) + geom_line(aes(color = Parameter)) + theme_classic() + geom_vline(aes(xintercept = value, color = key), lty = 2, data = Selected)


Reoptim <- select(Accuracy, ID, N, z, contains(as.character(Selected$value[1])))
colnames(Reoptim)[4] <- c("Threshold")


temp <- univariate[[1]]
values(temp) <- NA
values(temp)[Reoptim$ID] <- Reoptim$Threshold
RasterAcc <- temp

Reoptim <- select(Accuracy, ID, N, z, contains(as.character(Selected$value[2])))
colnames(Reoptim)[4] <- c("Threshold")
temp <- univariate[[1]]
values(temp) <- NA
values(temp)[Reoptim$ID] <- Reoptim$Threshold
RasterKappa <- temp

Reoptim <- select(Accuracy, ID, N, z, contains(as.character(Selected$value[3])))
colnames(Reoptim)[4] <- c("Threshold")
temp <- univariate[[1]]
values(temp) <- NA
values(temp)[Reoptim$ID] <- Reoptim$Threshold
RasterTSS <- temp



plot(SolutionsStack, colNA = "black")


OptimBinRast <- as.factor(OptimBinRast)
rat <- levels(OptimBinRast)[[1]]
rat[["decision"]] <- c("Not buy", "buy")
levels(OptimBinRast) <- rat

RasterAcc <- as.factor(RasterAcc)
rat <- levels(RasterAcc)[[1]]
rat[["decision"]] <- c("Not buy", "buy")
levels(RasterAcc) <- rat

RasterKappa <- as.factor(RasterKappa)
rat <- levels(RasterKappa)[[1]]
rat[["decision"]] <- c("Not buy", "buy")
levels(RasterKappa) <- rat

RasterTSS <- as.factor(RasterTSS)
rat <- levels(RasterTSS)[[1]]
rat[["decision"]] <- c("Not buy", "buy")
levels(RasterTSS) <- rat

library(rasterVis)
myTheme <- BTCTheme()
myTheme$panel.background$col = 'black'

SolutionsStack <- stack(OptimBinRast, RasterAcc, RasterKappa, RasterTSS)

names(SolutionsStack) <- c("Binary", "Accuracy", "Kappa", "TSS")

levelplot(SolutionsStack, col.regions=terrain.colors(2), xlab="", ylab="", par.settings = myTheme)

CostFun <- data.frame(ID = unique(Reoptim$ID), Costo = values(Cost)[unique(Reoptim$ID)])

Reoptim <- left_join(Reoptim, CostFun)


names(univariate) <- as.character(paste("Year",seq(2000, 2070, by = 10)))

library(tidyr)
