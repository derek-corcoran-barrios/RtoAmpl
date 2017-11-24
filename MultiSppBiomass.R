multisppBiomass <- function(Stacklist, maxalpha = 5, maxbiomass = 2, maxcapacidad = 10, name = "Stack", Dist = 500000, Threshold, costlayer, bf = c(1,2)){

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
  }

  Alpha <-  do.call(rbind, Alphas)
  Alpha$Alpha <- Alpha$Alpha/max(Alpha$Alpha)*maxalpha
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
    Biomasa[[i]] <- data.frame(Sp = rep(Spps[i], times = length(values(Stacklist[[i]][[1]]))), ID = 1:length(values(Stacklist[[i]][[1]])), Biomasa = (values(Stacklist[[i]][[1]])/max(values(Stacklist[[i]][[1]]), na.rm = T))*maxbiomass)
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
      temp <- data.frame(Capacidad = values(Stacklist[[j]][[i]]), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i-1)
      Capacidad[[i]] <- temp[complete.cases(temp),]
    }
    Capacidades[[j]]<- do.call("rbind", Capacidad)
    Capacidades[[j]]$Spp <- names(Stacklist)[j]
  }

  Capacidad <-  do.call(rbind, Capacidades)
  Capacidad$Capacidad <- Capacidad$Capacidad/max(Capacidad$Capacidad)*maxcapacidad
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

  BF <- data.frame(Spp = Spps, BF = bf, Space = "\n")


  Beta <- list()
  for(i in 1:length(Stacklist)){
    Beta[[i]] <- DistConect(Stacklist[[i]][[1]], Distance = Dist, Time = nlayers(Stacklist[[i]]))
    Beta[[i]] <- dplyr::filter(Beta[[i]], from %in% Nodos, to %in% Nodos)
    Beta[[i]]$Sp <- Spps[i]
  }

  Beta <- do.call(rbind, Beta)

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










##########################################################
#########################################################
library(RtoAmpl)
data(Cost)
data("univariate")
data("bivariate")

Stacklist <- list(univariate, bivariate)
names(Stacklist) <- c("SPA", "SPB")

setwd("/home/derek/Documents/PostdocPablo/AMPL")

library(beepr)
for (i in 1:5){
  multisppBiomass(Stacklist = Stacklist, costlayer = Cost, Threshold = 0.4, name = paste0("Multispp","_", i), bf = c((2*i),i), maxalpha = 4)
  print(i)
  beep(i)
}
beep(8)

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

DF <- list()
for(i in 1:length(Costs)){
  DF[[i]] <- read_delim(Costs[i],  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(DF[[i]]) <- c("bf", "cost")
  DF[[i]]$SP <- c("A", "B")
}

DF <- do.call(rbind, DF)

library(ggplot2)
library(dplyr)
DF <- filter(DF, cost > 0)
DF$Buy <- "Continuous"
ggplot(DF, aes(x = cost, y = bf)) + geom_line(aes(color = SP)) + geom_point(aes(color = SP)) + theme_classic()


animation::saveGIF(for(i in 1:length(bla)){plot(bla[[i]], colNA = "black", main = paste("Sp B,", "Final biomass =", i))}, movie.name = "Continuous.gif")

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

library(ggplot2)
library(dplyr)
DFBIN <- filter(DFBIN, cost > 0)
DFBIN$Buy <- "Binary"

ggplot(DFBIN, aes(x = cost, y = bf)) + geom_line(aes(color = SP)) + geom_point(aes(color = SP)) + theme_classic()


DF <- rbind(DF, DFBIN)
DF <- dplyr::filter(DF, SP == "A")
ggplot(DF, aes(x = cost, y = bf)) + geom_line(aes(color = Buy)) + geom_smooth(aes(color=Buy)) + theme_classic()


animation::saveGIF(for(i in 1:length(blaBIN)){plot(blaBIN[[i]], colNA = "black", main = paste("Sp B,", "Final biomass =", i))}, movie.name = "Binary.gif")
years <- seq(2000, 2070, by = 10)
#animation::saveGIF(for(i in 1:nlayers(univariate)){plot(univariate[[i]], colNA = "black", main = paste("Sp A,", "year", years[i]))}, movie.name = "SpeciesA.gif")
#animation::saveGIF(for(i in 1:nlayers(bivariate)){plot(bivariate[[i]], colNA = "black", main = paste("Sp B,", "year", years[i]))}, movie.name = "SpeciesB.gif")


###############################################################################
##################Redondeo Normal##############################################
###############################################################################


setwd("/home/derek/Documents/PostdocPablo/AMPL")

SOLS1 <- list.files(pattern = "Multispp_[0-9].txt")
SOLS2 <-  list.files(pattern = "Multispp_[0-9][0-9].txt")
SOLS <- c(SOLS1, SOLS2)
library(readr)

blaRound <- list()
for (i in 1:length(SOLS)){
  DF <- read_delim(SOLS[i], " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

  colnames(DF) <- c("ID", "z")

  temp <- univariate[[1]]

  values(temp) <- NA
  values(temp)[DF$ID] <- round(DF$z)
  blaRound[[i]] <- temp
}

animation::saveGIF(for(i in 1:length(blaRound)){plot(blaRound[[i]], colNA = "black", main = paste("Sp B,", "Final biomass =", i))}, movie.name = "RedondeoNormal.gif")


Costs1 <-  list.files(pattern = "Multispp_[0-9]cost.txt")
Costs2 <-  list.files(pattern = "Multispp_[0-9][0-9]cost.txt")
Costs <- c(Costs1, Costs2)

DF <- list()
for(i in 1:length(Costs)){
  DF[[i]] <- read_delim(Costs[i],  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(DF[[i]]) <- c("bf", "cost")
  DF[[i]]$SP <- c("A", "B")
}

DF <- do.call(rbind, DF)

library(ggplot2)
library(dplyr)
DF <- filter(DF, cost > 0)
DF$Buy <- "Continuous"
ggplot(DF, aes(x = cost, y = bf)) + geom_line(aes(color = SP)) + geom_point(aes(color = SP)) + theme_classic()


###############################################################################
##################Redondeo Techo##############################################
###############################################################################


setwd("/home/derek/Documents/PostdocPablo/AMPL")

SOLS1 <- list.files(pattern = "Multispp_[0-9].txt")
SOLS2 <-  list.files(pattern = "Multispp_[0-9][0-9].txt")
SOLS <- c(SOLS1, SOLS2)
library(readr)

blaRoundC <- list()
for (i in 1:length(SOLS)){
  DF <- read_delim(SOLS[i], " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

  colnames(DF) <- c("ID", "z")

  temp <- univariate[[1]]

  values(temp) <- NA
  values(temp)[DF$ID] <- ceiling(DF$z)
  blaRoundC[[i]] <- temp
}

animation::saveGIF(for(i in 1:length(blaRoundC)){plot(blaRoundC[[i]], colNA = "black", main = paste("Sp B,", "Final biomass =", i))}, movie.name = "RedondeoTecho.gif")




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
ForGraph <- dplyr::arrange(ForGraph, desc(Accuracy)) %>% gather(key = Parameter, value = value, -Threshold)

library(ggplot2)

ggplot(ForGraph, aes(x = Threshold, y = value)) + geom_line(aes(color = Parameter)) + theme_classic() + geom_vline(xintercept = 0.5, lty = 2, color = "red")

ROCCURVE <- data.frame(Sensitivity = rep(NA, 101), Specificity = NA)

for(i in 1:length(Values)){
  ROCCURVE$Sensitivity[i] <- confusionMatrix(as.numeric(data.matrix(Accuracy[,(4 + i)])),reference =  Accuracy$z, positive = "1")$byClass[1]
  print(i)
  ROCCURVE$Specificity[i] <-confusionMatrix(as.numeric(data.matrix(Accuracy[,(4 + i)])),reference =  Accuracy$z, positive = "1")$byClass[2]
}

ROCCURVE$FalsePositive <- 1 - ROCCURVE$Specificity

ggplot(ROCCURVE, aes(x = FalsePositive, y = Sensitivity)) + geom_line() + geom_abline(slope = 1, lty =2, color = "red") + theme_classic() + xlim(c(0,1))

Reoptim <- select(Accuracy, ID, BF, z, contains(as.character(ForGraph$Threshold[1])))

CostFun <- data.frame(ID = unique(Reoptim$ID), Costo = values(Cost)[unique(Reoptim$ID)])

Reoptim <- left_join(Reoptim, CostFun)

library(tidyr)

CurvaCostos <- Reoptim %>% group_by(BF) %>% summarise(Binary = sum(Costo*z), Optimo = sum(Costo*Threshold_0.16)) %>% gather(key = Modelo, value = Costo, -BF)

ggplot(CurvaCostos, aes(x = Costo, y = BF)) + geom_line(aes(color = Modelo)) + theme_classic()

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

