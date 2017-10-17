multisppBiomass <- function(Stacklist, maxalpha = 5, maxbiomass = 2, maxcapacidad = 10, name = "Stack", Dist = 500000, Threshold, costlayer){

  Masklayer <- costlayer
  values(Masklayer) <- ifelse(is.na(values(Masklayer)), NA, 1)
  for (i in 1:length(Stacklist)){
    Stacklist[[i]] <- Stacklist[[i]] * Masklayer
    TempStack <- Stacklist[[i]]

    values(TempStack)[values(TempStack) < Threshold] = 0
    values(TempStack)[values(TempStack) >= Threshold] = 1

    TempRaster <- sum(TempStack)

    TempRaster[values(TempRaster) > 0] = 1
    TempRaster[values(TempRaster) == 0] = NA

    Stacklist[[i]] <- Stacklist[[i]]*TempRaster
  }


  Alphas <- list()
  for(j in 1:length(Stacklist)){
    Alpha <- list()
    for (i in 1:nlayers(Stacklist[[j]])){
      temp <- data.frame(Alpha = values(Stacklist[[j]][[i]]), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i)
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
      temp <- data.frame(Capacidad = values(Stacklist[[j]][[i]]), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i)
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

  cost <-values(costlayer)[Nodos]

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
  cat(c("set c :=", round(cost,4), ";"))
  cat("\n")
  sink()
  return(list(Nodos = Nodos, Biomasa = Biomasa, Alphas = Alphas, Alpha = Alpha, Cost = Cost))
}










##########################################################
  #########################################################
Stacklist <- list(univariate, bivariate)
names(Stacklist) <- c("spA", "spB")
maxbiomass <- 2
maxcapacidad <- 10
costlayer <- Cost
Dist = 500000
  Masklayer <- Cost
  values(Masklayer) <- ifelse(is.na(values(Masklayer)), NA, 1)
  for (i in 1:length(Stacklist)){
    Stacklist[[i]] <- Stacklist[[i]] * Masklayer
    TempStack <- Stacklist[[i]]

    values(TempStack)[values(TempStack) < 0.3] = 0
    values(TempStack)[values(TempStack) >= 0.3] = 1

    TempRaster <- sum(TempStack)

    TempRaster[values(TempRaster) > 0] = 1
    TempRaster[values(TempRaster) == 0] = NA

    Stacklist[[i]] <- Stacklist[[i]]*TempRaster
  }


    Alphas <- list()
    for(j in 1:length(Stacklist)){
      Alpha <- list()
      for (i in 1:nlayers(Stacklist[[j]])){
        temp <- data.frame(Alpha = values(Stacklist[[j]][[i]]), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i)
        Alpha[[i]] <- temp[complete.cases(temp),]
      }
      Alphas[[j]]<- do.call("rbind", Alpha)
      Alphas[[j]]$Spp <- names(Stacklist)[j]
    }

    Alpha <-  do.call(rbind, Alphas)
    Alpha$Alpha <- Alpha$Alpha/max(Alpha$Alpha)*5
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
    ###alpha LISTO
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

#Biomasa listo
    Capacidades <- list()
    for(j in 1:length(Stacklist)){
      Capacidad <- list()
      for (i in 1:nlayers(Stacklist[[j]])){
        temp <- data.frame(Capacidad = values(Stacklist[[j]][[i]]), ID = 1:length(values(Stacklist[[j]][[i]])), Time = i)
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

    cost <-values(costlayer)[Nodos]
    ###Costo y capacidad listos
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
        setNames(data.frame(paste("[",temp[[i]][["Sp"]],temp[[i]][["from"]], ",", temp[[i]][["to"]], ",*","]", sep = ""), temp[[i]]["Time"], temp[[i]]["beta"]),
                 c("V", paste("T", i, sep = ""), i-1))
      } else if (i == length(temp)){
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["beta"], rep("\n", NROW(temp[[i]]))),
                 c(paste("T", i, sep = ""), i-1, "line"))
      } else {
        setNames(data.frame(temp[[i]]["Time"], temp[[i]]["beta"]),
                 c(paste("T", i, sep = ""), i-1))
      }
    }))


multisppBiomass(Stacklist = Stacklist, costlayer = Cost, Threshold = 0.2)