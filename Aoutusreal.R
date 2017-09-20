library(RtoAmpl)

e <- new("Extent"
         , xmin = -67.8977128959309
         , xmax = -49.4861624007285
         , ymin = -31.813667063965
         , ymax = -10.9555025824394
)

data("Aotus")
data("Costo")

miniaotus <- crop(Aotus, e)
costito <- resample(Costo, miniaotus, method = "ngb")

system.time(a <-RasterToAmplDat(Stack = miniaotus, maxalpha = 5, maxbiomass = 24000, maxcapacidad = 24000, name ="aotusito800", Threshold = 800, costlayer = costito, Dist = 20000))

system.time(a <-RasterToAmplDat(Stack = miniaotus, maxalpha = 5, maxbiomass = 24000, maxcapacidad = 24000, name ="aotusito700", Threshold = 700, costlayer = costito, Dist = 20000))


temp <- miniaotus[[1]]

library(readr)


Miniaotus500 <- read_delim("~/Documents/PostdocPablo/AMPL/Miniaotus500.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

colnames(Miniaotus500) <- c("ID", "Values")

values(temp) <- ifelse(is.na(values(temp)), NA, 0)

values(temp)[Miniaotus500$ID] <- Miniaotus500$Values

plot(miniaotus, colNA = "black")
plot(costito, colNA = "black")
plot(temp, colNA = "black")

plot(miniaotus[[1]], colNA = "black")
plot(temp, colNA = "black", add = TRUE, legend = FALSE, alpha = 0.5)

plot(miniaotus[[8]], colNA = "black")
plot(temp, colNA = "black", add = TRUE, legend = FALSE, alpha = 0.5)


temp2 <- miniaotus[[1]]
values(temp2) <- ifelse(is.na(values(temp2)), NA, 0)
values(temp2)[a$Biomasa$ID] <- a$Biomasa$Biomasa

pa <- costito

values(pa) <- ifelse(values(pa) == 0, 1, NA)

plot(pa)


plot(temp, colNA = "black")
plot(pa, colNA = "black")

sol <- merge(Miniaotus500, a$Cost)

sol <- dplyr::filter(sol, Values != 0)


###################



OneMill <- miniaotus[[1]]

library(readr)


Miniaotus1M <- read_delim("~/Documents/PostdocPablo/AMPL/Miniaotus500.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

colnames(Miniaotus1M) <- c("ID", "Values")

values(OneMill) <- ifelse(is.na(values(OneMill)), NA, 0)

values(OneMill)[Miniaotus1M$ID] <- Miniaotus1M$Values

plot(miniaotus, colNA = "black")
plot(costito, colNA = "black")
plot(OneMill, colNA = "black")

##############

FiveMill <- miniaotus[[1]]

library(readr)


Miniaotus5M <- read_delim("~/Documents/PostdocPablo/AMPL/Miniaotus500.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

colnames(Miniaotus5M) <- c("ID", "Values")

values(FiveMill) <- ifelse(is.na(values(FiveMill)), NA, 0)

values(FiveMill)[Miniaotus5M$ID] <- Miniaotus5M$Values

plot(miniaotus, colNA = "black")
plot(costito, colNA = "black")
plot(FiveMill, colNA = "black")

plot(miniaotus[[1]], colNA = "black")
plot(FiveMill, colNA = "black", add = TRUE, legend = FALSE, alpha = 0.5)

plot(miniaotus[[8]], colNA = "black")
plot(FiveMill, colNA = "black", add = TRUE, legend = FALSE, alpha = 0.5)

hist(Miniaotus5M$Values)


sol <- merge(Miniaotus5M, a$Cost)
sol <- dplyr::filter(sol, Values != 0 & Cost < 5000)
Miniaotus5M

#####


##############

ThreeMill <- miniaotus[[1]]

library(readr)


Miniaotus3M <- read_delim("~/Documents/PostdocPablo/AMPL/Miniaotus500.txt", " ", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE)

colnames(Miniaotus3M) <- c("ID", "Values")

values(ThreeMill) <- ifelse(is.na(values(ThreeMill)), NA, 0)

values(ThreeMill)[Miniaotus3M$ID] <- Miniaotus3M$Values

plot(miniaotus, colNA = "black")
plot(costito, colNA = "black")
plot(ThreeMill, colNA = "black")

plot(miniaotus[[1]], colNA = "black")
plot(ThreeMill, colNA = "black", add = TRUE, legend = FALSE, alpha = 0.5)

plot(miniaotus[[8]], colNA = "black")
plot(ThreeMill, colNA = "black", add = TRUE, legend = FALSE, alpha = 0.5)

hist(Miniaotus3M$Values)


sol <- merge(Miniaotus3M, a$Cost)
sol <- dplyr::filter(sol, Values != 0)
hist(sol$Cost)

#biomasa

Miniaotus3Mb <- read_delim("~/Documents/PostdocPablo/AMPL/Miniaotus500y.txt",  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

colnames(Miniaotus3Mb) <- c("ID", "Time", "Biomass")

times <- unique(Miniaotus3Mb$Time)


StackBiomass <- list()

for (i in 1:length(times)){
  temp <- miniaotus[[1]]
  values(temp) <- ifelse(is.na(values(temp)), NA, 0)
  DF <- dplyr::filter(Miniaotus3Mb, Time == times[i])
  values(temp)[DF$ID] <- DF$Biomass
  StackBiomass[[i]] <- temp
}

StackBiomass <-do.call(stack, StackBiomass)

names(StackBiomass) <- c("2000", "2010", "2020","2030", "2040", "2050", "2060", "2070")

rasterVis::levelplot(StackBiomass)

library(animation)

brks <- seq(0, 24000, length.out = 11)
nb <- length(brks)-1
colors <- rev(terrain.colors(nb))

saveGIF(for(i in 1:nlayers(StackBiomass)){plot(StackBiomass[[i]], main = names(StackBiomass)[i], col = colors, breaks = brks)}, "BiomassAotus.gif")
