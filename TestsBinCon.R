library(readr)
library(RtoAmpl)

data("univariate")


costuniv <- univariate
values(costuniv) <- 1

RasterToAmplDat(univariate, Dist = 300000, costlayer = costuniv, name = "TESTUNIV", Threshold = 0.6)

#Conservation goal 68, cost = 12.80918391
TESTUNIV <- read_delim("~/Documents/PostdocPablo/AMPL/TESTUNIV.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(TESTUNIV) <- c("ID", "Value")

TESTUNIVROUND <- TESTUNIV
TESTUNIVROUND$Value <- ceiling(TESTUNIVROUND$Value)

UNI <- univariate[[1]]
values(UNI) <- 0

UNIROUND <- univariate[[1]]
values(UNIROUND) <- 0

values(UNI)[TESTUNIV$ID] <- TESTUNIV$Value

plot(UNI)

values(UNIROUND)[TESTUNIVROUND$ID] <- TESTUNIVROUND$Value

plot(UNIROUND)

#Conservation goal 68, cost = 21
TESTUNIVB <- read_delim("~/Documents/PostdocPablo/AMPL/TESTUNIVB.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(TESTUNIVB) <- c("ID", "Value")

UNIB <- univariate[[1]]

values(UNIB) <- 0

values(UNIB)[TESTUNIVB$ID] <- TESTUNIVB$Value

plot(UNIB)

##################More tests

RasterToAmplDat(univariate, Dist = 300000, costlayer = costuniv, name = "TESTUNIV8", Threshold = 0.6, maxalpha = 8)

#Conservation goal 56, cost = 11.95

TESTUNIV8 <- read_delim("~/Documents/PostdocPablo/AMPL/TESTUNIV8.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(TESTUNIV8) <- c("ID", "Value")

UNI8 <- univariate[[1]]

values(UNI8) <- 0

values(UNI8)[TESTUNIV$ID] <- TESTUNIV$Value

plot(UNI8)

#Conservation goal 56, cost = 22
TESTUNIVB8 <- read_delim("~/Documents/PostdocPablo/AMPL/TESTUNIVB8.txt", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(TESTUNIVB8) <- c("ID", "Value")

UNIB8 <- univariate[[1]]

values(UNIB8) <- 0

values(UNIB8)[TESTUNIVB8$ID] <- TESTUNIVB8$Value

plot(UNIB8)

RasterToAmplDat(univariate, Dist = 300000, costlayer = costuniv, name = "TESTUNIV6", Threshold = 0.6, maxalpha = 6)
RasterToAmplDat(univariate, Dist = 300000, costlayer = costuniv, name = "TESTUNIV4", Threshold = 0.6, maxalpha = 4)
RasterToAmplDat(univariate, Dist = 300000, costlayer = costuniv, name = "TESTUNIV2", Threshold = 0.6, maxalpha = 2)
