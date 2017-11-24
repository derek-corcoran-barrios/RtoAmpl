library(gdistance)
library(raster)

SPA <- raster(nrows=3, ncols=3, xmn = -10, xmx = -4, ymn = 4, ymx = 10)

values(SPA) <- c(0.1, 0.4, 0.6, 0, 0.2, 0.4, 0, 0.1, 0.2)

B <- structure(c(-9, -7, -5, -9, -7, -5, -9, -7, -5, 9, 9, 9, 7, 7,
            7, 5, 5, 5), .Dim = c(9L, 2L), .Dimnames = list(NULL, c("x",
                                                                    "y")))

h16  <- transition(SPA, transitionFunction=function(x){1},16,symm=FALSE)
h16   <- geoCorrection(h16, scl=FALSE)


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


connections <- data.frame(from = rep(1:nrow(B), each = nrow(B)),to = rep(1:nrow(B), nrow(B)), dist =as.vector(apply(B,1, accCost2, x = h16)))


accept


There are several syntax issues in your code.

This code works for me.

library("parallel")

accCost_wrap <- function(x){accCost2(h16,x)}
#Instead of including h16 in the parRapply function,
#just get it in the node environment

cl = makeCluster(3)

clusterExport(cl, c("h16", "accCost2"))
#B will be "sent" to the nodes through the parRapply function.

clusterEvalQ(cl, {library(gdistance)})
#raster is a dependency of gdistance, so no need to include raster here.

pp <- parRapply(cl, x=B, FUN=accCost_wrap)

stopCluster(cl)

connections <- data.frame(from = rep(1:nrow(B), each = nrow(B)),
                          to = rep(1:nrow(B), nrow(B)),
                          dist = as.vector(pp))

