corrValue <- function(counts, location){
  require(spdep)
  coords = as.matrix(location)
  k1 = knn2nb(knearneigh(coords))
  all.linked <- max(unlist(nbdists(k1, coords)))
  col.nb <- dnearneigh(coords, 0, all.linked, 
                       row.names=row.names(coords), 
                       longlat = FALSE)
  col.S <- nb2listw(col.nb, style = "S")
  MI.S <- apply(counts, 1, moran, listw = col.S, n = length(col.nb), S0 = Szero(col.S))
  f <- function(x){x = x[[1]]}
  MI.S <- lapply(MI.S, f)
  MI <- NULL
  for(i in 1:length(MI.S)){
    MI <- c(MI, MI.S[[i]])
  }
  names(MI) <- row.names(counts)
  
  return(MI)
}