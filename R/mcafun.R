mcafun <- function(XO, Burt, np, idr, idc, nmod){
   
     cost <- XO[[1]]
     tot <- sum(XO[[1]])

     for(i in 2:np) {
          xo <- cbind(cost, XO[[i]])
          cost <- xo
     }

     nr <- nrow(xo)
     nc <- ncol(xo)
     uni <- matrix(1, nr, 1)
     uni1 <- diag(rep(1, nr))
     uni2 <- diag(rep(1, nc))
     tot <- nr * np
     xc <- t(xo/(tot)) %*% uni 
     xc <- matrix(xc, nc, 1)
     dj <- diag(c(xc))
     di <- diag(c((uni * 1)/nr))
     Burt <- t((1/np * xo)/sqrt(nr)) %*% (1/np * xo)/sqrt(nr)
     gdj <- solve(dj)

     ####################################################################
     #                                                                  # 
     # Using indicator matrices                                         #
     #                                                                  # 
     ####################################################################

     pcZN <- (xo/(np * sqrt(nr))) %*% gdj %*% t(xo/(np * sqrt(nr)))
     rispcZN <- eigen(pcZN)#
     autovetn <- rispcZN$vectors
     valuesn <- rispcZN$values
     pc1 <- ((1/np * xo)/sqrt(nr)) %*% solve(sqrt(dj))
     pc2 <- solve(sqrt(dj)) %*% t((1/np * xo)/sqrt(nr))
     pc0 <- solve(sqrt(dj)) %*% (Burt) %*% solve(sqrt(dj))
     sing <- svd(pc1)$d#
     rispc <- eigen(pc0)
     autovet <- solve(sqrt(dj)) %*% rispc$vectors
     values <- rispc$values
     cordr <- ((1/np * xo)/sqrt(nr)) %*% solve(sqrt(dj)) %*% autovet 
     cordc <- solve(sqrt(dj)) %*% t((1/np * xo)/sqrt(nr)) %*% autovetn
     dimnames(cordr) <- list(idr, NULL)
     dimnames(cordc) <- list(idc, NULL)
     dimnames(autovet) <- list(idc, NULL)
     totin <- sum(values)
     list(xo = xo, xc = xc, autovet = autovet, autovetn = autovetn, 
          values = sing, valuesn = valuesn, pc1 = pc1, pc0 = pc0, pc2 = 
          pc2, dj = dj, totin = totin, tot = tot, sing = sing, nr = nr, 
          Burt = Burt, Raxes = autovet, Caxes = autovetn[,-1], mu = 
          sing[-1], valuesn = valuesn, R = pc1, C = pc2, Rweights = uni1, 
          Cweights = uni2)
}
