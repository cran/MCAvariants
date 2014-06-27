mcabasic <- function(xo, np = 2, nmod = c(5, 5), tmod = 10, rows = 511, 
     idr, idc, idcv){

     xodisg <- XO <- nc <- mj1 <- mj2 <- list()
     nr <- nrow(xo)
     np <- np
     idc2 <- NULL

          for (i in 1:ncol(xo)){
               XO[[i]] <- insertval2(xo[, i], nmod[i])
               XO[[i]] <-matrix(unlist(XO[[i]]), nr,nmod[[i]])
               dimnames(XO[[i]]) <- list(idr, paste(idc[[i]], 
                    1:nmod[[i]], sep = ""))
               idc2 <- c(idc2, dimnames(XO[[i]])[2])
          }
          idc2 <- unlist(idc2)
      
     cBURT <- aBURT<-list()
     k <- 1

     for(j in 1:np) {
          for(i in 1:np) {
               cBURT[[k]] <- t(XO[[j]]) %*% XO[[i]]
               aBURT[[k]] <- cBURT[[k]] # Tables of absolute frequencies
               cBURT[[k]] <- cBURT[[k]]/sum(cBURT[[k]])
               k <- k + 1
          }
     }

     BURT1 <- list()
     k <- 1
     num1 <- 2
     num2 <- np
     num3 <- 1

     for(j in 1:np){
          for(i in num1:num2){
               BURT1[[j]] <- cbind(cBURT[[num3]], cBURT[[i]])
               cBURT[[num3]] <- BURT1[[j]]
          }
          k <- k + 1
          num1 <- np + num1
          num2 <- k * np
          num3 <- num3 + np
     }

     for(i in 2:np){
          Burt <- rbind(BURT1[[1]], BURT1[[i]])
          BURT1[[1]] <- Burt
     }

     cost <- XO[[1]]
     tot <- sum(XO[[1]])

     for(i in 2:np){
          xo <- cbind(cost, XO[[i]])
          cost <- xo
     }

     nr <- nrow(xo)
     nc <- ncol(xo)
     uni <- matrix(1, nr, 1)
     uni1 <- diag(rep(1,nr))
     uni2 <- diag(rep(1,nc))
     tot <- nr * np
     xc <- t(xo/(tot)) %*% uni     # Column marginals
     xc <- matrix(xc, nc, 1)
     dj <- diag(c(xc))
     dj2 <-diag(c(xc^-1/2))
     di <- diag(c((uni * 1)/nr))
     di2 <- diag(c((uni * 1)/sqrt(nr)))
     Burt <- t((1/np * xo)/sqrt(nr)) %*% (1/np * xo)/sqrt(nr)
     gdj <- solve(dj)
     pcZN <- (xo/(np * sqrt(nr))) %*% gdj %*% t(xo/(np * sqrt(nr)))
     rispcZN <- eigen(pcZN)
     autovetn <- rispcZN$vectors
     valuesn <- rispcZN$values
     pc1 <- ((1/np * xo)/sqrt(nr)) %*% solve(sqrt(dj))
     pc2 <- solve(sqrt(dj)) %*% t((1/np * xo)/sqrt(nr))
     pc0 <- solve(sqrt(dj)) %*% (Burt) %*% solve(sqrt(dj))
     sing <- svd(pc1)$d
     rispc <- eigen(pc0)
     autovet <- solve(sqrt(dj)) %*% rispc$vectors 
     values <- rispc$values
     dimnames(pc2) <- list(idc2, NULL)
     dimnames(autovet) <- list(idc2, NULL)
     idj2 <- solve(sqrt(dj))
     mcabasic <- new("mcabasicresultsclass", RX = pc1, CX = pc2, Rweights = 
          uni1, Cweights = idj2, nmod = nmod, tmod = tmod, np = np, Raxes 
          = autovet, Caxes = autovetn,  mu = sing, dj = dj, xo = xo, BURT  
          = aBURT)
}
