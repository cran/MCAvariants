omcabasic <- function(xo){
     xodisg <- XO <- nc <- mj1 <- mj2 <- list()
     nr <- nrow(xo)
     np <- ncol(xo)
     nmod <- vector()

     for (i in 1:np){
          nmod[i] <- max(xo[,i])
     }

     tmod <- sum(nmod)
     idr <- dimnames(xo)[[1]]
     idc <- dimnames(xo)[[2]]
     idcv <- dimnames(xo)[[2]]
     idc2 <- NULL

          for(i in 1:np){
               XO[[i]] <- insertval2(xo[, i], nmod[i])
               XO[[i]] <- matrix(unlist(XO[[i]]),nr,nmod[[i]])
               dimnames(XO[[i]]) <- list(idr, paste(idc[[i]], 
                    1:nmod[[i]], sep = ""))
               idc2 <- c(idc2, dimnames(XO[[i]])[2])
          }
          idc3 <- idc2
          idc2 <- unlist(idc2)
      

     cBURT <- aBURT <- list()
     k <- 1

     for(j in 1:np){
          for(i in 1:np){
               cBURT[[k]] <- t(XO[[j]]) %*% XO[[i]]
               aBURT[[k]] <- cBURT[[k]]
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
          for(i in num1:num2) {
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

     nmod1 <- nmod[1]
     rismca <- mcafun(XO, Burt, np = np, idr = idr, idc = idc2, nmod = 
          nmod1)
     xc <- rismca$xc
     nr <- rismca$nr
     tot <- rismca$tot
     idc <- rismca$idc
     xo <- rismca$xo
     dj <- rismca$dj * np
     dj2 <- rismca$dj
     autovetn <- rismca$autovetn
     autovet <- rismca$autovet
     uni1 <- rismca$Rweights
     values <- rismca$values#eigenvalues
     Burt <- rismca$Burt

     ####################################################################
     #                                                                  # 
     # Computation of polynomials                                       #
     #                                                                  # 
     ####################################################################

     nburt <- dim(Burt)
     Superpoly <- matrix(0, nburt[1], nburt[2])
     num1 <- 1
     num11 <-1
     k <- 1
     listBpoly <- list()
     num2 <- nmod[1]
     num22 <- (nmod[1]-1)

     for(i in 1:np){
          Bpoly <- orthopoly(c(diag(dj))[num1:num2], c(1:nmod[i]))
               listBpoly[[i]]<-Bpoly[,-1] 
               # Orthogonal polynomial without the trivial component
          dimnames(listBpoly[[i]]) <- list(idc3[[i]], paste("ax", 
               1:(nmod[i] - 1), sep = ""))
          Bpoly2 <-  Bpoly[,-1]
          Superpoly[num1:num2, num11:num22] <- Bpoly2
          num1 <- sum(nmod[1:i]) + 1
          num11 <- num11 + (nmod[i] - 1)
          num2 <- (num1 - 1) + nmod[i]
          num22 <- (num11 - 1) + (nmod[i] - 1)
     }

     Z <- t(autovetn) %*% (xo/sqrt(nr * np)) %*% Superpoly
     tZ<-t(Z)
     Coordi <- autovetn %*% Z
     dimnames(Z) <- list(idr, NULL)
     percentage <- miocount(Coordi[, 1])
     LinearPercentage <- round(percentage$conta/nr * 100, digits = 1)
     idj2 <- solve(sqrt(dj))
     omcabasicresults <- new("mcabasicresultsclass", RX = Z, CX = tZ, Rweights 
          = uni1, Cweights = idj2, nmod = nmod[1], tmod = tmod, np = np, 
          Raxes = Superpoly, Caxes = autovetn,  mu = values, dj = dj, xo 
          = xo, listBpoly = listBpoly, LinearPercentage = 
          LinearPercentage, BURT = aBURT)
}
