omcabasic <- function(xo,np , nmod , tmod , rows, idr, idc, idcv,vordered){
     xodisg <- XO <- nc <- mj1 <- mj2 <- list()
     idc2 <- NULL
          for(i in 1:np){
               XO[[i]] <- insertval2(xo[, i], nmod[i])
               XO[[i]] <- matrix(unlist(XO[[i]]),rows,nmod[[i]])
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
     rows <- rismca$nr
     tot <- rismca$tot
     idc <- rismca$idc
     xo <- rismca$xo
     dj <- rismca$dj * np
     dj2 <- rismca$dj
     autovetn <- rismca$autovetn
     autovet <- rismca$autovet
     uni1 <- rismca$Rweights
     mu <- rismca$values #eigenvalues
     Burt <- rismca$Burt
     ####################################################################
     #                                                                  # 
     # Computation of polynomials                                       #
     #                                                                  # 
     ####################################################################
     nburt <- dim(Burt)
     Superpoly <- matrix(0, nburt[1], nburt[2])
     numr <- 1
     k <- 1
     listBpoly <- list()
numc <- nmod[1]
 for(i in 1:np){
          Bpoly <- orthopoly(c(diag(dj))[numr:numc], c(1:nmod[i]))
               listBpoly[[i]]<-Bpoly[,-1]   # Orthogonal polynomial without the trivial component
          dimnames(listBpoly[[i]]) <- list(idc3[[i]], paste("ax", 1:(nmod[i] - 1), sep = ""))
          Bpoly2 <-  Bpoly[,-1]
#browser()
if (vordered[i]==TRUE){
          Superpoly[numr:numc, numr:numc] <- Bpoly}
else {Superpoly[numr:numc, numr:numc] <-autovet[numr:numc, numr:numc]}
#browser()
numr<-numr+nmod[i]
numc<-numc+nmod[i+1]         
}
     Z <- t(autovetn) %*% (xo/sqrt(rows * np)) %*% Superpoly
     tZ<-t(Z)
     Coordi <- autovetn %*% Z
     dimnames(Z) <- list(idr, NULL)
percentage<-list()
nvordered<-0
#if (all(vordered)==T) {nvordered=nvordered+1}
percenta=matrix(0,max(nmod),tmod)
ncluster=list()
j=1
for (i in 1:tmod){   
percentage[i] <- miocount(Coordi[, i]) #percentage of linear compnents of the np variables
#percenta[1:length(percentage[i]),i]<-c(unlist(percentage[i]))
#ncluster[i]<-sort(unique(Coordi[,i])) #number of unique values
j=j+1
#browser()
}
#LinearPerc=matrix(0,max(nmod[vordered]),length(nmod[vordered]))
LinearPerc=list()
cost=sum(nmod[!vordered])+length(nmod[!vordered])
for (i in 1:length(nmod[vordered])){
LinearPerc[i]=percentage[cost]
cost=cost+nmod[vordered][[i]]
}
#percmean=apply(as.matrix(LinearPerc),1,mean)
#browser()
percmean=apply(matrix(unlist(LinearPerc), max(nmod[vordered]),length(nmod[vordered])),1,mean)
    LinearPercentage <- round(percmean/rows * 100, digits = 1)
     idj2 <- solve(sqrt(dj))
 omcabasicresults <- list(RX = Z, CX = tZ, Rweights 
          = uni1, Cweights = idj2, nmod = nmod, tmod = tmod, np = np, 
          Raxes = Superpoly, Caxes = autovetn,  mu = mu, dj = dj, xo 
          = xo, listBpoly = listBpoly, LinearPercentage = 
          LinearPercentage, BURT = aBURT)
}
