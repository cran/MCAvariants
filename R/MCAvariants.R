MCAvariants <- function(Xtable, catype = "omca",np = 5, vordered=c(TRUE,TRUE,TRUE,TRUE,TRUE)) { 
#--------------------------------------------------------------------------------------------------------
#the ordered variables should be indicated using the flag parameter vordered and should be close each other
#happymca<-tableconvert(happy)
#dimnames(happymca)=list(paste("i",1:1517,sep=""),c("happy","sib","edu"))
#library(R.utils)
#sourceDirectory("C:\\Users\\Amministratore\\Dropbox\\Rnsca\\MCAvariants\\R")
#-------------------------------------------------------------------------------------------------------
     X <- as.matrix(Xtable)
          rowlabels <- rownames(Xtable)
          collabels <- colnames(Xtable) 
     idcv <- list()
     rows <- dim(X)[1]
     nmod <- vector()
     comp <- vector()
     comppvalue1 <- vector()
     degreef <- numeric()
     for (i in 1:np){
          nmod[i] <- max(X[,i])
          idcv[[i]] <- paste(collabels[i], 1:nmod[i], sep = "")
     }

     tmod <- sum(nmod)
     n <- sum(X)
     cols <- tmod
     xw <- rep(1,rows) #weight vector
     maxaxes <- min(rows - 1, cols - 1)
     S <- switch(catype, "mca" = mcabasic(X, np = np, nmod = nmod, tmod 
          = tmod, rows = rows, idr = rowlabels, idc = collabels, idcv = 
          idcv), "omca" = omcabasic(X,np = np, nmod = nmod, tmod 
          = tmod, rows = rows, idr = rowlabels, idc = collabels, idcv = 
          idcv,vordered=vordered))

     ####################################################################
     #                                                                  #
     # (Classic) Multiple Correspondence Analysis                       #
     #                                                                  #
     ####################################################################
    
     if(catype == "mca"){
          Fmat2 <- ((S$xo/np)/sqrt(rows))  %*% S$Raxes 
                              # To reconstruct the total inertia
          Fmat<- Fmat2[,-1]   # Row principal coordinates
          Gmat <- solve(S$dj) %*% t((S$xo/np)/sqrt(rows)) %*%  S$Caxes[,1:tmod] %*% diag(1/S$mu)
          Gmat<-Gmat[,-1]     # Column principal coordinates
          degreef <- (1/2*((tmod - np)^2 - np*(nmod[1] - 1)^2))
          comp <- rows*(S$RX[-1,-1])^2
          comppvalue1 <- 1 - pchisq(sum(S$mu^2), degreef) 
                              # p-value  of chi-square inertia
     }
     ####################################################################
     #                                                                  #
     # Ordered Multiple Correspondence Analysis                         #
     #                                                                  #
     ####################################################################

     if((catype == "omca")||(catype == "OMCA")){
          Fmat <- S$Caxes %*%  S$RX
          Gmat <- S$Raxes %*% S$CX      # Variable coordinates
          Gmat <- Gmat[,-1]
          if ((S$RX[2,1] < 0) & (S$RX[2,2] > 0)|(S$RX[2,1] < 0) & 
                                                  (S$RX[2,2] > 0)){
               Fmat[,3] < -(-1) * Fmat[,3]
               Fmat[,2] < -(-1) * Fmat[,2]
          } 
          degreef <- 1/2 * ((tmod - np)^2 - np * (nmod[1] - 1)^2)
          comp <- diag(t(S$RX) %*% (S$RX) %*% t(S$RX) %*% 
               (S$RX))[1:(tmod-np)]     # The total inertia
          comp <- comp*rows
          comppvalue1 <- 1 - pchisq(comp, degreef/(tmod-np)) 
          comppvalue1 <- round(as.matrix(comppvalue1),digits=4)
     }

     ####################################################################
     #                                                                  # 
     # OTHER CALCULATIONS                                               #
     #                                                                  #
     ####################################################################
       
     inertia <- (S$mu[-1])^4 
     m <- 0
     benz <- S$mu[-1]^2

     for (i in 1:(tmod - 1)){
          if (benz[i] >= 1/np){
               m<-m+1
          } else { 
               m <- m 
          }
     }

     inertiaB <- round(((np/(np - 1))^2*(S$mu[-1]^2 - 1/np)^2), digits 
          = 3)
     inertiaBsum <- sum(((np/(np - 1))^2*(S$mu[-1]^2 - 1/np)^2)[1:m])           
                              # Benzecri's adjusted value   
     inertiapc <- (((np/(np - 1))^2*(S$mu[-1]^2 - 1/np)^2)/inertiaBsum) * 100
     inertiaBurt <- round(inertia, digits = 4) 
     inertiaBurtsum <- sum(inertia) 
     inertiaX <- round(S$mu[-1]^2, digits = 4)
     inertiaXsum <- sum(S$mu[-1]^2)
     cuminertiapc <- cumsum(inertiapc)
     inertiapc <- round(inertiapc,digits=1)
     cuminertiapc <- round(cuminertiapc,digits=1)
     inertias <- cbind(inertiaX, inertiaBurt)
     inertiasAdjusted <- cbind(inertiaB[1:m], inertiapc[1:m], 
          cuminertiapc[1:m])
     Xstd <- S$xo/sum(S$xo)
     rownames(X) <- rowlabels
     colnames(X) <- collabels
     collabels2 <- dimnames(S$xo)[[2]]
#browser()
#     mcacorporateris<- new("mcacorporaterisclass", br = S, Xtable = S$xo, 
 #         rows = rows, cols = cols,  rowlabels = rowlabels, collabels = 
  #        collabels2, Rprinccoord = round(Fmat, digits = 4), Cprinccoord 
   #       = round(Gmat, digits = 4), inertiaXsum = inertiaXsum, 
    #      inertiaBurtsum = inertiaBurtsum, inertias = inertias, 
     #     inertiasAdjusted = inertiasAdjusted, catype = catype,  maxaxes = maxaxes, comp = comp, componentpvalue1 = 
      #    comppvalue1, degreef = degreef)
#-----------------------------------------------------------------------------------

resultMCA<-list( Xtable = S$xo, BURT=S$BURT,listBpoly=S$listBpoly,LinearPercentage=S$LinearPercentage, 
Rweights=S$Rweights, Cweights=S$Cweights,Raxes=S$Raxes,Caxes=S$Caxes,np=np,
rows = rows, cols = cols, nmod=nmod,tmod=tmod, np=np,rowlabels = rowlabels, collabels = 
          collabels2, Rprinccoord = Fmat, Cprinccoord 
          = Gmat, inertiaXsum = inertiaXsum, 
          inertiaBurtsum = inertiaBurtsum, inertias = inertias, inertiapc=inertiapc,
          inertiasAdjusted = inertiasAdjusted, catype = catype,  maxaxes = maxaxes, comp = comp, componentpvalue1 = 
          comppvalue1, degreef = degreef)
class(resultMCA)<-"MCAvariants"
return(resultMCA)

}
