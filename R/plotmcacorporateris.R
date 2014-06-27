plotmcacorporateris <-
function(x, catype = "mca", firstaxis = 1, 
     lastaxis = 2, cex = 0.8, prop = 1){

     if ((firstaxis < 1)|(firstaxis>x@maxaxes-1)){ 
          stop(paste("incorrect first axis =",firstaxis, "\n\n")) 
     }
     if (lastaxis > x@maxaxes){
          stop(paste("incorrect last axis = ", lastaxis, "\n\n")) 
     }
     if (firstaxis >= lastaxis){
          stop(paste("last axis must be greater than first axis\n\n"))
     }

     cord2 <- x@Rprinccoord
     cord1 <- x@Cprinccoord
     row.names(cord2) <- x@rowlabels
     row.names(cord1) <- x@collabels
     main <- "Classical Plot"

     graph1(cord1[, c(firstaxis, lastaxis)], cex = cex, a1 = firstaxis, 
          a2 = lastaxis, inertiapc = x@inertiasAdjusted[,2], main = 
          "Category Plot", prop = prop, Inames = row.names(cord1)) 

     if (catype == "mca"){
          {Inames<-row.names(cord2)}
     } else {
          Inames<-paste("C",1:x@br@nmod,sep="")
     }
     dev.new()
     graph1(cord2[, c(firstaxis, lastaxis)], cex = cex, a1 = firstaxis, 
          a2 = lastaxis, inertiapc = x@inertiasAdjusted[,2], main = 
          "Individual Plot", prop = prop, Inames = Inames) 
} 