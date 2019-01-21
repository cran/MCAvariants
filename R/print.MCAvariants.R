print.MCAvariants <-
 function(x, printdims=2,...) {
     d <- min(printdims, x$maxaxes)
     axnames <- character(length = d)

     for (i in 1:d) { 
          axnames[i] <- paste("Axis",i) 
     } 

     cat("\n    RESULTS for Variants of Multiple Correspondence    
          Analysis:\n",x$catype, "\n")
     cat("\n BURT table \n")
     print(x$BURT)
     cat("\n Rows in principal coordinates: the first 10 \n")
     printwithaxes(data.frame(x$Rprinccoord[ 1:10,1:d],  row.names=x$rowlabels[1:10]), axnames)
     cat("\n Columns in principal coordinates \n")
#browser()
     printwithaxes(data.frame(x$Cprinccoord[ ,1:d], row.names=x$collabels), axnames)
          if (x$catype == "omca"){
          cat("\n Polynomial functions of each variable \n")
          print(x$listBpoly)
          clusterlabels <- paste("C",1:x$nmod[1],sep="")
          cat("\n Linear Percentage of Clusters \n")
          print(x$LinearPercentage)
          cat("\n  Polynomial Components  of Total Inertia \n")
          print(x$comp)
          cat("\n p-values of Polynomial Components of Total Inertia \n")
          print(x$componentpvalue1[1:(x$tmod - x$np)])
          cat("Degree of Freedom of Polynomial Component", x$degreef
               /(x$tmod - x$np), "\n")
     }
     cat("\n Inertia values of super-indicator and Burt table\n")
     print(round(x$inertias,digits=3))
     cat("\n Benzecri's Inertia values, percentages and cumulative 
          percentages \n")
     print(round(x$inertiasAdjusted, digits = 3))
     cat("Total Degree of Freedom", x$degreef, "\n")
     cat("Total inertia of Super-Indicator table\n")
     print(x$inertiaXsum)
     cat("Total inertia of BURT table\n")
     print(sum(x$inertiaBurtsum))
     cat("Chi-square values of  BURT Inertia \n")
     print(round(x$inertias[,2]*x$rows),dig=3)
     cat("Chi-square value of Total inertia of BURT\n")
     print(round(sum(x$inertiaBurtsum)*x$rows),dig=3) 
}