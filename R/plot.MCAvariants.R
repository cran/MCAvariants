plot.MCAvariants <-
function(x, catype = "mca", firstaxis = 1, 
     lastaxis = 2, cex = 0.8, cex.lab=0.8,prop=1,M=2,...){

     if ((firstaxis < 1)|(firstaxis>x$maxaxes-1)){ 
          stop(paste("incorrect first axis =",firstaxis, "\n\n")) 
     }
     if (lastaxis > x$maxaxes){
          stop(paste("incorrect last axis = ", lastaxis, "\n\n")) 
     }
     if (firstaxis >= lastaxis){
          stop(paste("last axis must be greater than first axis\n\n"))
     }

     cord1 <- x$Rprinccoord
    cord2 <- x$Cprinccoord
     rowlab <- x$rowlabels
     collab <- x$collabels
nrow <- x$rows
ncol <- x$tmod
#-------------------------------------------------------------------------------------------
#library(scales)
#library (ggplot2)
#library(ggrepel)
#library(gridExtra)
categ<-NULL
frows <- data.frame(coord=cord1, labels=rowlab, categ=rep("rows", nrow)) # build a dataframe to be used as input for plotting via ggplot2
gcols <- data.frame(coord=cord2, labels=collab, categ=rep("cols", ncol)) # build a dataframe to be used as input for plotting via ggplot2
#-------------------------------------------------------------
#FGcord <- rbind(frows, gcols)                                       # build a dataframe to be used as input for plotting via   
FGcord <- gcols                                     # build a dataframe to be used as input for plotting via   

xmin <- min(FGcord[,firstaxis],FGcord[,lastaxis])
 xmax <- max(FGcord[,firstaxis],FGcord[,lastaxis])
ymin <- min(FGcord[,lastaxis],FGcord[,firstaxis])
ymax <- max(FGcord[,lastaxis],FGcord[,firstaxis])
 MCAplot <- ggplot(FGcord, aes(x=FGcord[,firstaxis], y=FGcord[,lastaxis]),type="b") + 
    geom_point(aes(color=categ, shape=categ), size=1.5) +
    geom_vline(xintercept = 0, linetype=2, color="gray") + 
    geom_hline(yintercept = 0, linetype=2, color="gray") + 
    labs(x=paste0("Principal Axis ",firstaxis,sep=" ", round(x$inertiapc[1],), "%"),y=paste0("Principal Axis ",lastaxis,sep= " ", round(x$inertiapc[2],0),"%"))  +  
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    theme(panel.background = element_rect(fill="white", colour="black")) + 
    scale_colour_manual(values=c("blue", "red")) + 
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    geom_text_repel(data=FGcord, aes(colour=categ, label = labels), size = 3) +
#geom_line(aes(group=categ,linetype=categ),lwd=.2)+
#scale_linetype_manual(values=c("rows"="dashed","cols"="blank"))+
#geom_segment(data=frows,aes(x=rep(0,c(nrow)),y=rep(0,c(nrow)),xend=frows[, firstaxis],yend=frows[, 
#lastaxis]),colour=rep("blue",nrow))+
theme(legend.position="none")+
   ggtitle(" ") 
  grid.arrange(MCAplot, ncol=1)
 
#------------------------------------------------------------------------------
#     graph1(cord1[, c(firstaxis, lastaxis)], cex = cex,cex.lab=cex.lab, a1 = firstaxis, 
#         a2 = lastaxis, inertiapc = x$inertiasAdjusted[,2], main = 
#         "Category Plot", prop = prop, Inames = row.names(cord1)) 
#
#    
#     dev.new()
#     graph1(cord2[, c(firstaxis, lastaxis)], cex = cex,cex.lab=cex.lab, a1 = firstaxis, 
#          a2 = lastaxis, inertiapc = x$inertiasAdjusted[,2], main = 
#          "Individual Plot", prop = prop, Inames = Inames) 
} 