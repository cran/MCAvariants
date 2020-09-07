caplot3d<-function(coordR,coordC,inertiaper,firstaxis=1,lastaxis=2,thirdaxis=3){
#------------------------------------
# 3Dim plot 
# coordR= row  coordinates
# coordC= column coordinates  
# inertiaper= percentage of explained inertia per dimension
#------------------------------------  
rc = round(coordR,2)
I<-dim(rc)[2]
cc = round(coordC,2)
iner<-paste(round(inertiaper,1),"%",sep="")
colnames(rc)<-paste("Dim",1:I,sep="")
iner1<-paste("(",iner,")",sep="")
namecol<-paste(colnames(rc),iner1,sep=" ")
colnames(rc)<-namecol
p = plot_ly(type="scatter3d",mode = 'text') 
p = add_trace(p, x = rc[,firstaxis], y = rc[,lastaxis], z = rc[,thirdaxis],
              mode = 'text', text = rownames(rc),
              textfont = list(color = "red"), showlegend = FALSE) 
p = add_trace(p, x = cc[,firstaxis], y = cc[,lastaxis], z = cc[,thirdaxis], 
              mode = "text", text = rownames(cc), 
              textfont = list(color = "blue"), showlegend = FALSE) 
p <- config(p, displayModeBar = FALSE)
p <- layout(p, scene = list(xaxis = list(title = colnames(rc)[1]),
                            yaxis = list(title = colnames(rc)[2]),
                            zaxis = list(title = colnames(rc)[3]),
                            aspectmode = "data"),
            margin = list(l = 0, r = 0, b = 0, t = 0))
p$sizingPolicy$browser$padding <- 0
my.3d.plot = p
my.3d.plot
}
