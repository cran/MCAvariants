graph1 <-
function(f, prop=1,cex = 1, cex.lab = 0.8, 
cols = c(2, 4), main="Plot",mar = c(5, 4, 4, 2)+0.1, 
a1,a2,inertiapc,Inames){
#########################################################################
# f and g are the row and column coordinates                            #
#########################################################################
par(pty = "s", mar = mar)
plot(0, 0, pch = " ", xlim = range(f)/prop, ylim = range(f)/prop,
xlab=paste("Axis", a1, " ", inertiapc[a1], "%", sep=""), 
ylab=paste("Axis", a2, " ", inertiapc[a2], "%", sep=""), cex = cex, 
cex.lab = cex.lab,main=main)
abline(h=0,v=0)
nv <- rep(0, nrow(f))
vec <- f[, c(1, 2)]
#if(arrow) {
#arrows(nv, nv, vec[, 1], vec[, 2], length = length)
#	}
text(f[, 1], f[, 2],  label=Inames,cex = cex,cex.lab=cex.lab,pos=1, 
col =cols[1])
}
