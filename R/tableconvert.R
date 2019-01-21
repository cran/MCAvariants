tableconvert <- function(N){

ncats = length(dim(N))
I = dim(N)[1]
J = dim(N)[2]

if (ncats == 3) K = dim(N)[3]

data = NULL

if (ncats == 3){

     for (i in 1:I){
          for (j in 1:J){
               for (k in 1:K){
                    if (N[i,j,k] > 0){
                         data = rbind(data, matrix(rep(c(i,j,k), times = N[i,j,k]), nrow = N[i,j,k], byrow = T))
                    }
               }
          }
     }
}

if (ncats == 2){

     for (i in 1:I){
          for (j in 1:J){
               if (N[i,j] > 0){
                    data = rbind(data, matrix(rep(c(i,j), times = N[i,j]), nrow = N[i,j], byrow = T))
               }
          }
     }
}

data

}