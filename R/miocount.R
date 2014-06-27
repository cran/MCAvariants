miocount <-
function (x){

     options(digits = 7)
     x <- matrix(sort(round(x, digits = 3)), length(x), 1)
     y <- unique(sort(x))
     ncat <- length(y)
     conta <- matrix(0, ncat, 1)
     i <- 1

     for (j in 1:ncat){
          for (i in 1:nrow(x)){
               if (x[i, 1] == y[j]){
                    conta[j, ] <- conta[j, ] + 1
               }
          }
     }

     list(conta = conta)
}
