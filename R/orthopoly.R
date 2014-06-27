orthopoly<-
function(marginals, scores = c(1:length(marginals))){
 
     L <- length(marginals)
     a <- matrix(0, nrow = L, ncol = L)
     A <- vector(mode = "numeric", length = L)
     B <- vector(mode = "numeric", length = L)
     C <- vector(mode = "numeric", length = L)

     for(u in 1:2) {
          A[u] <- 0
          B[u] <- 0
          C[u] <- 0
     }

     for(l in 1:L) {
          a[1, l] <- 1
          a[2, l] <- ((marginals %*% (scores^2) - (marginals %*% 
               scores)^2)^-0.5)  * (scores[l] - marginals %*% scores)
     }

     if (L > 3){
          for(u in 3:L) {
               for(l in 1:L) {
          
               B[u] <- sum(marginals * scores * (a[u - 1,  ]^2))
               C[u] <- sum(marginals * scores * a[u - 1,  ] * 
                    a[u - 2,  ])
               A[u] <- (sum(marginals * (scores^2) * (a[u - 1,  ]^2)) 
                    - (B[u]^2) - (C[u]^2))^-0.5
               a[u, l] <- A[u] * ((scores[l] - B[u]) * a[u - 1, l] 
                    - C[u] * a[u - 2, l])
               }
          }
     }

     ortho <- t(a)
     return(ortho)
}
