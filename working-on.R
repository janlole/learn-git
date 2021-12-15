library(plotrix)
topol <- c('core', 'socket', 'node')
net <- c('gig', 'infin')
pml <- c('ob1', 'ucx')
btl <- c('vader','tcp' ,'openib')

howmany <- 100 
iter <- howmany -1
for  ( to in topol ){
  for (ne in net ){
    for (pm in pml ){
      for (bt in btl ){
        if (((pm=='ob1')&&(bt=='openib'))||
            ((to=='node')&&(pm=='ob1')&&(bt=='vader'))){ next 
        }else{
          setwd(paste("~/Pubblici/learn-git/data_col_mpi/topo-",
                      to,"/net-",ne,"/pml-",pm,"/btl-",bt,sep=''))
          mat <- array(data = 0, dim = c(howmany,30,2))
          y.dat <- c()
          x.dat <- c()
          z.dat <- c()
          for (y in 0:iter){
            tab <- read.csv(paste('results-',y,'.csv', sep = ''))
            mat[y+1,,1] <- tab[,3]
            mat[y+1,,2] <- tab[,4]
            y.dat <- c(y.dat, tab[,3])
            x.dat <- c(x.dat, tab[,1])
            z.dat <- c(z.dat, tab[,4])
          }
          modello <- nls(y.dat~ a + x.dat/b, start = list (a=2, b = 12000))
          latency <- coef(modello)[1]
          bandwidth <- coef(modello)[2]
          
          mat.mean <- matrix(data = 0, nrow = 30, ncol = 2)
          mat.var <- matrix(data = 0, nrow = 30, ncol = 2)
          for (y in 1:30){
            mat.mean[y,1] <- mean(mat[,y,1])
            mat.mean[y,2] <- mean(mat[,y,2])
            mat.var[y,1] <- var(mat[,y,1])
            mat.var[y,2] <- var(mat[,y,2])
            tab[y,5] <- latency + tab[y,1]/bandwidth
            tab[y,6] <- tab[y,1]/tab[y,5]
          }
          tab[,3] <- mat.mean[,1]
          tab[,4] <- mat.mean[,2]
          
          write.csv2(round(tab,digits = 2), 
                     file = paste('~/Pubblici/learn-git/data_col_mpi/summary--topo-',to,
                                  "--net-",ne,
                                  "--pml-",pm,
                                  "--btl-",bt,'.csv',sep=''),
                     quote = FALSE ,row.names = FALSE)

        }
      }
    }
  }
}







