topol <- c('core', 'socket', 'node')
howmany <- 100 
iter <- howmany -1
for  ( x in topol ){
  setwd(paste("~/Pubblici/learn-git/data/topo-",x,sep=''))
  mat <- array(data = 0, dim = c(howmany,30,2))
  for (y in 0:iter){
    tab <- read.csv(paste('results-',x,'-',y,'.csv', sep = ''))
    mat[y+1,,1] <- tab[,3]
    mat[y+1,,2] <- tab[,4]
  }
  mat.mean <- matrix(data = 0, nrow = 30, ncol = 2)
  mat.var <- matrix(data = 0, nrow = 30, ncol = 2)
  for (y in 1:30){
    mat.mean[y,1] <- mean(mat[,y,1])
    mat.mean[y,2] <- mean(mat[,y,2])
    mat.var[y,1] <- var(mat[,y,1])
    mat.var[y,2] <- var(mat[,y,2])
  }
  tab[,3] <- mat.mean[,1]
  tab[,4] <- mat.mean[,2]
  write.csv2(tab,'results-',x,'-summary.csv',
            quote = FALSE ,row.names = FALSE)
}

