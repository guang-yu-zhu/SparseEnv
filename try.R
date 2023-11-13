source('sourceDir.R')
sourceDir('R/')

data(Berkeley)
X = Berkeley$X
Y = Berkeley$Y[, c(1, 2, 21, 23)]
spenvbase(X,Y,1,lambda = 0.1)
