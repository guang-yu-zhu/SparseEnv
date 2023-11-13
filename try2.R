#devtools::install_github('guang-yu-zhu/SparseEnv')
library(SparseEnv)
data(SAT)
X <- SAT[, 1:4]
Y <- SAT[, 5:7]
choose_spxenv(X, Y)
