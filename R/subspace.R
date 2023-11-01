#   SUBSPACE(A,B) finds the angle between two subspaces specified by the
#   columns of A and B.
#
#   If the angle is small, the two spaces are nearly linearly dependent.
#   In a physical experiment described by some observations A, and a second
#   realization of the experiment described by B, SUBSPACE(A,B) gives a
#   measure of the amount of new information afforded by the second
#   experiment not associated with statistical errors of fluctuations.
#
#   Class support for inputs A, B:
#      float: double, single

#   The algorithm used here ensures that small angles are computed
#   accurately, and it allows subspaces of different dimensions following
#   the definition in [2]. The first issue is crucial.  The second issue is
#   not so important; but since the definition from [2] coinsides with the
#   standard definition when the dimensions are equal, there should be no
#   confusion - and subspaces with different dimensions may arise in
#   problems where the dimension is computed as the numerical rank of some
#   inaccurate matrix.

#   References:
#   [1] A. Bjorck & G. Golub, Numerical methods for computing
#       angles between linear subspaces, Math. Comp. 27 (1973),
#       pp. 579-594.
#   [2] P.-A. Wedin, On angles between subspaces of a finite
#       dimensional inner product space, in B. Kagstrom & A. Ruhe (Eds.),
#       Matrix Pencils, Lecture Notes in Mathematics 973, Springer, 1983,
#       pp. 263-285.

#   Thanks to Per Christian Hansen
#   Copyright 1984-2007 The MathWorks, Inc. 
#   $Revision: 5.10.4.3 $  $Date: 2007/09/18 02:15:38 $

# Compute orthonormal bases, using SVD in "orth" to avoid problems
# when A and/or B is nearly rank deficient.


subspace=function(A,B){
  A = qr.Q(qr(A));
  B = qr.Q(qr(B));

  #Check rank and swap
  if(ncol(A)<ncol(B))  {tmp = A; A = B; B = tmp}
  # Compute the projection the most accurate way, according to [1].
  for(k in 1:ncol(A)){
    B = B - A[,k]%*%(t(A[,k])%*%B)
  }
  # Make sure it's magnitude is less than 1.
  theta = asin(min(1,base::norm(B,type="2")))*180/pi
  theta
}


# orth<-function(A){
#   # ORTH   Orthogonalization.
#   #   Q = ORTH(A) is an orthonormal basis for the range of A.
#   #   That is, Q'*Q = I, the columns of Q span the same space as 
#   #   the columns of A, and the number of columns of Q is the 
#   #   rank of A.
#   #
#   #   Class support for input A:
#   #      float: double, single
#   #
#   #   See also SVD, RANK, NULL.
#   
#   #   Copyright 1984-2011 The MathWorks, Inc. 
#   #   $Revision: 5.11.4.3 $  $Date: 2012/01/19 16:59:31 $
#   A=as.matrix(A)
#   tmp = svd(A) #S is always square.  
#   Q=tmp$v
#   S=tmp$d
#   S = diag(S)
#   tol = max(dim(A) * S[1] * 2.2204e-16)
#   r = sum(S > tol)
#   Q = as.matrix(Q[,1:r])
#   Q
# }