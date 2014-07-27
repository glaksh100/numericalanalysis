from __future__ import division
import numpy as np
import scipy as sp
import numpy.linalg as la
import math

def CholeskyFactorization(M):
    A=M.copy()
    n=len(A)
 
    for j in range(0,n):
        
        for k in range(0,j):
            for i in range(j,n):
                A[i][j]=A[i][j]-A[i][k]*A[j][k]
                
        A[j][j]=math.sqrt((A[j][j]))
        
        for k in range(j+1,n):
            A[k][j]=A[k][j]/A[j][j]
    
    for i in range(0,n):
         for j in range(i+1,n):
             A[i][j]=0

    return A        

def rand_spd(n):
        A = np.random.rand(n,n)
        return np.dot(A, A.T)

X=np.array([[5,0,2.5],[0,2.5,0],[2.5,0,2.125]])
Y=CholeskyFactorization(X)
z=la.norm(np.dot(Y,Y.T)-X)/la.norm(X)
print "Test Case 1"

X=rand_spd(20)

Y=CholeskyFactorization(X)

relerror=la.norm(np.dot(Y,Y.T)-X)/la.norm(X)
c=la.cond(X)
print "Relative Error"
print relerror
print "Condition Number"

print c

print "Test Case 2"

X=rand_spd(20)

Y=CholeskyFactorization(X)

relerror=la.norm(np.dot(Y,Y.T)-X)/la.norm(X)
c=la.cond(X)
print "Relative Error"
print relerror
print "Condition Number"
print c
print "Test Case 3"

X=rand_spd(20)

Y=CholeskyFactorization(X)

relerror=la.norm(np.dot(Y,Y.T)-X)/la.norm(X)
c=la.cond(X)
print "Relative Error"
print relerror
print "Condition Number"
print c
