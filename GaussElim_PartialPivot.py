from __future__ import division

import numpy as np
import numpy.linalg as la
def column(matrix, i):
    return [row[i] for row in matrix]
def back_subst(A,b):
 
    n=len(A)
    x=np.zeros((n,1))
    for j in range(n-1,-1,-1):
        if A[j][j]==0:
            break
        x[j]=b[j]/A[j][j]
        
        for i in range(j-1,-1,-1):
            b[i]=b[i]-A[i][j]*x[j]
    return x
def partial_pivot(A,b):
    n=len(A)
    #print n
    
    A=np.concatenate((A,b),axis=1)
    for j in range(0,n):
        c=column(A,j)
        
        for i in range(0,len(A)):
            if c[i]<0:
                c[i]=c[i]*-1
        imax=np.argmax(c);
        if(imax!=j):
            A[[j,j]]=A[[imax,j]]
       
        pivot=A[j][j]
        if pivot==0:
            break
        for i in range(j+1,n):
            #print i
            factor=A[i][j]/pivot
            
            for k in range(0,n+1):
                A[i][k]=A[i][k]-factor*A[j][k]
                #print A[i][k]
    #print A
    
    l=column(A,n)
    R=A[0:n, 0:n]
    sol=back_subst(R,l)
    return sol
n=100   

###Testing case 1
New=np.random.randn(n,n)
v=np.random.randn(n,1)
b1=np.dot(New,v)
x1=partial_pivot(New,v)
c1=la.cond(New)
rr= (la.norm(v-x1)/la.norm(b1))
err=la.norm(v-x1)
relerr=err/la.norm(v)
print "####Test Case 1"
print "Condition number of A"
print c1
print "Relative Residue"
print rr
print "Relative Error"
print relerr

###Testing case 2
New=np.zeros((n,n))
for i in range(0,n):
    for j in range(0,n):
        if((i+2)%n==j):
            New[i][j]=5
        else:
            New[i][j]=1e-3
v=np.random.randn(n,1)
b1=np.dot(New,v)
x1=partial_pivot(New,v)
c1=la.cond(New)
print "####Test Case 2"
print "Condition number of A"
print c1

rr= (la.norm(v-x1)/la.norm(b1))
err=la.norm(v-x1)
relerr=err/la.norm(v)


print "Relative Residue"
print rr
print "Relative Error"
print relerr

####Testing case 3

New=np.zeros((n,n))
for i in range(0,n):
    for j in range(0,n):
        
            New[i][j]=1/(1+np.abs(((i+1)%n)-j))**4

v=np.random.randn(n,1)
b1=np.dot(New,v)
x1=partial_pivot(New,v)
c1=la.cond(New)
rr= (la.norm(v-x1)/la.norm(b1))
err=la.norm(v-x1)
relerr=err/la.norm(v)
print "####Test Case 3"
print "Condition number of A"
print c1
print "Relative Residue"
print rr
print "Relative Error"
print relerr