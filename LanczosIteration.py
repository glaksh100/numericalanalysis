from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
def lanczos_iteration(A):
    n=A.shape[1]
    #print n
    Q=np.zeros((n,n))
    H=np.zeros((n,n))
    beta=np.zeros(n)
    x0=np.random.randn(n)
    #print x0
   
    Q[:,0]=x0/la.norm(x0)
    
    beta[0]=0
    for k in range(0,n-1):
        qk=Q[:,k]
        u=np.dot(A,qk)
        
        alpha=np.dot(qk.T,u)
     
        u=u-np.dot(beta[k-1],Q[:,k-1])-np.dot(alpha,Q[:,k])
        beta[k]=la.norm(u)
        H[k][k]=alpha
        H[k+1][k]=H[k][k+1]=beta[k]
        if beta[k]==0:
            break
        Q[:,k+1]=u/beta[k]
        qk1=u/beta[k]
        alpha1=np.dot(qk1.T,np.dot(A,qk1))
    H[n-1][n-1]=alpha1
    #print alpha1
    
    return Q,H
n = 32
B=np.random.randn(n,n)
Q1,R1=la.qr(B)
D1=np.arange(1,33,1)
D=np.diag(D1)
A=np.dot(np.dot(Q1,D),Q1.T)

Q,H=lanczos_iteration(A)
ritz=[]

#pt.spy(H)
j=1
for j in range(1,n):
    R=np.zeros((j,j))
    R=H[:j,:j]
    ritz.append(la.eigvals(R))
for i, rv in enumerate(ritz):
    pt.plot([i] * len(rv), rv, "x")
pt.show()

