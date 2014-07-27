from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


def modified_GS(M):
    A=M.copy()
    n=len(A.T)
    m=len(A)
    q=np.zeros((m,n))
    r=np.zeros((n,n))
    for k in range(0,n):
        a=A.T[k]
        
        r[k][k]=la.norm(a)
        
        q.T[k]=a/r[k][k]
        
        c=q.T[k]
        for j in range(k+1,n):
            
            r[k][j]=np.dot(c.T,A.T[j])
            #print r[k][j]
            A.T[j]=A.T[j]-np.dot(r[k][j],c)
            
    return q,r
def back_subsitute(U, bb):
    n = U.shape[1]
    x = np.zeros(n)
    for j in range(n - 1, -1, -1): # loop backwards over columns
        if U[j, j] == 0:
            raise RuntimeError("singular matrix")
        x[j] = bb[j] / U[j, j]
        for i in range(0, j):
            bb[i] -= U[i, j] * x[j]
    return x
x=np.arange(0,345)
#Degree1_GS
A=np.zeros((345,2))
for i in range(0,345):
    
        A[i][0]=1
        A[i][1]=i
        
v = np.genfromtxt("/Users/lrao/Documents/UIUC_stuff/UIUC_Acads/CS450/HW2_LAKSHMI_GURURAJA_RAO/Price_of_Gasoline.txt", delimiter="\n")


C0,D0= modified_GS(A)
B0=np.dot(C0.T,v)
sol0=back_subsitute(D0,B0)
rr1=la.norm(np.dot(A,sol0)-v)/la.norm(v)
print "Degree 1 Relative Residual"
print rr1
y0=sol0[0]+sol0[1]*x
plt.plot(v)
plt.plot(x,y0)

#Quadratic_GS
A=np.zeros((345,3))
for i in range(0,345):
    
        A[i][0]=1
        A[i][1]=i
        A[i][2]=i**2


C,D= modified_GS(A)
B=np.dot(C.T,v)
sol=back_subsitute(D,B)
y=sol[0]+sol[1]*x+sol[2]*(x**2)
plt.plot(x,y)
rr2=la.norm(np.dot(A,sol)-v)/la.norm(v)
print "Degree 2 Relative Residual"
print rr2

#Degree3_GS
A=np.zeros((345,4))
for i in range(0,345):
    
        A[i][0]=1
        A[i][1]=i
        A[i][2]=i**2
        A[i][3]=i**3



C1,D1= modified_GS(A)
B1=np.dot(C1.T,v)

sol1=back_subsitute(D1,B1)

y1=sol1[0]+sol1[1]*x+sol1[2]*(x**2)+sol1[3]*(x**3)
plt.plot(x,y1)
rr3=la.norm(np.dot(A,sol1)-v)/la.norm(v)
print "Degree 3 Relative Residual"
print rr3

#Degree4_GS
A=np.zeros((345,5))
for i in range(0,345):
    
        A[i][0]=1
        A[i][1]=i
        A[i][2]=i**2
        A[i][3]=i**3
        A[i][4]=i**4



C2,D2= modified_GS(A)
B2=np.dot(C2.T,v)
sol2=back_subsitute(D2,B2)
y2=sol2[0]+sol2[1]*x+sol2[2]*(x**2)+sol2[3]*(x**3)+sol2[4]*(x**4)
plt.plot(x,y2)
rr4=la.norm(np.dot(A,sol2)-v)/la.norm(v)
print "Degree 4 Relative Residual"
print rr4

#Degree5_GS
A=np.zeros((345,6))
for i in range(0,345):
    
        A[i][0]=1
        A[i][1]=i
        A[i][2]=i**2
        A[i][3]=i**3
        A[i][4]=i**4
        A[i][5]=i**5


C3,D3= modified_GS(A)
B3=np.dot(C3.T,v)
sol3=back_subsitute(D3,B3)
y3=sol3[0]+sol3[1]*x+sol3[2]*(x**2)+sol3[3]*(x**3)+sol3[4]*(x**4)+sol3[5]*(x**5)
plt.plot(x,y3)
plt.show()
rr5=la.norm(np.dot(A,sol3)-v)/la.norm(v)
print "Degree 5 Relative Residual"
print rr5