from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt

def cubic_spline(x,y):
	n=len(x)
	var=4*(n-1)
	print var/2
	A=np.zeros([var,var])
	k=0
	i=0
	j=0
	while(i<var/2):
		p=3
		while(p>=0):
			A[i][j]=x[k]**p
			p-=1
			j+=1
		if(i%2==0):
			k+=1
			j=j-4
		i+=1
	i=var/2
	k=1
	j=0
	while(i<14):

			A[i][j]=3*x[k]**2
			A[i][j+1]=2*x[k]
			A[i][j+2]=1
			A[i][j+4]=-3*x[k]**2
			A[i][j+5]=-2*x[k]
			A[i][j+6]=-1
			i+=1
			j+=4
			k+=1
			#print j
	i=14
	k=1
	j=0
	while(i<18):

			A[i][j]=6*x[k]
			A[i][j+1]=2
			A[i][j+4]=-6*x[k]
			A[i][j+5]=-2
			i+=1
			j+=4
			k+=1
	A[18][0]=6*x[0]
	A[18][1]=2
	A[19][16]=6*x[5]
	A[19][17]=2
	b=np.zeros([20,1])
	b[0][0]=y[0]
	b[9][0]=y[5]
	c=1
	k=1
	while(c<9):
		b[c][0]=y[k]
		if(c%2==0):
			k+=1
		c+=1
	sol=la.solve(A,b)
	return sol
def polynomial_ret(M,x,index):
	return (M[index[0]]*(x**3)+M[index[1]]*(x**2)+M[index[2]]*(x)+M[index[3]])

x1=np.random.ranf(6)
y1=np.random.ranf(6)
m=np.sort(x1)
B=cubic_spline(m,y1)
ind1=np.arange(0,4,1)
ind2=np.arange(4,8,1)
ind3=np.arange(8,12,1)
ind4=np.arange(12,16,1)
ind5=np.arange(16,20,1)
x_1=np.arange(m[0],m[1],0.01)
x_2=np.arange(m[1],m[2],0.01)
x_3=np.arange(m[2],m[3],0.01)
x_4=np.arange(m[3],m[4],0.01)
x_5=np.arange(m[4],m[5],0.01)
p1=polynomial_ret(B,x_1,ind1)
p2=polynomial_ret(B,x_2,ind2)
p3=polynomial_ret(B,x_3,ind3)
p4=polynomial_ret(B,x_4,ind4)
p5=polynomial_ret(B,x_5,ind5)
pt.plot(x_1,p1,label='p1')
pt.plot(x_2,p2,label='p2')
pt.plot(x_3,p3,label='p3')
pt.plot(x_4,p4,label='p4')
pt.plot(x_5,p5,label='p5')
pt.plot(m,y1,label='data points')
pt.xlabel("x values")
pt.ylabel("y values")
pt.legend()
pt.show()
