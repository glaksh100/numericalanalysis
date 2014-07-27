from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
def generate_equispaced(n,x_1=-1,x_n=1):
    x=np.linspace(x_1,x_n,n)
    return x
def generate_chebyshev(n,x_1=-1,x_n=1):
    x=[]
    for i in range(0,n):
        val=np.cos((2*i+1)*np.pi/(2*n))
        x.append(val)
    return np.array(x)
def generate_monomial_vander(nodes):
    n=len(nodes)
    V=np.zeros([n,n])
    for i in range(0,n):
        for j in range(0,n):
            V[i][j]=nodes[i]**j
    return V
def generate_cheby_vander(nodes):
    n=len(nodes)
    V=np.zeros([n,n])
    for i in range(0,n):
        for j in range(0,n):
            val=nodes[i]
            V[i][j]=np.cos(j*np.arccos(val))
    return V
masternodes=np.linspace(5,100,20)
es=[]
ch=[]
V1=[]
c1=[]
V2=[]
c2=[]
V3=[]
c3=[]
V4=[]
c4=[]
for i in range(0,len(masternodes)):
    es.append(generate_equispaced(masternodes[i]))
    ch.append(generate_chebyshev(int(masternodes[i])))
    V1.append(generate_monomial_vander(es[i]))
    c1.append(np.log(la.cond(V1[i])))
    V2.append(generate_monomial_vander(ch[i]))
    c2.append(np.log(la.cond(V2[i])))
    V3.append(generate_cheby_vander(es[i]))
    c3.append(np.log(la.cond(V3[i])))
    V4.append(generate_cheby_vander(ch[i]))
    c4.append(np.log(la.cond(V4[i])))

pt.plot(masternodes,c1,label='Equispaced nodes monomial')
pt.plot(masternodes,c2,label='Chebyshev nodes monomial')
pt.plot(masternodes,c3,label='Equispaced nodes Chebyshev poly')
pt.plot(masternodes,c4,label='Chebyshev nodes Chebyshev poly')
pt.xlabel("Nodes")
pt.ylabel("Log of Cond")
pt.legend(loc=4)
pt.show()