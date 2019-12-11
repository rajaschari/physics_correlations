# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from __future__ import division
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import cmath
from random import *

m1 = 0.11;
m2 = 0.12;
m3 = 0.13;
m4 = 0.14;
m5 = 0.15;


J = 1;

mu = 1;
d = 0.2*J;

Bx = 0.12;
By = 0.12;
Bz = 0.12;

Ex = 0;
Ey = 0;
Ez = 0;

n1 = np.array([1/2, math.sqrt(3)/2]);
n2 = np.array([-1/2, math.sqrt(3)/2]);


ex = [-1/2, math.sqrt(3)/2];
ey = [-1/2, -math.sqrt(3)/2];
ez = [1,0];

Qx = Bz*(m1*(Ez - Ex) + m2*(Ez + Ex) + m5*(Ex)) + By*(m1*(Ex - Ey) 
    - m2*(Ex + Ey) - m5*Ex);
         
Qy = Bx*(m1*(Ex - Ey) + m2*(Ex + Ey) + m5*(Ey)) + Bz*(m1*(Ey - Ez) 
    - m2*(Ey + Ez) - m5*Ey);

Qz = By*(m1*(Ey - Ez) + m2*(Ey + Ez) + m5*(Ez)) + By*(m1*(Ez - Ex) 
    - m2*(Ez + Ex) - m5*Ez);         

sh1=random()*10**(-8);
sh2=random()*10**(-8);
                  
sh=[sh1,sh2];

b1=np.array([(2*math.pi)+sh1,(2*math.pi/math.sqrt(3))+sh2]);
    
b2=np.array([(-2*math.pi)+sh1,(2*math.pi/math.sqrt(3))+sh2]);
    
c1=0;
c2=0;

#k=(c1/20)*b1+(c2/20)*b2;
L=20;
   
def k(c1,c2,L):
    return (c1/L)*b1+(c2/L)*b2;
#    return [x + y + z for x, y, z in zip(map(lambda x: x * float((c1/L)), b1), map(lambda x: x * float((c2/L)), b1),map(lambda x: x , sh))];
             
             
kap = 2*Bx*By*Bz*(mu**3)/(d**2);
                 
j=complex(0,1)
               
#f[k_] := 2*J (Exp[I*k.n1] + Exp[I*k.n2] + 1)
#\[CapitalDelta][k_] := 
# 4*\[Kappa] (Sin[k.n1] + Sin[k.(-n2)] + Sin[k.(n2 - n1)])
#G[k_] := (4*\[Mu])*(Subscript[Q, x]*Sin[k.ex] + 
#     Subscript[Q, y]*Sin[k.ey] + Subscript[Q, z]*Sin[k.ez])/\[Delta]
#\[Epsilon][k_] := G[k] - fe[k];
#d[k_] := \[CapitalDelta][k] + fo[k];
 
# k = [k1,k2];

#n2_=(map(lambda x: x * (-1), n2));    
#n21=[x - y for x, y in zip(n2, n1)]

def f(c1,c2,L):
     return (2*J* (cmath.exp(j*(np.dot(k(c1,c2,L),n1))) + cmath.exp(j*np.dot(k(c1,c2,L),n2)) + 1));
 
def D(c1,c2,L):
     return 4*kap*(cmath.sin(np.dot(k(c1,c2,L),n1)) + cmath.sin(np.dot(k(c1,c2,L),-n2)) + cmath.sin(np.dot(k(c1,c2,L),n2-n1)));

def g(c1,c2,L):
    return 4*mu*(Qx*cmath.sin(np.dot(k(c1,c2,L),ex)) + Qy*cmath.sin(np.dot(k(c1,c2,L),ey)) + Qz*cmath.sin(np.dot(k(c1,c2,L),ez)));

def fe(c1,c2,L):
    return (f(c1,c2,L)+f(-c1,-c2,L))/2;

def fo(c1,c2,L):
    return (f(c1,c2,L)-f(-c1,-c2,L))/2;

def E(c1,c2,L):
    return g(c1,c2,L)-cmath.sqrt(fe(c1,c2,L)**2 - fo(c1,c2,L)**2 + 4*(D(c1,c2,L)**2))/2;
 
#def u(c1,c2,L):
#    return (fe(c1,c2,L) + cmath.sqrt((fe(c1,c2,L)**2) - (fo(c1,c2,L)**2) + 
#        4*(D(c1,c2,L)**2)))/((fo(c1,c2,L) - 2*D(c1,c2,L)+sh1)* cmath.sqrt(1 + abs((fe(c1,c2,L)
#        + cmath.sqrt((fe(c1,c2,L)**2)- (fo(c1,c2,L)**2) + 4*(D(c1,c2,L)**2)))/(sh1+fo(c1,c2,L) - 2*D(c1,c2,L)))));
#    

def u(c1,c2,L):
    return (fe(c1,c2,L) + 2*(g(c1,c2,L)-E(c1,c2,L)))/(cmath.sqrt(abs(fo(c1,c2,L) - 2*D(c1,c2,L))**2 + abs((fe(c1,c2,L)
        + 2*(sh1+g(c1,c2,L)-E(c1,c2,L))))));
    
#np.angle(fo(c1,c2,L)-2*D(c1,c2,L))* 

#def v(c1,c2,L):
#    return (fe(c1,c2,L) - cmath.sqrt((fe(c1,c2,L)**2) - (fo(c1,c2,L)**2) + 
#        4*(D(c1,c2,L)**2)))/((fo(c1,c2,L) - 2*D(c1,c2,L)+sh1)* cmath.sqrt(1 + abs((fe(c1,c2,L)
#        - cmath.sqrt((fe(c1,c2,L)**2)- (fo(c1,c2,L)**2) + 4*(D(c1,c2,L)**2)))/(sh1+fo(c1,c2,L) - 2*D(c1,c2,L)))));
    
def v(c1,c2,L):
    return (fe(c1,c2,L) - 2*(g(c1,c2,L)-E(c1,c2,L)))/( cmath.sqrt(abs(fo(c1,c2,L) - 2*D(c1,c2,L))**2  + abs(fe(c1,c2,L)
        - 2*(sh1+g(c1,c2,L)-E(c1,c2,L)))));

#np.angle(fo(c1,c2,L)-2*D(c1,c2,L))* 

         
x=[];
y=[];
z=[];

#  
#for c1 in range(0,20):
#    
#    for c2 in range(0,20):
#        
#        x.append(c1);
#        y.append(c2);
##        z.append(abs(cmath.sin(np.dot(k(c1,c2),n2_))));
##        plt.plot(x,y)
#        z.append(abs(u(c1,c2,10)));
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#ax.scatter( x, y, z)

G1=0;
G2=0;

c11=0;
c12=0;

def X01(c1,c2,G1,G2,L):
    
    c11=c1-L*(G1+math.sqrt(3)*G2)/(4*math.pi);
    c12=c2-L*(-G1+math.sqrt(3)*G2)/(4*math.pi);    
    
    return v(c11,c12,L).conjugate()*u(c1,c2,L)*(u(-c11,-c12,L).conjugate()*v(-c1,-c2,L) -v(c11,c12,L)*u(c1,c2,L).conjugate())*np.heaviside(-abs(E(-c11,-c12,L)),0)*np.heaviside(-abs(E(c1,c2,L)),0)/(E(-c11,-c12,L) + E(c1,c2,L));

def X02(c1,c2,G1,G2,L):
    
    c11=c1-L*(G1+math.sqrt(3)*G2)/(4*math.pi);
    c12=c2-L*(-G1+math.sqrt(3)*G2)/(4*math.pi);    
    
    return v(c11,c12,L).conjugate()*v(c1,c2,L)*(u(-c11,-c12,L).conjugate()*u(-c1,-c2,L) - v(c1,c2,L).conjugate()*v(c11,c12,L))*np.heaviside(-abs(E(-c11,-c12,L)),0)*np.heaviside(abs(E(-c1,-c2,L)),0)/(E(-c11,-c12,L) + E(-c1,-c2,L));
 
def X03(c1,c2,G1,G2,L):
    
    c11=c1-L*(G1+math.sqrt(3)*G2)/(4*math.pi);
    c12=c2-L*(-G1+math.sqrt(3)*G2)/(4*math.pi);    
    
    return u(c11,c12,L).conjugate()*u(c1,c2,L)*(v(-c11,-c12,L).conjugate()*v(-c1,-c2,L) -  u(c11,c12,L)*u(c1,c2,L).conjugate())*np.heaviside(abs(E(c11,c12,L)),0)*np.heaviside(-abs(E(c1,c2,L)),0)/(E(c11,c12,L) + E(c1,c2,L))

def X04(c1,c2,G1,G2,L):
    
    c11=c1-L*(G1+math.sqrt(3)*G2)/(4*math.pi);
    c12=c2-L*(-G1+math.sqrt(3)*G2)/(4*math.pi);    
    
    return u(c11,c12,L).conjugate()*v(c1,c2,L)*(v(-c11,-c12,L).conjugate()*u(-c1,-c2,L) - u(c11,c12,L)*v(c1,c2,L).conjugate())*np.heaviside(abs(E(c11,c12,L)),0)*np.heaviside(abs(E(-c1,-c2,L)),0)/(E(c11,c12,L) + E(-c1,-c2,L));

def X0(c1,c2,G1,G2,L):   
             
    return X01(c1,c2,G1,G2,L) + X02(c1,c2,G1,G2,L) + X03(c1,c2,G1,G2,L) + X04(c1,c2,G1,G2,L);


def X(G1,G2,L):
    var=0;
    for c1 in range(0,L):
        for c2 in range(0,L):
            
            var=var+X0(c1,c2,G1,G2,L)

    return var/(L**2)

x1=[];
x2=[];
x3=[];
pp=[];
   
#             
#pp=np.linspace(0,2*np.pi,20);
#
#for x in range(len(pp)):
#    
#    G2=pp[x];
#    
#    x1.append(G2);
#    x2.append(X(0,G2,10))
#
#plt.plot(x1,x2)

pp=np.linspace(0,4*np.pi,40);

for x in range(len(pp)):
    G1=pp[x];
    for y in range(len(pp)):
        G2=pp[y];
             
        x1.append(G1);
        x2.append(G2);
#        z.append(abs(cmath.sin(np.dot(k(c1,c2),n2_))));
#        plt.plot(x,y)
        x3.append(abs(X(G1,G2,5)));
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#ax.scatter( x1, x2, x3)

ax.plot_trisurf(x1,x2,x3)

#Axes3D.plot_trisurf(x1,x2,x3)




#(v(c11,c12).conjugate()*u(c1,c2)*(u(-c11,-c12).conjugate()*v(-c1,-c2) - 
#            v(c11,c12)*u(c1,c2).conjugate())*np.heaviside(-E(-c11,-c12),0)*np.heaviside(-E(c1,c2))/(E(-c11,-c12) + E(c1,c2)) +
#            (v(c11,c12).conjugate()*v(c1,c2))*(u(-c11,-c12).conjugate()*u(-c1,-c2) - 
#            v(c1,c2).conjugate()*v(c11,c12))*np.heaviside(-E(-c11,-c12),0)*np.heaviside(E(-c1,-c2))/(E(-c11,-c12) + E(-c1,-c2)) +
#            (u(c11,c12).conjugate()*u(c1,c2))*(v(-c11,-c12).conjugate()*v(-c1,-c2) - 
#            u(c11,c12)*u(c1,c2).conjugate())*np.heaviside(E(c11,c12),0)*np.heaviside(-E(c1,c2),0)/(E(c11,c12) + E(c1,c2)) +
#            (u(c11,c12).conjugate()*v(c1,c2))*(v(-c11,-c12).conjugate()*u(-c1,-c2) - u(c11,c12)*v(c1,c2).conjugate())*np.heaviside(E(c11,c12))*np.heaviside(E(-c1,-c2),0)/(E(c11,c12) + E(-c1,-c2));
 
 
