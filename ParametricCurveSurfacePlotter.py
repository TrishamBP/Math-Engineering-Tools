# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 13:04:46 2020

@author: willi
"""


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D




t0=0.;tf=2*np.pi

tt=np.linspace(t0,tf,1000)




# def r(t):
#     R=4.
#     omega=1.
#     a=1.
#     return np.array([R*np.cos(omega*t),R*np.sin(omega*t),a*t])




def R(u,v):
    R=2.
    r=1.
    return np.array([(R+r*np.cos(v))*np.cos(u),(R+r*np.cos(v))*np.sin(u),r*np.sin(v)])


uu=np.linspace(0,2*np.pi,1000)
vv=np.linspace(0,2*np.pi,1000)

UU,VV=np.meshgrid(uu,vv)



def r(t):
    omega=3.
    beta=4.
    return np.array([R(omega*t,beta*t)[0],R(omega*t,beta*t)[1],R(omega*t,beta*t)[2]])



def v(t,h):
    return (r(t+h/2.)-r(t-h/2.))/h

def a(t,h):
    return (v(t+h/2.,h)-v(t-h/2,h))/h





fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


rr=r(tt)

ax.plot(rr[0,:],rr[1,:],rr[2,:],color='k')

RR=R(UU,VV)

ax.plot_surface(RR[0,:,:],RR[1,:,:],RR[2,:,:],alpha=.5)





t_0=np.pi/2

H=.001;

ax.quiver(0,0,0,r(t_0)[0],r(t_0)[1],r(t_0)[2],color='b')

ax.quiver(0,0,0,r(t_0+H/2.)[0],r(t_0+H/2.)[1],r(t_0+H/2.)[2],color='g')
ax.quiver(0,0,0,r(t_0-H/2.)[0],r(t_0-H/2.)[1],r(t_0-H/2.)[2],color='g')


ax.quiver(r(t_0)[0],r(t_0)[1],r(t_0)[2],v(t_0,H)[0],v(t_0,H)[1],v(t_0,H)[2],color='r')



L=3;

ax.set_xlabel('x')
ax.set_xlim(-L,L)
ax.set_ylabel('y')
ax.set_ylim(-L,L)
ax.set_zlabel('z')
ax.set_zlim(-L,L)


#plt.axis('equal');
