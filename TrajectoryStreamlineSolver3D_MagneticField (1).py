
# Import Axes3D for 3d plotting capabilities
from mpl_toolkits.mplot3d import axes3d
#Import numpy for matlab-like array capabilities
import numpy as np
#Import matplotlib.pyplot for dope plotting capabilities
import matplotlib.pyplot as plt
#Import scipy.integrate for use of ODE numerical solvers
import scipy.integrate as scpI


#Sets up the system of ODEs to be solved
# dx/dt=P, dy/dt=Q, dz/dt=R, 
def odes(t, xi):
    x, y, z= xi
    w=x**2+y**2-1.
    l=z**2+w**2
    return [(-2.*x*z)/l,(-2.*y*z)/l,w/l ]

# Sets the initial and final time. 
T0=0.;Tf=10; 


#Set domain for initial conditions. Centered at point xC,yC,zC in a cube of L**3
L=1./(np.sqrt(3));
xC=1.5;yC=1.5;zC=1.5;

# L=2.;
# xC=0;yC=0.;zC=0.;


l1=xC-L/2.;r1=xC+L/2.;
l2=yC-L/2.;r2=yC+L/2.;
l3=zC-L/2.;r3=zC+L/2.;


n=3;#Number of points sampled along x, y, and z, n^3 is the total number of streamlines plotted be careful not to raise too high

x = np.linspace(l1,r1,n); y = np.linspace(l2,r2,n); z = np.linspace(l3,r3,n)
X, Y, Z = np.meshgrid(x, y,z, indexing='ij', sparse=False)

#This line removes all of the ``bad'' initial conditions where the vector field is not defined.
shp=np.where(Z**2+(X**2+Y**2-1.)**2 > 0)#shp=np.where(X**2+Y**2>.001)
X=X[shp];Y=Y[shp];Z=Z[shp]


#Reshaping the previous arrays into 1D vectors
x_0=np.reshape(X,X.size); y_0=np.reshape(Y,Y.size);z_0=np.reshape(Z,Z.size)




#Solve IVP using scipys solve_ivp program. Similar to Matlabs ode45
# Numerical Algothim used is Runge-Kutta-Felberg 45 (or RKF45) and
# can be changed if the problem is stiff
sols = [scpI.solve_ivp(odes, [T0, Tf], [x_0[i], y_0[i], z_0[i]], method='RK45', 
                dense_output=False,atol=1e-12,rtol=1e-12) for i in np.arange(x_0.size)]

#Store Solution Data for easy access
R=[np.array([sols[i].y[0,:],sols[i].y[1,:],sols[i].y[2,:]]) for i in np.arange(x_0.size)]






#Plot
fig = plt.figure('Trajectories')
plt.clf()

ax = fig.gca(projection='3d')


#This creates the color palette for the curves and guarantees that each streamline has a unique color
plt.gca().set_prop_cycle(plt.cycler('color',plt.cm.jet(np.linspace(0,1,n**3))))

#This plots the steamlines
LW=1.5;AL=.3;
Trajectories=[ax.plot(R[i][0],R[i][1],R[i][2],linewidth=LW,alpha=.8) for i in np.arange(x_0.size)]


# #This plots the vector field on the streamlines
# numdiv=3;
# pointsI=[np.arange(0,R[i][0].size,int(R[i][0].size/numdiv)) for i in np.arange(x_0.size)]
# TangVect=[odes(0.0,[R[i][0,pointsI[i]],R[i][1,pointsI[i]],R[i][2,pointsI[i]]]) for i in np.arange(x_0.size)]
# VectorField=[ax.quiver(R[i][0,pointsI[i]],R[i][1,pointsI[i]],R[i][2,pointsI[i]],TangVect[i][0],TangVect[i][1],TangVect[i][2],color='k',length=.5,arrow_length_ratio = 0.1) for i in np.arange(x_0.size)]
# ICplot= ax.scatter(x_0,y_0,z_0,s=15,c='r',alpha=1)


#This plots the initial conditions and the endpoints
ICplot= ax.scatter(x_0,y_0,z_0,s=15,c='r',alpha=1)
EndPlot=[ax.scatter(R[i][0,-1],R[i][1,-1],R[i][2,-1],s=15,c='magenta',alpha=1) for i in np.arange(x_0.size)]



# This plots the circle where the vector field is undefined.
tt=np.linspace(0,2*np.pi,1000);
ax.plot(np.cos(tt),np.sin(tt),np.zeros(tt.size),color='k',linewidth=5,alpha=.3)



#This sets labels on the axes
ax.set_xlabel('x');ax.set_ylabel('y');ax.set_ylabel('y')

plt.show()
