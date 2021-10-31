# Trajectory force solver
# Import Axes3D for 3d plotting capabilities
from mpl_toolkits.mplot3d import axes3d
#Import numpy for matlab-like array capabilities
import numpy as np
#Import matplotlib.pyplot for dope plotting capabilities
import matplotlib.pyplot as plt
#Import scipy.integrate for use of ODE numerical solvers
import scipy.integrate as scpI


#Set Fundamental constants that are parameters in ODE
G=6.67430*10**(-11)
m_Sun=1.9891*10**(30)


#Set fundamental length scale L and fundamental time scale T
#L is the length of the aphelion of the Earths orbit or
#   the shortest distance it is away from to the sun in meters
L=(147.09*10**6)*10**3
#T is the number of seconds in 1 earth day
T=86400.



# The fundamental system of odes is a vector equation
#             m{a}=(GMm)(-1/r){u}
#      where I have put the vector quantities in brackets
#An equivalent nondimensionalized form for this equation is
#             {a*}=(GMmT^2/L^3)(-1/r*){u}
# or
#             {a*}=(gamma)(-1/r*){u} where gamma=GMmT^2/L^3
#This is the system that is solved in the code

#Sets up the gamma to be used in the ODE
ARGS=(G*m_Sun*T**2/L**3,)

#Sets up the system of ODEs to be solved
# dx/dt=v1, dy/dt=v2, dz/dt=v3, 
# dv1/dt=(gamma/r^2)(-x/r), dv3/dt=(gamma/r^2)(-y/r), dv3/dt=(gamma/r^2)(-z/r)
def odes(t, xi, gamma=ARGS[0]):
    x, y, z, v1, v2, v3= xi
    r=np.sqrt(x**2+y**2+z**2)
    mag=(gamma/r**2)
    return [v1,v2,v3,mag*(-x/r),mag*(-y/r),mag*(-z/r) ]



# Sets the initial and final time. Measured in days
T0=0.;

# Tf=87.969; #Period of Mercury
# Tf=224.7; #Period of Venus
Tf=365; #Period of Earth
# Tf=687; #Period of Mars


#Set the initial starting points of each planet plotted. 
#Relative to orbital plane of Earth (Equatorial plane)
#Factors in angle of inclination from equatorial plane

#Note: Code assumes that the perihelions line up, which is not true!
#       Thus, the code does not factor in
#         shift due to the argument of periapsis, however,
#        this can be easily done with slight modification to ICs
#      Does not factor in precesion of Mercury either. 
#          This can also be addressed using reletavistic correction methods 


#Parameters

#Mercury Paramaters
r0_Merc=(46.00*10**6)*10**3
theta_Merc=(7./180)*np.pi

x0_Merc=r0_Merc*np.cos(theta_Merc)/L
z0_Merc=r0_Merc*np.sin(theta_Merc)/L

v0_Merc=(T/L)*(58.98*10**3)


#Venus Paramaters
r0_Venus=(107.48*10**6)*10**3
theta_Venus=(3.39/180)*np.pi

x0_Venus=r0_Venus*np.cos(theta_Venus)/L
z0_Venus=r0_Venus*np.sin(theta_Venus)/L

v0_Venus=(T/L)*(35.26*10**3)

#Earth Parameters

x0_Earth=1.
v0_Earth=(T/L)*(30.29*10**3)


# #Mars Parameters
# r0_Mars=(206.7*10**6)*10**3
# theta_Mars=(3.39/180)*np.pi

# x0_Mars=r0_Mars*np.cos(theta_Mars)/L
# z0_Mars=r0_Mars*np.sin(theta_Mars)/L

# v0_Mars=(T/L)*(26.5*10**3)


#Set initial conditions

x_0=np.array([x0_Merc,x0_Venus,x0_Earth]); y_0=np.array([0.,0.,0.]); z_0=np.array([z0_Merc,z0_Venus,0.])
v1_0=np.array([0.,0.,0.]); v2_0=np.array([v0_Merc,v0_Venus,v0_Earth]); v3_0=np.array([0.,0.,0.])



# x_0=np.array([x0_Merc,x0_Venus,x0_Earth,x0_Mars]); y_0=np.array([0.,0.,0.,0.]); z_0=np.array([z0_Merc,z0_Venus,0.,z0_Mars])
# v1_0=np.array([0.,0.,0.,0.]); v2_0=np.array([v0_Merc,v0_Venus,v0_Earth,v0_Mars  ]); v3_0=np.array([0.,0.,0.,0.])






#Solve IVP using scipys solve_ivp program. Similar to Matlabs ode45
# Numerical Algothim used is Runge-Kutta-Felberg 45 (or RKF45) and
# can be changed if the problem is stiff
sols = [scpI.solve_ivp(odes, [T0, Tf], [x_0[i], y_0[i], z_0[i],v1_0[i],v2_0[i],v3_0[i]], method='RK45', args=ARGS,
                dense_output=False,atol=1e-12,rtol=1e-12) for i in np.arange(x_0.size)]

#Store Solution Data for easy access
R=[np.array([sols[i].y[0,:],sols[i].y[1,:],sols[i].y[2,:]]) for i in np.arange(x_0.size)]




#Plot
fig = plt.figure('Trajectories')
plt.clf()

ax = fig.gca(projection='3d')

#Plot the Trajectories
LW=1.5;AL=.3;
Trajectories=[ax.plot(R[i][0],R[i][1],R[i][2],linewidth=LW,alpha=AL) for i in np.arange(x_0.size)]


# numdiv=10;
# pointsI=[np.arange(0,R[i][0].size,int(R[i][0].size/numdiv)) for i in np.arange(x_0.size)]
# NumericalPoints=[ax.scatter(R[i][0,pointsI[i]],R[i][1,pointsI[i]],R[i][2,pointsI[i]],s=10,alpha=.3) for i in np.arange(x_0.size)]


ICplot= ax.scatter(x_0,y_0,z_0,s=15,c='r',alpha=1)

EndPlot=[ax.scatter(R[i][0,-1],R[i][1,-1],R[i][2,-1],s=15,c='magenta',alpha=1) for i in np.arange(x_0.size)]

ax.scatter(0,0,0,s=30,c='orange')


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_ylabel('y')


ax.set_zlim([-.2,.2])


plt.show()


