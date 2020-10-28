'''
Eduard Larrañaga
Computational Astrophysics 
2020

Gravitational N-Body Code
'''


import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as Axes3D

# Relevant Constants in the units
# (years, AU, Solar_masses)
# Newtonian Gravitational Constant
G = 4.*np.pi**2

# Conversion factors
m_Sun = 1.98855e30 # Solar mass in kg
au_in_meters = 1.49598261e11 # 1 au in meters
year_in_seconds = 3600*24*365



def RK4(dt, t0, q0, mass):
    '''
    ------------------------------------------
    RK4(h, t0, q0)
    ------------------------------------------
    4th Order Runge-Kutta method for solving 
    a system of ODEs.
    Arguments:
    ODE: function defining the system of ODEs
    h: stepsize for the iteration
    t0: independent parameter initial value
    q0: numpy array with the initial values of
        the functions in the ODEs system
    ------------------------------------------
    '''
    k1 = dt*ODE(t0, q0, mass)
    k2 = dt*ODE(t0 + dt/2., q0 + k1/2., mass)
    k3 = dt*ODE(t0 + dt/2., q0 + k2/2., mass)
    k4 = dt*ODE(t0 + dt, q0 + k3, mass)
    q1 = q0 + (k1 + 2.*k2 + 2.*k3 + k4)/6.
    return q1


def ODE(t0, q0, mass):
    
    [INSERT RHS HERE]
    return Q

def TotalEnergy(q0,mass):
	[INSERT HERE ENERGY CONSERVATION]
	return Ekin + Egrav




# Read the initial data
initial_data_file = "sun_earth.asc"
(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)

# Convert from SI units to (years, AU, Solar_Mass) units
x = x/au_in_meters
y = y/au_in_meters
z = z/au_in_meters
vx = vx*year_in_seconds/au_in_meters
vy = vy*year_in_seconds/au_in_meters
vz = vz*year_in_seconds/au_in_meters
mass = mass/m_Sun

# Creation of the time grid (in years)
t_0 = 0.
t_f = 5.

# Number of steps in the grid
n = 40000

t_grid = np.linspace(t_0, t_f, n)



# Constant stepsize defined by the number of steps in the grid
dt = (t_f - t_0)/n

# Array with the initial conditions for all particles
q = np.array((x,y,z,vx,vy,vz)).transpose()

energy = np.zeros(n)

# Limits for the plot
xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(q[:,0],q[:,1],q[:,2], color=('red'), marker='.')
plt.draw()

for i in range(n):
	energy[i] = TotalEnergy(q,mass)
	ax.scatter(q[:,0],q[:,1],q[:,2], color=('black'), marker='.')
	plt.draw()
	q = RK4(dt, t_grid[i], q, mass)
	print(f'{t_grid[i]*100/t_f:.1f} %')

energychange = (energy[n-1]-energy[0])/energy[0]
print(f'The relative change in energy is: {energychange:.2f} %')
ax.set_title(r'Sun-Earth System')
ax.set_xlim((-rmax,rmax))
ax.set_ylim((-rmax,rmax))
ax.set_zlim((-rmax,rmax))
ax.text(-5, 0, -8.5, f'The relative change in energy is: {energychange:.2f} %   dt = {dt:.1E}')
plt.savefig('NBody-output.jpg')
#plt.show()

	
