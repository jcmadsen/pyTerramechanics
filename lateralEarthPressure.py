# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 20:40:14 2012
@brief lateral earth pressure in a homogenous soil that
 causes incipient plastic failure (e.g., the soil strength
 to bulldozing)
 
Assuming constant force along width

Fb_max = gamma*Z^2*nPhi/2 + 2*c*Z*nPhi^0.5
    -gamma and c are soil material properties, won't vary a lot
    -sinkage and width will vary the most
    -phi, and thus nPhi will vary, but less

@lit Bekker, "Theory of Land Locomotion" (1956)
@author: Justin Madsen, 2013
"""
import pylab as py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def F_per_len(z,nPhi,c,gamma):
    """ multiply by width applied to get total force at incipient failure,
    assuming force constant with width.
    Input:
        z:  lug height
        nPhi: flow value, tan(pi/4+ phi/2)^2, phi = soil internal friction 
        c:  soil cohesion
        gamma: soil unit weight
    
    """
    return (gamma*nPhi*z**2)/2.0 + 2.*c*z*py.sqrt(nPhi)

def F_per_len_sur(z,nPhi,c,gamma,sur):
    """ same as above, but if there is a surcharge/pile of extra soil in blade
    Input:
        sur: average depth of surcharge on the surface of the bulldozing wedge
        len: longitudinal length of bulldozing wedge. Forms an angle of 
                pi/4- phi/2 with the horizontal
    """
    return (gamma*nPhi*z**2)/2.0 + 2.*c*py.sqrt(nPhi)*z + nPhi*gamma*sur*z
    
def F_blade(z,nPhi,c,gamma,sur,beta):
    """ a blade with inclination angle from the vertical of beta
    Input:
        beta: angle from vertical
    Return:
        [Fp, Fn, Ft] the total, normal and tangent components of the blade force
    """
    return (gamma*nPhi*z**2)/2.0 + 2.*c*py.sqrt(nPhi)*z + nPhi*gamma*sur*z

def Ft_blade(z,nPhi,c,gamma,sur,beta)

# default font size
font = {'size' : 15}
plt.rc('font', **font)

# n_phi = tan(pi/4 + [0,pi/8])^2
nphi = py.tan(py.pi/4.+py.arange(0,py.pi/8.,py.pi/64.))**2
py.plot(py.arange(0,py.pi/4,py.pi/32.),nphi)
py.xlabel(r'$\phi$ [rad]')
py.ylabel(r'N$\phi$')
py.grid(True)

# soil properties
gamma = 100./(12.**3.)  # soil unit weight [lb/in3]
c = 2.9 # cohesion, [psi]
# variables to plot. b and Z can vary a lot, from 0 to ~1 foot
# range of phi, nphi is smaller
b = py.arange(1.0,10.0,0.2)    # blade width [in]
Z = py.arange(1.0,10.0,0.2)    # blade depth [in]
phi = py.arange(0,py.pi/4.,py.pi/128.)   # internal friction angle
nphi = py.tan(py.pi/4.+phi/2.0)**2  # flow value

b0 = 4.0    # constant width [in]
phi0 = py.pi/8.0    # [-]   constant phi = pi/8 = 22.3 [deg]
nphi0 = py.tan(py.pi/4.+ phi0/2.0 )**2    # constant flow value
z0 = 4.0    # constant Z[in]

# q = [phi, b, z]
# plot the bulldozing force while varying f(phi, b) for constant Z
nphiMG, bMG = py.meshgrid(nphi,b)
q12 = bMG*F_per_len(z0,nphiMG,c,gamma)
# (gamma*nphiMG*z0**2/2.+2.*c*z0*py.sqrt(nphiMG))
fig = plt.figure()
#---- First subplot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# I'd actually rather see phi, b as the axis variables, since nPhi = f(phi)
phiMG, bMG = py.meshgrid(phi,b)
surf = ax.plot_surface(phiMG, bMG, q12, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)                 
# ax.set_zlim3d(-1.01, 1.01)
ax.set_xlabel(r'$\phi$ [rad]')
ax.set_ylabel('b [in]')
ax.set_zlabel('F [lb]')
gamma_str = "%4.3f" % gamma
ax.set_title(r'$\gamma$ = '+ gamma_str +'[lb/in3], z = ' +str(z0) + '[in]')
fig.colorbar(surf, shrink=0.5, aspect=10)

# plot bulldozing force while varying f(b, z) for constant phi
bMG, zMG = py.meshgrid(b, Z)
q23 = bMG*F_per_len(zMG, nphi0, c, gamma)
# (gamma*nphi0*zMG**2 / 2.0 + 2.0*c*zMG*py.sqrt(nphi0))
fig2 = plt.figure()
#---- second
ax2 = fig2.add_subplot(1,1,1, projection='3d')
surf2 = ax2.plot_surface(bMG, zMG, q23, rstride=1, cstride=1, cmap=cm.cool,
                       linewidth=0, antialiased=False)
ax2.set_xlabel('b [in]')
ax2.set_ylabel('z [in]')
ax2.set_zlabel('F [lb]')
phi_str = "%4.3f" % phi0
ax2.set_title(r'$\gamma$ = '+ gamma_str +'[lb/in3], $\phi$ = ' + phi_str + '[in]')
fig2.colorbar(surf2, shrink=0.5, aspect=10)

# F_b = f(phi, z)
nphiMG, zMG = py.meshgrid(nphi,Z)
q13 = b0 * F_per_len(zMG, nphiMG, c, gamma)
# (gamma*nphiMG*zMG**2 /2.0 + 2.0*c*zMG*py.sqrt(nphiMG)) 
# --- third plot
fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1, projection='3d')
# again, I'd rather plot phi directly rather than nphi
phiMG, zMG = py.meshgrid(phi, Z)
surf3 = ax3.plot_surface(phiMG, zMG, q13, rstride=1, cstride=1, cmap=cm.seismic,
                         linewidth = 0, antialiased=False)
ax3.set_xlabel(r'$\phi$ [rad]')
ax3.set_ylabel('z [in]')
ax3.set_zlabel('F [lb]')
ax3.set_title(r'$\gamma$ = '+ gamma_str +'[lb/in3], b = ' + str(b0) + '[in]')
fig3.colorbar(surf3, shrink=0.5, aspect=10)


plt.show()