# -*- coding: utf-8 -*-
"""
Plot mohr's static friction strength, vs. mohr coulomb soil strength.
Take the minimum of the two, as the slip failure point transitions from 
the tire/soil surface to the soil. 

Also, take slip rate into account with a term, acts as a transition between the
case where shear strength is a static term (mu*pn) or determined by the soil shear strength
kap = 1-py.exp(-s_dot/ks)
then tau_max = kap*tau_soil + (1-kap)*tau_static

Created on Fri Apr 12 11:34:39 2013

@author: Justin Madsen
"""

import pylab as py
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import logging as lg

# tau_max, soil
def tau_max_soil(p,c,phi):
    return c + p*py.tan(phi)

# just based on soil strength
def tau_soil(p,c,phi,j,k):
    t_m = tau_max_soil(p,c,phi)
    tS = 1 - py.exp(-j/k)
    t_out = t_m*tS
    return t_out
    
# just based on static friction
def tau_mu(p,mu,j,k):
    tS = 1 - py.exp(j/k)
    t_out = p*mu*tS
    return t_out

# From harnisch, return the min. of the two shear values
def tau_a(p,mu,c,phi,j,k):
    tS = 1 -  py.exp(-j/k)
    t1 = c*p*py.tan(phi)
    t2 = mu*p
    t_out = min(t1,t2)*tS
    return t_out

# from my PhD thesis, scales "in" lug/treads w.r.t.
# slip velocity, resulting in a kappa in the range from 0-1
# 0 = no soil interlock, 1 = complete soil interlock
def tau_b(p,mu,c,phi,j,k,kappa):
    tS = 1 - py.exp(-j/k)
    t1 = c*p*py.tan(phi)
    t2 = mu*p
    t_out = (kappa*t1 + (1-kappa)*t2)*tS
    return t_out
    
def plotshearMixed(mu,c,phi,pcrit):
    # @brief: plot the shear for a tire/soil combination
    pRange1 = py.arange(0,pcrit, pcrit/100.)
    pRange2 = py.arange(pcrit,2.*pcrit, pcrit/100.)
    # calculate both shears for both ranges
    
    tau_s_1 = tau_max_soil(pRange1,c,phi)
    tau_s_2 = tau_max_soil(pRange2,c,phi)
    tau_mu_1 = pRange1*mu
    tau_mu_2 = pRange2*mu
    
    fig = plt.figure()
    ax = plt.subplot(111)
    fig.canvas.set_window_title('tau_max')

    ax.plot(pRange1,tau_s_1,'b--',pRange2,tau_s_2,'b-',linewidth=2.5,)
    ax.plot(pRange1,tau_mu_1,'r-',pRange2,tau_mu_2,'r--',linewidth=2.5)
    ax.legend((r'$\tau,soil $',r'$ \tau,soil $',r'$ \tau_\mu $',r'$\tau_\mu $'), loc=2)
    py.xlabel(r'$ \sigma [psi]$',size=20)
    py.ylabel(r'$ \tau,max [psi]$',size=20)
    py.grid('on')
    py.title(r'$\mu =  %.2f $' %mu + ', c= %.2f'%c + r', $\phi =$ %.2f'%phi)
    

def plot_kappa(ks,s_dot_max):
    # @brief kappa scales tau_max based on slip velocity, constant ks
    s_dot = py.arange(0,s_dot_max,s_dot_max/100.)
    kap = 1-py.exp(-s_dot/ks)
    fig = plt.figure()
    ax = plt.subplot(111)    
    ax.plot(s_dot,kap,linewidth=1.5)
    py.xlabel('s_tb')
    py.ylabel(r'$\kappa $')
    py.grid('on')
    py.title('mod. max shear for slip rate, k = %.2f'%ks)
    
def plot_jp_tmax_surf(mu,c,phi,pmax,smax,ks):
    # @brief tau max based on varying slip rate, normal pressure
    s_dot = py.arange(0,smax,smax/100.)
    prange = py.arange(0,pmax,pmax/100.)
    kap = 1-py.exp(-s_dot/ks)   #  kappa
    
    TMAX = py.zeros((len(kap),len(prange)))
    tphi = py.tan(phi)  # keep tan(phi) handy
    for k_i in range(0,len(kap)):
        k_tmp = kap[k_i]    
        for p_j in range(0,len(prange)):
            p_tmp = prange[p_j]
            TMAX[k_i][p_j] = k_tmp*(c+p_tmp*tphi) + (1-k_tmp)*p_tmp*mu
            
    fig = plt.figure()           
    ax = fig.add_subplot(121)
    # should be ok to plot the surface
    S, P = py.meshgrid(s_dot, prange)
    CS = plt.contour(S,P,TMAX,8,colors='k',linewidths=1.5)
    plt.clabel(CS,inlne=1,fontsize=16)
    img = plt.imshow(TMAX, interpolation='bilinear', origin='lower',
                     cmap=cm.jet,extent=(min(s_dot),max(s_dot),min(prange),max(prange)))
    CBI = plt.colorbar(img, orientation='vertical',shrink=0.8)
    CBI.set_label(r'$\tau ,max $[psi]')
    ax.set_title(r'$\tau ,max = f(\sigma,\kappa), ks=%.2f $'%ks)
    ax.set_xlabel('slip rate [in/sec]')
    ax.set_ylabel(r'$\sigma_z $',size=24)
    
    # use twice ks, re-calc what's necessary, then replot
    ks2 = ks * 2
    kap2 = 1-py.exp(-s_dot/ks2)
    
    TMAX2 = py.zeros((len(kap2),len(prange)))
    # tphi = py.tan(phi)  # keep tan(phi) handy
    for k_i in range(0,len(kap2)):
        k2_tmp = kap2[k_i]    
        for p_j in range(0,len(prange)):
            p_tmp = prange[p_j]
            TMAX2[k_i][p_j] = k2_tmp*(c+p_tmp*tphi) + (1-k2_tmp)*p_tmp*mu
            
    
    #fig = plt.figure()           
    ax = fig.add_subplot(122)
    # should be ok to plot the surface
    # S, P = py.meshgrid(s_dot, prange)
    CS2 = plt.contour(S,P,TMAX2,8,colors='k',linewidths=1.5)
    plt.clabel(CS2,inlne=1,fontsize=16)
    img2 = plt.imshow(TMAX2, interpolation='bilinear', origin='lower',
                     cmap=cm.jet,extent=(min(s_dot),max(s_dot),min(prange),max(prange)))
    CBI2 = plt.colorbar(img2, orientation='vertical',shrink=0.8)
    CBI2.set_label(r'$\tau ,max $[psi]')
    ax.set_title(r'$\tau ,max = f(\sigma,\kappa), ks=%.2f $'%ks2)
    ax.set_xlabel('slip rate [in/sec]')
    ax.set_ylabel(r'$\sigma_z $',size=24)
    
if __name__ == '__main__':
 # logger
    lg.basicConfig(fileName = 'logFile.log', level=lg.WARN, format='%(message)s')
    # default font size
    font = {'size' : 18}
    matplotlib.rc('font', **font)
    
    ## params
    mu = 0.605    # static friction coef.
    c = 3.0 # cohesion, psi
    phiRad = 17.0 * py.pi/180.    # deg. to rad
    max_pressure = 10.0 # max vertical pressure, psi
    plotshearMixed(mu,c,phiRad,max_pressure)
    
    # now plot tau_max, function of slip rate, normal pressure
    slip_rateMax = 5.0  # in/sec
    ks = 1.1    # approx. inverse slope of kappa modifier
    plot_jp_tmax_surf(mu,c,phiRad,max_pressure,slip_rateMax,ks)    
    
    py.show()