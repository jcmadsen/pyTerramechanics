# -*- coding: utf-8 -*-
"""
Created on Tue May 01 22:40:23 2012

@author: JustinFast

This runs tests to replicate some tests on Boussinesq, Cerruti, etc.
that can be found in the literature
"""

import pylab as py
import logging as lg
import BousCerr as BC
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colorbar


if __name__ == '__main__':
    lg.basicConfig(fileName = 'logFile.log',
                   level=lg.INFO,
                   format='%(message)s')
   
   

# Wong (2001) Theory of Ground Vehicles
#   Fig 2.10
#   Compare analytical values of Boussinesq 
    
    contour3D = False   # create a 3-D contour plot? (2D is plotted by default)
    # 3 = metal, 4 = firm, 5 = medium, 6 = soft
    nuVect = py.array([3.0, 4.0, 5.0, 6.0])
    Fvect = py.array([1650./4.0, 0.0, 1650])
    platesize = 137.5             # p0 = Fz/platesize
    # minimum sum_stress to consider the contribution of cerruti as a percent of the sum
    min_sum_stress_cutoff = 2.0
    # in this case, p0 = 12 psi
    # set the plotting limits. Think of [r,z] as [x,y], ...
    rlim = py.array([0.2,12.])
    zlim = py.array([1.0,12.0]) # z is positive down in VTTI
    rz_rez = 0.2                # grid dim of SMASH
    R = py.arange(rlim[0],rlim[1]+rz_rez,rz_rez)
    Z = py.arange(zlim[0],zlim[1]+rz_rez,rz_rez)
    Rm,Zm = py.meshgrid(R,Z)    # get a meshgrid rep of R,Z          
    angDeg = 45 * py.pi/180.0   # cone angle off the horizontal, in degrees
    slope = py.tan(angDeg)
    for i in range(0, len(nuVect)):
        # current value of Frolich paramter, nu
        nu = nuVect[i]
        # Compute sigma_z for boussinesq and Cerruti
        BousMat = BC.get_BousPlateMat(Fvect,R,Z,nu,platesize,slope)
        # Question: is sigma_z negative when cos(theta) for the applied
        #   surface force, H, is > pi [rad] ???
        CerrMat = BC.get_CerrMat(Fvect,R,Z,platesize,slope)
        # figure out what sign Cerr. is to be, and superimpose
        # with bous
        superPosMat = BousMat + CerrMat
        # cerruti as a percent of the total stress. Note, avoid entries
        # with zero total stress
        cerrPercentMat = py.zeros((len(Z),len(R)))
        for row in range(0,len(Z)):
            for col in range(0,len(R)):
                curr_sum = superPosMat[row][col]
                # arbitrary limit to compute cerruti percent of total stress
                if(curr_sum > min_sum_stress_cutoff):
                    cerrPercentMat[row][col] = (CerrMat[row][col] / curr_sum) * 100.0
                else:
                    cerrPercentMat[row][col] = 0.0
        
        # 1) Either create a 3-d contour plot...
        if( contour3D ):
            # create new figures for each of Boussinesq, cerruti
            fig_b = plt.figure()
            ax_b = fig_b.gca(projection='3d')
            fig_c = plt.figure()
            ax_c = fig_c.gca(projection='3d')
            fig_sum = plt.figure()
            ax_sum = fig_sum.gca(projection='3d')
            # 3d contour surface for boussinesq            
            surf_b = ax_b.plot_surface(Rm,Zm,BousMat, cmap=cm.RdBu,rstride=1,cstride=1,linewidth=0,antialiased=False)   
            # 3d contour surface for cerruti            
            surf_c = ax_c.plot_surface(Rm,Zm,CerrMat, cmap=cm.RdBu,rstride=1,cstride=1,linewidth=0,antialiased=False)
            surf_sum = ax_sum.plot_surface(Rm,Zm,superPosMat, cmap=cm.RdBu,rstride=1,cstride=1,linewidth=0,antialiased=False)            
            # colorbars
            cbar_b = fig_b.colorbar(surf_b,shrink=0.5,aspect=5)
            cbar_c = fig_c.colorbar(surf_c,shrink=0.5,aspect=5)
            cbar_sum = fig_sum.colorbar(surf_sum,shrink=0.5,aspect=5)
        else:
            # 2) ... or use a 2-d representation,
            # BOUS.
            fig_b = plt.figure()
            fig_b.canvas.set_window_title('Boussinesq, nu ='+str(nu))
            ax_b = fig_b.add_subplot(111)
            # Zm = -Zm
            CS_b = plt.contour(Rm,Zm,BousMat)
            plt.clabel(CS_b, inline=1,fontsize=10)
            # if you want a colorbar added:
            plt.flag()
            img_b = plt.imshow(BousMat, interpolation = 'bilinear', origin='lower',
                             cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))                
            CBI_b = plt.colorbar(img_b, orientation='horizontal', shrink = 0.5)            
            # plot stresses downwards
            ax_b.invert_yaxis()
            
            # Cerruti
            fig_c = plt.figure()
            fig_c.canvas.set_window_title('Cerruti, nu ='+str(nu))
            ax_c = fig_c.add_subplot(111)
            CS_c = plt.contour(Rm,Zm,CerrMat)
            plt.clabel(CS_c, inline=1,fontsize=10)
            img_c = plt.imshow(CerrMat, interpolation = 'bilinear', origin='lower',
                 cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))             
            
            CBI_c = plt.colorbar(img_c, orientation='horizontal', shrink = 0.5)
            ax_c.invert_yaxis()
            
            # sum of Cerruti and Bous
            fig_sum = plt.figure()
            fig_sum.canvas.set_window_title('sum, nu = ' +str(nu))
            ax_sum = fig_sum.add_subplot(111)
            CS_sum = plt.contour(Rm,Zm,superPosMat)
            plt.clabel(CS_sum, inline=1,fontsize=10)
            img_sum = plt.imshow(superPosMat, interpolation = 'bilinear', origin='lower',
                 cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))             
            
            CBI_sum = plt.colorbar(img_sum, orientation='horizontal', shrink = 0.5)
            ax_sum.invert_yaxis() 

            # cerruti as a percentage of the total stress
            # only consider stresses whose sum is greater than a minimum cutoff
            fig_cper = plt.figure()
            fig_cper.canvas.set_window_title('cerruti %, nu = '+str(nu))
            ax_cper = fig_cper.add_subplot(111)
            CS_cper = plt.contour(Rm,Zm,cerrPercentMat)
            plt.clabel(CS_cper,inline=1,fontsize=10)
            img_cper = plt.imshow(cerrPercentMat, interpolation = 'bilinear', origin='lower',
                                  cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
            CBI_cper = plt.colorbar(img_cper, orientation='horizontal', shrink = 0.5)
            ax_cper.invert_yaxis()                     
                        
        # set the titles and axes for Bousinesq
        ax_b.set_title('Bous, nu = ' + str(nu))
        ax_b.set_xlabel('R [in]')
        ax_b.set_ylabel('Z [in]')
        # for Cerruti also
        ax_c.set_title('Cerr, nu = ' + str(nu))
        ax_c.set_xlabel('R [in]')
        ax_c.set_ylabel('Z [in]')
        # for the superposition
        ax_sum.set_title('sum, nu = ' + str(nu))
        ax_sum.set_xlabel('R [in]')
        ax_sum.set_ylabel('Z [in]')
        # for cerruti as a percentage of the sum stress
        ax_cper.set_title('cerruti %, nu = ' + str(nu) + ', cutoff ' + str(min_sum_stress_cutoff))
        ax_cper.set_xlabel('R [in]')
        ax_cper.set_ylabel('Z [in]')

    '''
    
# Compare the plate and point load versions of Boussinesq, using same input
#   values from the Wong (2001) Fig 2.10 example
    # 3 = metal, 4 = firm, 5 = medium, 6 = soft
    nu = 5
    
    Fvect = py.array([21.2, 0.0, 80.0])
    platesize = 4.0  # p0 = Fz/platesize
    # in this case, p0 = 12 psi
    # set the plotting limits. Think of [r,z] as [x,y], ...
    rlim = py.array([-4.0,4.0])
    zlim = py.array([0.0,6.0]) # z is positive down in VTTI
    rz_rez = 0.05                # grid dim of SMASH
    Rvect = py.arange(rlim[0],rlim[1]+rz_rez,rz_rez)
    Zvect = py.arange(zlim[0],zlim[1]+rz_rez,rz_rez)
    Rm,Zm = py.meshgrid(Rvect,Zvect)    # get a meshgrid rep of R,Z          
    angDeg = 45.0 * py.pi/180.0   # cone angle off the horizontal, in degrees
    slope = 2.0 # slope = py.tan(angDeg)

    # calculate subsoil stresses via Boussinesq (point and plate versions)
    BousMat = BC.get_BousMat(Fvect,Rvect,Zvect,nu,platesize,slope)
    BousPlateMat = BC.get_BousPlateMat(Fvect,Rvect,Zvect,nu,platesize,slope)
    CerrMat = BC.get_CerrMat(Fvect,Rvect,Zvect, platesize, slope)
    BousMatTire = BC.get_BousMatTire(Fvect,Rvect,Zvect,nu, platesize, slope)


    # Plot 1) boussinesq, 2) cerruti, 3) combination
     # 1) BOUSSINESQ
    fig = plt.figure()
    ax = fig.add_subplot(131)
    levels = [0.5, 1.0, 2.0, 5.0, 10.0, 15.]
    CS = plt.contour(Rm,Zm,BousPlateMat,levels,colors='k',linewidths=1.5)
    plt.clabel(CS,inline=1,fontsize=18)
    # plt.flag()
    img = plt.imshow(BousPlateMat, interpolation = 'bilinear', origin='lower',
                     cmap=cm.jet,extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.6,ticks=levels)
    CBI.set_label(r'$\sigma_z$ [psi]',fontsize=20)
    # set the titles and axes
    ax.set_title(' Boussinesq, P(y) = ' + str(Fvect[2]/platesize) + ' psi',fontsize=18)
    ax.set_xlabel('R [in]',fontsize=16)
    ax.set_ylabel('Z [in]',fontsize=16)
     

    ax = fig.add_subplot(132)
    CS = plt.contour(Rm,Zm,BousPlateMat,5,colors='k',linewidths=1.5)
    plt.clabel(CS,inline=1,fontsize=18)
    # plt.flag()
    img = plt.imshow(BousPlateMat, interpolation = 'bilinear', origin='lower',
                     cmap=cm.jet,extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.4,ticks=[0.,2.,4.,6.,8.,10.])
    CBI.set_label(r'$\sigma_z$ [psi]',fontsize=20)
    # set the titles and axes
    ax.set_title(r'$ Boussinesq Plate, \sigma_z = $' + str(Fvect[2]/platesize) + ' psi',fontsize=18)
    ax.set_xlabel('R [in]',fontsize=16)
    ax.set_ylabel('Z [in]',fontsize=16)    
    
    ax = fig.add_subplot(133)
    CS = plt.contour(Rm,Zm,BousMatTire,5,colors='k',linewidths=1.5)
    plt.clabel(CS,inline=1,fontsize=18)
    # plt.flag()
    img = plt.imshow(BousMatTire, interpolation = 'bilinear', origin='lower',
                     cmap=cm.jet,extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.4,ticks=[0.,2.,4.,6.,8.,10.])
    CBI.set_label(r'$\sigma_z$ [psi]',fontsize=20)
    # set the titles and axes
    ax.set_title(r'$ Boussinesq Tire, \sigma_z = $' + str(Fvect[2]/platesize) + ' psi',fontsize=18)
    ax.set_xlabel('R [in]',fontsize=16)
    ax.set_ylabel('Z [in]',fontsize=16)
    
    # 2) Cerruti
    ax = fig.add_subplot(132)
    CS = plt.contour(Rm,Zm,CerrMat,12,colors='k',linewidths=1.5)
    plt.clabel(CS,inline=1,fontsize=16)
    img = plt.imshow(CerrMat, interpolation = 'bilinear', origin='lower',
                     cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8,ticks=[-2.,-1.,0.,1.,2.0])
    CBI.set_label(r'$\sigma_z$ [psi]',fontsize=20)
    # set the titles and axes
    ax.set_title('Cerruti, tau = ' + str(Fvect[0]/platesize) +' psi' )
    ax.set_xlabel('R [in]',fontsize=16)
    ax.set_ylabel('Z [in]',fontsize=16)
    
    # 3) Combined
    ax = fig.add_subplot(133)
    CS = plt.contour(Rm,Zm,BousMat + CerrMat,8,colors='k',linewidths=1.5)
    plt.clabel(CS,inline=1,fontsize=16)
    img = plt.imshow(BousMat + CerrMat, interpolation = 'bilinear', origin='lower',
                     cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    CBI.set_label(r'$\sigma_z$ [psi]',fontsize=20)
    # set the titles and axes
    ax.set_title('Bous + Cerr, nu =  ' + str(nu))
    ax.set_xlabel('R [in]',fontsize=16)
    ax.set_ylabel('Z [in]',fontsize=16)
    
    # Plot the subsoil stresses according to a rigid wheel load at the surface
    # BC.plotSuperPos_withTire(R_vect=Rvect, Z_vect=Zvect, tire_d=37.0, sink=4.0, max_p=12.0, n_pts=5)

    '''
    
    
    '''
    # plot the point force subsoil vertical stresses, BOUSSINESQ
    fig = plt.figure()
    ax = fig.add_subplot(131)
    CS = plt.contour(Rm,Zm,BousMat,10,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    # plt.flag()
    img = plt.imshow(BousMat, interpolation = 'bilinear', origin='lower',
                     cmap=cm.jet,extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    # set the titles and axes
    ax.set_title('Bous, Cone constraint, nu = ' + str(nu))
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]')
    
    # plot the plate subsoil vertical stresses
    ax = fig.add_subplot(132)
    CS = plt.contour(Rm,Zm,BousPlateMat,10,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    img = plt.imshow(BousPlateMat, interpolation = 'bilinear', origin='lower',
                     cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    # set the titles and axes
    ax.set_title('BousPlate, nu = ' + str(nu))
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]')
    
    # plot difference between plate, cone constraint eqs
    ax = fig.add_subplot(133)
    CS = plt.contour(Rm,Zm,BousPlateMat-BousMat,6,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    img = plt.imshow(BousPlateMat-BousMat, interpolation = 'bilinear', origin='lower',
                     cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    # set the titles and axes
    ax.set_title('BousPlate - Bous, nu = ' + str(nu))
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]')

    # plot the point force subsoil vertical stress, CERRUTI
    fig = plt.figure()
    ax = fig.add_subplot(111)
    CS = plt.contour(Rm,Zm,CerrMat,6,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    # plt.flag()
    img = plt.imshow(CerrMat, interpolation = 'bilinear', origin='lower',
                     cmap=cm.jet,extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.5)
    # set the titles and axes
    ax.set_title('Cerruti, conic constraint')
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]')
    
    '''
    
    '''
    # plot the difference between the two
    ax = fig.add_subplot(133)
    stressDiff = BousPlateMat - BousMat
    CS = plt.contour(Rm,Zm,stressDiff,6,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    img = plt.imshow(stressDiff, interpolation='bilinear', origin='lower',
                     cmap=cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    # set the titles and axes
    ax.set_title('BousPlate - Bous, nu = ' + str(nu))
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]') 
    '''
    
    '''
    # Let's see the cumulative effect on the soil response model
    [rhoMat, rhoDeltaMat, delta_zMat, eMat] = test.getRhoSinkageE_Mat(BousMat)
    [rhoMat2, rhoDeltaMat2, delta_zMat2, eMat2] = test.getRhoSinkageE_Mat(BousPlateMat)
    # calculate the differences between the two versions    
    rho_diff = rhoMat2 - rhoMat
    z_diff = delta_zMat2 - delta_zMat
    e_diff = eMat2 - eMat
        
    
    # plot the bulk density difference
    fig = plt.figure()
    ax = fig.add_subplot(131)
    CS = plt.contour(Rm,Zm,rho_diff,6,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    # plt.flag()
    img = plt.imshow(rho_diff, interpolation = 'bilinear', origin='lower',
                     cmap=cm.jet,extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    # set the titles and axes
    ax.set_title('Rho difference, nu = ' + str(nu))
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]')
    
    # plot the plate subsoil vertical stresses
    ax = fig.add_subplot(132)
    CS = plt.contour(Rm,Zm,z_diff,6,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    img = plt.imshow(z_diff, interpolation = 'bilinear', origin='lower',
                     cmap = cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    # set the titles and axes
    ax.set_title('delta_z difference, nu = ' + str(nu))
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]')

    # plot the difference between the two
    ax = fig.add_subplot(133)
    CS = plt.contour(Rm,Zm,e_diff,6,colors='k')
    plt.clabel(CS,inline=1,fontsize=10)
    img = plt.imshow(e_diff, interpolation='bilinear', origin='lower',
                     cmap=cm.jet, extent=(rlim[0],rlim[1],zlim[0],zlim[1]))
    ax.invert_yaxis()
    CBI = plt.colorbar(img, orientation='horizontal', shrink=0.8)
    # set the titles and axes
    ax.set_title('Energy difference, nu = ' + str(nu))
    ax.set_xlabel('R [in]')
    ax.set_ylabel('Z [in]')     
    '''
    
        
    
    py.show()