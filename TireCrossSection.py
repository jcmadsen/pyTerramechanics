# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 08:37:03 2013

A class for parametric generation of pneumatic tire cross sections, based on
a few important parameters. Specify the max diameter and width, then the tire
cross section is generated as a set of two swept arcs. First, the arc of the bottom
section of the tire, OA, is swept a given arc length, using a curvatuve specified by r_OA.
Section AB is generated by specifying a second radius of curvature, r_AB, OR, an xy offset
that defines the center of the circle that will be used to sweep AB, r_sw_XYoff

@author: Justin Madsen
"""
import pylab as py
import matplotlib.pyplot as plt
# for changing legend text size
import matplotlib.font_manager as fm

class TireCrossSection:
    '''
    @class: to investigate different tire cross sectional geometry and its effect
            on the JTM terrain database response, specifically for comparing my results
            with those shown in the literature (experimental data)
    '''
    def __init__(self, dia=56.5, width=20.5, r_OA=56.5/2., r_AB=20.5/2., ds = 0.5,
                 sw_sweep_ang = 3.14159/2.0, r_sw_XYoff=[-1,-1] ):
        '''
        Input:
            dia = tire diameter [in.]
            width = tire width [in.]
            r_Base  = radius of circle that sweeps OA
            r_sw    = radius of circle that sweeps AB
            ds  = increment between nodes for each arc segment
            r_sw_XYoff(optional): the x/y offsets for the SW swept arc. THESE
                    OVERRIDE r_sw, naturally.
        '''
        self._r_tire = dia/2.0
        self._w_tire = width
        self._r_OA = r_OA
        self._ds = ds
        # arrays that contain the XY data for the tire cross sections
        self._xOA = [0,0]
        self._yOA = [0,0]
        self._xAB = [0,0]
        self._yAB = [0,0]
        self._OA_nodes = 0
        self._AB_nodes = 0
        # some aux. data FOR plotting
        self._sw_sweep_ang = sw_sweep_ang
        self._phi0 = -1
        self._th_incr = -1
        self._phi_incr = -1
        
        # if we specify either an x OR y offset for the SW swept circle, set
        # both values as nonzero. NOTE: you'll want to set both values to nonzero
        if( r_sw_XYoff[0] != -1 or r_sw_XYoff[1] != -1):
            self._useABoffset = True
            if( r_sw_XYoff[0] != -1):
                self._AB_x_off = r_sw_XYoff[0]
            else:
                self._AB_x_off = 0.0
            if( r_sw_XYoff[1] != -1):
                self._AB_y_off = r_sw_XYoff[1]
            else:
                self._AB_y_off = 0.0
            # since we specificed the center of the SW circle, we'll calculate
            # the radius of AB when we find the crossSectionNodes
            self._r_AB = 0.0
        else:
            # otherwise, use the input radius, and find the x,y offsets for the
            # SW circle from that.
            self._useABoffset = False
            self._r_AB = r_AB          
            
        # now, calculate the tire nodal locations in the x-y plane
        self.calcCrossSectionNodes()
        
    def reset(self, d=-1,w=-1,r1=-1,r2=-1,ds=-1):
    #_r_tire, w = self._w_tire, r1 = self._r_OA, r2 = self._r_AB, ds = self._ds):
        '''
        Input: new values for the tire cross section, if they are specificed
        '''
        if( d != -1):
            self._r_tire = d/2.
        if( w != -1):    
            self._w_tire = w
        if( r1 != -1):            
            self._r_OA = r1
        if( r2 != -1):    
            self._r_AB = r2
        if( ds != -1):
            self._ds = ds

    def calcCrossSectionNodes(self):
        '''
        @brief: find the (x,y) nodal locations in the tire cross section, where
            x is the longitudinal (forward), y is the lateral dirs.
        Appends:
            _OA_nodes: number of nodes in arc OA
            _th_incr: last rotation along arc OA
            _r_AB:  if set XYoffsets for AB, we'll be finding the radius,
            _phi0: the angle of the last node in the swept arc OA
            _AB_nodes: number of nodes in arc AB
            _phi_incr: increment of the sweep angle along arc AB
            _xOA, _yOA, _xAB, _yAB: the nodal locations
            
        '''
        w = self._w_tire
        ys = self._ds
        
        # *** OA ***
        OA_ymax = w/2.0  # how far in the y-dir does section OA extend?
        OA_rad = self._r_OA          # radius for OA curve
        OA_x_offset = self._r_tire - OA_rad    # if we change the radius of the circle swept, need
                                    # to move center of circle when finding coords
        OA_nodes = int(OA_ymax / ys) + 1
        self._OA_nodes = OA_nodes
        self._xOA = py.zeros((2,OA_nodes))
        self._yOA = py.zeros((2,OA_nodes))
        th = th_old = th_incr = 0.0    # rotation angle
        # just move each node over ys in the x-dir
        for i in range(0, OA_nodes):
            self._yOA[0,i] = ys * i     # right side of curve
            self._yOA[1,i] = -ys * i    # left side of curve
            th_old = th
            th = py.arcsin(ys*i / OA_rad)   # find theta based on y_i
            th_incr = th - th_old       # theta increment, this step
            self._xOA[0,i] = self._xOA[1,i] = OA_x_offset + OA_rad*py.cos(th) # x is also on the circle
        
        # keep the last th_incr, for plotting
        self._th_incr = th_incr
        # *** AB ***
        # this is how I usually do it, but use the input XY offset for AB instead
        if(self._useABoffset):
            # AB_y_off = w/2.0 / 2.0  # essentially the CM of circle that defines AB, use as offsets
            AB_y_off = self._AB_y_off            
            # AB_x_off = r - w/2.0
            AB_x_off = self._AB_x_off
            AB_rad = py.sqrt((self._yOA[0,OA_nodes-1]-AB_y_off)**2 + (self._xOA[0,OA_nodes-1]-AB_x_off)**2)  # new radius of curve to sweep, AB
            self._r_AB = AB_rad     # set the AB_rad, now that we know it
            phi0 = py.arcsin((self._yOA[0,OA_nodes-1]-AB_y_off) / AB_rad)
            self._phi0 = phi0
            arc_len = AB_rad*(self._sw_sweep_ang - phi0)
            AB_nodes = int( arc_len / ys )
            self._AB_nodes = AB_nodes
            self._xAB = py.zeros((2,AB_nodes))
            self._yAB = py.zeros((2,AB_nodes))
            # angle increment based on sweeping about 20% past 90 degrees from x-plane
            phi_incr = (self._sw_sweep_ang - phi0) / AB_nodes  
            self._phi_incr = phi_incr
            for j in range(0,AB_nodes):
                # don't re-do first node on arc, we're starting from last swept OA node
                phi_j = phi0 + (1+j)*phi_incr
                self._yAB[0,j] = AB_y_off + AB_rad * py.sin(phi_j)
                self._yAB[1,j] = -AB_y_off - AB_rad * py.sin(phi_j)
                self._xAB[0,j] = self._xAB[1,j] = AB_x_off + AB_rad * py.cos(phi_j)
        # I'll figure this out later
        else:
            arg = 2
            
            
    def plotCrossSection(self,withArrows=True,withNotes=True):
        '''
        @brief: plot the tire cross section, with or w/out figure labels
        Usage: myTire = TireCrossSection([15.0,5.0])
        myTire.plotCrossSection(withArrows=False)
        '''
        if( self._OA_nodes == 0 or self._AB_nodes == 0):
            print 'cant plot a crossSection when there are no Nodes!\n'
            return
            
        else:
            # OK to plot, set some vals. that will be used
            ys = self._ds
            OA_x_off = self._r_tire - self._r_OA
            AB_x_off = self._AB_x_off 
            AB_y_off = self._AB_y_off
            OA_nodes = self._OA_nodes
            AB_nodes = self._AB_nodes
            phi_incr = self._phi_incr
            r_tire = self._r_tire
            AB_rad = self._r_AB
            th_incr = self._th_incr * (self._r_OA / r_tire)
            # now do the plotting
            prop = fm.FontProperties(size=18)
            fig = plt.figure()
            ax = plt.subplot(111,aspect='equal')
            fig.canvas.set_window_title('Tire Cross section, s_bar = ' + str(ys))
            # plot OA
            ax.plot(self._yOA[0,:],self._xOA[0,:],'k-*',self._yOA[1,:],self._xOA[1,:],'k-*')
            # plot AB
            ax.plot( self._yAB[0,:],self._xAB[0,:], 'r-*',self._yAB[1,:],self._xAB[1,:],'r-*')
            
            
            # plot r_OA
            if(OA_x_off < 0):
                arrow_y0 = (self._yOA[0,OA_nodes/2]/(self._xOA[0,OA_nodes/2]-OA_x_off))*abs(OA_x_off)
            else:
                arrow_y0 = OA_x_off
            # ax.plot( (0.0,self._yOA[1,0]),(0.0,self._xOA[0,0]),'b--o',linewidth=1.5)
            # ax.arrow(arrow_y0,0.0, self._yOA[0,OA_nodes/2]-arrow_y0,self._xOA[0,self._OA_nodes/2],
            #          head_width = 0.5, head_length = 0.5, fc='b',ec='b',length_includes_head=True) 
            ax.plot( (0.0,self._yOA[1,0]),(OA_x_off,self._xOA[0,0]),'b--o',linewidth=1.5)
            ax.arrow(0.0, OA_x_off, self._yOA[0,OA_nodes/2],self._xOA[0,self._OA_nodes/2]-OA_x_off,
                      head_width = 0.5, head_length = 0.5, fc='b',ec='b',length_includes_head=True)

            #plot r_AB
            ax.plot( (AB_y_off,self._yOA[0,OA_nodes-1]),(AB_x_off,self._xOA[0,OA_nodes-1]),'r--o',linewidth=1.5)
            ax.arrow( AB_y_off,AB_x_off, self._yAB[0,AB_nodes/2]-AB_y_off, self._xAB[0,AB_nodes/2]-AB_x_off,
                     head_width = 0.5, head_length = 0.5, fc='r',ec='r',length_includes_head=True)
            # plot the spin axis
            ax.plot( -self._w_tire/2. + self._w_tire*py.arange(0,30)/29., py.zeros(30),'k--',linewidth=1.5)
            leg = ax.legend(('sec. OA','sec. OA','AB','AB',
                             r'$\vec r_\overline{OA}$',r'$\vec r_\overline{AB}$',
                            'wheel axis'),
                            loc=3,prop=prop)
            leg.draggable()
            
            # plot the angles, annotate them
            ax.plot( r_tire*py.sin(th_incr*py.arange(0,OA_nodes/2-1)), r_tire*(py.cos(th_incr*py.arange(0,OA_nodes/2-1)))-r_tire/2.0 )
            ax.annotate(r'$\theta$', xy=(0.5, r_tire/2.0), xytext=(-25,30),size=22,textcoords='offset points', arrowprops=dict(arrowstyle="-|>",color='black'))
            AB_dx = self._xOA[0,OA_nodes-1] - AB_x_off
            AB_dy = self._yOA[0,OA_nodes-1] - AB_y_off
            phi0_idx = (int)((py.pi/2.0 - py.arctan(AB_dx/AB_dy) ) / phi_incr)+1;
            py.plot( AB_y_off + (AB_rad/2.0)*py.sin(phi_incr*py.arange(phi0_idx,AB_nodes)), AB_x_off + (AB_rad/2.0)*(py.cos(phi_incr*py.arange(phi0_idx,AB_nodes))),'r--',linewidth=2. )

            ax.annotate(r'$\phi$', xy=(AB_y_off+AB_rad/3.0, AB_x_off + AB_rad/3.0), xytext=(-5,-55),size=22,textcoords='offset points', arrowprops=dict(arrowstyle="-|>",color='black'))
            
            ax.annotate('O= '+str([r_tire,0.0]), xy=(0.0,r_tire), xytext=(25,45),size=16,textcoords='offset points', arrowprops=dict(arrowstyle="-|>",color='black'))
            ax.annotate("A= [%.2f" % self._xOA[0,OA_nodes-1] + ", %.2f" % self._yOA[0,OA_nodes-1] +']', xy=(self._yOA[0,OA_nodes-1],self._xOA[0,OA_nodes-1]),
                        xytext=(60,25),size=16,textcoords='offset points', arrowprops=dict(arrowstyle="-|>",color='black'))
            ax.annotate("B= [%.2f" % self._xAB[0,AB_nodes-1]+ ", %.2f" % self._yAB[0,AB_nodes-1] +']', xy=(self._yAB[0,AB_nodes-1],self._xAB[0,AB_nodes-1]),
                        xytext=(40,35),size=16,textcoords='offset points', arrowprops=dict(arrowstyle="-|>",color='black'))
            py.xlabel('y-dir [inches]',size=14)
            py.ylabel('x-dir [inches]',size=14)
            py.grid('on')
    
    def summary(self):
        print 'nodes in each cross section:\n'
        print 'r_OA, r_AB = ' + str(self._r_OA) +', ' + str(self._r_AB) +'\n'
        print 'OA_nodes = ' + str(self._OA_nodes) +'\n'
        print 'AB_nodes = ' + str(self._AB_nodes) +'\n'
        
if __name__ == '__main__':
    # r = 56.5/2.0
    # sw_xy off = py.array([r-w/2.0, w/4.0])
    # w = 20.5
        
    # Fig 4.3, PhD doc
    r = 37.0/2.0    # physical dia
    w = 12.5        # physical width
    incr = 0.5      # arc length increment in cross section arcs
    sw_xy_off = [r-w/3.0, w/3.0]
    myTire = TireCrossSection(dia =2*r, width =w, ds =incr, sw_sweep_ang= py.pi/2.0+0.45,
                              r_sw_XYoff = sw_xy_off, r_OA=56.5/2.)
                              
    
    '''
    # Fig2
    myTire = TireCrossSection(ds=1.0,sw_sweep_ang=py.pi/2.0+0.2,r_sw_XYoff=[r-w/3.0, w/2.4],
                              r_OA=565./2.)
    '''
    # other Figs
    
    
    
    # plot the cross section, the summary, and show the plot
    myTire.plotCrossSection()
    myTire.summary()

    plt.show() 
    print 'done'