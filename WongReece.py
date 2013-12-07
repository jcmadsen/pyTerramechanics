# -*- coding: utf-8 -*-
"""
Created on Tue Oct 01 09:20:51 2013

@author: newJustin

Calculate the drawbar pull for a towed wheel, based on 
North Gower claey loam, which was used in the JTM Terramechanics runs
Caluclations for a towed rigid wheel in sand are from Wong/Reece[67]

"""
# I already wrote functions to find pressure/shear from Bekker, Reece
import pylab as py
import scipy.optimize as sci_opt
import scipy.integrate as sci_int
import logging as lg
import matplotlib.pyplot as plt

def degToRad(degrees):
    return degrees* py.pi / 180.0
    
def radToDeg(radians):
    return radians * 180.0 / py.pi

# find the slip displacement, towed wheel, in region AC, fig 2
def j1(th,th0,th1,r,slip):
    Kv = (1.0/(1.0+slip)) * ((1.0+slip)*(py.sin(th1)-py.sin(th0))/(th1-th0) - 1.0)
    slip_out = r*((th1-th)*(1.0+Kv*(1.0+slip)) - (1.0+slip)*(py.sin(th1)-py.sin(th)) )
    return slip_out

# find the slip displacement, towed wheel, in the region AE, fig 2
def j2(th,th0,r,slip):
    slip_out = r*((th0-th) - (1.0+slip)*(py.sin(th0)-py.sin(th)) )
    return slip_out

# slip displacement, driven wheel
def jdriven(th,th1,r,slip):
    j_out = r*( (th1-th) - (1.0-slip)*(py.sin(th1)-py.sin(th)) )
    return j_out
    
# normal stress, front region (AC, fig. 2). Works for both driven and towed case
def sig_1(th,th1,r,b,n,k1,k2):
    sigma_out = ((py.cos(th) - py.cos(th1))**n) *(k1+k2*b)*(r/b)**n
    return sigma_out

# normal stress, bottom region (AE, fig. 2)
# can be used for driven wheels; replace th0 with th_m  
def sig_2(th,th0,th1,th2,r,b,n,k1,k2):
    sigma_out = ((py.cos(th1- (th-th2)*(th1-th0)/(th0-th2)) - py.cos(th1))**n) *(k1+k2*b)*(r/b)**n
    return sigma_out

# towed shear stress, front region (AC)
def tau_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip):
    j_disp = j1(th,th0,th1,r,slip)
    if j_disp < 0:
        j_disp = -j_disp
    tau_out = (c + sig_1(th,th1,r,b,n,k1,k2)*py.tan(phi))*(1.0-py.exp(-j_disp/K) )
    return tau_out    

# towed shear stress, bottom region (AE)    
def tau_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
    j_disp = j2(th,th0,r,slip)
    if j_disp < 0:
        j_disp = -j_disp
    tau_out = (c+ sig_2(th,th0,th1,th2,r,b,n,k1,k2,)*py.tan(phi))*(1.0-py.exp(-j_disp/K) )
    return tau_out

# driven shear stress, front section
def tau_d1(th,th1,r,b,n,k1,k2,c,phi,K,slip):
    j_disp = jdriven(th,th1,r,slip)
    if (j_disp < 0):
        j_disp = -j_disp
    tau_out = (c+sig_1(th,th1,r,b,n,k1,k2) * py.tan(phi) )*(1.0 - py.exp(-j_disp/K))
    return tau_out

# driven shear stress, bottom section
def tau_d2(th,th_m,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
    j_disp = jdriven(th,th1,r,slip)
    if( j_disp < 0):
       j_disp = -j_disp
    tau_out = (c+sig_2(th,th_m,th1,th2,r,b,n,k1,k2)*py.tan(phi))*(1.0-py.exp(-j_disp/K))
    return tau_out
    
# for a driven wheel, find the inflection angle
def theta_m(th1,c1,c2,slip):
    outval = (c1+c2*slip)*th1
    return outval

# Include Bekker's formulas for comparison
def sig_bek(th,th1,r,b,n,k1,k2):
    sig_out = ((py.cos(th) - py.cos(th1))**n) *(k1/b + k2)*(r)**n
    return sig_out

def tau_bek(th,th1,r,b,n,k1,k2,coh,phi,K,slip):
    j_disp = jdriven(th,th1,r,slip)
    tau_out = (coh +sig_bek(th,th1,r,b,n,k1,k2) *py.tan(phi) ) *(1.0 -py.exp(-j_disp /K) )
    return tau_out

class WongReece:
    '''
    A class for running the wong/reece equilibrium model
    '''
    
    def __init__(self,coh,phi,k1,k2,n,K,wheel_r,wheel_b,weight_W,skid,
                 c1=0.43,c2=0.32,gen_plots=False,units='ips'):
        '''
        Purpose:
            requires the user to input the necessary soil, wheel constants
        Input:
            coh = soil "cohesion" constant, [lb/in2]
            phi = soil internal friction angle, [-]
            k1 = first (cohesive) pressure-sinkage constant, [psi]
            k2 = second (frictional) pressure-sinkage constant, [lb-in]
            n = exponent [-]
            K = shear exponent constant [in], "Shear deformation modulus"
            wheel_r = wheel radius [in]
            wheel_b = wheel width [in]
            weight_W = wheel vertical weight [lb]
            skid = skid ratio, NOTE: should be larger than 0.05
        Append:
            _coh
            _phi
            _k1

            _k2
            _n
            _K
            _radius
            _width
            _weight
            _skid
            _plots
        '''
        self._coh = coh
        self._phi = phi
        self._k1 = k1
        self._k2 = k2
        self._n = n
        self._K = K
        self._radius = wheel_r
        self._width = wheel_b
        self._weight = weight_W
        self._skid = skid
        self._plots = gen_plots
        # for a driven wheel
        self._c1 = c1
        self._c2 = c2
        self._units = units # unit system, assumed to be inch/pound/second.
        # if specified else, assume we're using SI, e.g. meter/kg/sec
        
        # do any other constants need to be pre-calculated?
        self._slip_arr = py.arange(0.1,0.85,0.05)
        # point of maximum stress, th0, can be found immediately for a towed wheel
        self._th0 = self.__eval_th0_towed(self._plots)
        # for towed wheels, generally the exit angle is zero
        self._th2 = 0.0
        
    def __eval_th0_towed(self,generate_plot=False):
        """ 
        Appends:
            th0_arr: for plotting ranges of skid/slip
        """
        # based on the internal friction angle, find the angle th0 where the max
        #   stress occurs, and where shear changes directions
        def contactAngleFunc(th0, i, phi):
            out = py.tan(degToRad(45.0)-phi/2.0) - (py.cos(th0) - (1.0/(1.0+i)) ) / py.sin(th0)
            return out
        # initial guess for theta 0
        th0_initial = degToRad(25.0)
        i0 = self._skid
        solve_output = sci_opt.fsolve(contactAngleFunc,th0_initial,args=(i0,self._phi),xtol=1E-7  )
        lg.info('skid rate i = ' + str(i0) + ' th0 = ' + str(radToDeg(solve_output)) + 'degrees')
        # plot a range of values for theta_0 ( skid ), and also some values of phi
        if(generate_plot):
            i_range = self._slip_arr
            phi_arr = py.array([degToRad(31.1),degToRad(33.3),degToRad(24.0),degToRad(self._phi)] )
            th0_arr = py.zeros((len(phi_arr),len(i_range)) )  # keep the output angles here
            for row in range(0,len(phi_arr)):
                phi_curr = phi_arr[row]
                for col in range(0,len(i_range)):
                    i_curr = i_range[col]
                    th0_out = sci_opt.fsolve(contactAngleFunc,solve_output,args=(i_curr,phi_curr),xtol=1E-6 )
                    th0_arr[row,col] = th0_out               
            fig=plt.figure()
            ax = fig.add_subplot(111,title=r'$\phi$ = ' + str(radToDeg(self._phi)) + ' degrees')
            ax.plot(i_range,radToDeg(th0_arr[0,:]),i_range,radToDeg(th0_arr[1,:]),i_range,radToDeg(th0_arr[2,:]),linewidth=1.5)
            ax.set_xlabel('skid ratio')
            ax.set_ylabel(r'$\theta_0$ [degrees]' )
            ax.set_xlim([0,i_range[len(i_range)-1]+.2])
            ax.legend((str(radToDeg(phi_arr[0])),str(radToDeg(phi_arr[1])),str(radToDeg(phi_arr[2])) ),loc=2 )
            ax.grid(True)
            
            # finally, hold onto these arrays for skid, th0
            # note: only keep the first th0_arr row
            self._th0_arr = th0_arr[len(phi_arr)-1,:]
            return solve_output
                
        return solve_output
    
    # from Wong/Reece's second 1967 paper, a towed wheel
    def eval_W_integral_towed(self,slip,figNum=3):
        # individual sigma and tau terms in the weight integral
        def w_func_t1(th,th1,r,b,n,k1,k2):
            outval = sig_1(th,th1,r,b,n,k1,k2)*py.cos(th)
            return outval
        def w_func_t2(th,th0,th1,th2,r,b,n,k1,k2):
            outval = sig_2(th,th0,th1,th2,r,b,n,k1,k2)*py.cos(th)
            return outval
        def w_func_t3(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip)*py.sin(th)
            return outval
        def w_func_t4(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip)*py.sin(th)
            return outval
        # weight function, used for fsolve
        def W_towed_func(th1,th0,th2,W,r,b,n,k1,k2,phi,slip,K,c):                
            term1 = sci_int.quad(w_func_t1,th0,th1,args=(th1,r,b,n,k1,k2) )
            term2 = sci_int.quad(w_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2))
            term3 = sci_int.quad(w_func_t3,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term4 = sci_int.quad(w_func_t4,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))            
            error = r*b*(term1[0] + term2[0] +term3[0] - term4[0]) - W
            return error            
        # weight function, returns all the individual terms    
        def W_towed_func_terms(th1,th0,th2,W,r,b,n,k1,k2,phi,slip,K,c):
            term1 = sci_int.quad(w_func_t1,th0,th1,args=(th1,r,b,n,k1,k2) )
            term2 = sci_int.quad(w_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2))
            term3 = sci_int.quad(w_func_t3,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term4 = sci_int.quad(w_func_t4,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))           
            error = r*b*(term1[0] + term2[0] +term3[0] - term4[0]) - W
            return [error, r*b*term1[0], r*b*term2[0], r*b*term3[0], r*b*term4[0] ]
         
        # end helper functions
        # solve for contact angle, th_1, using equlibrium of vertical forces
        th0 = self._th0
        th2 = self._th2
        k1 = self._k1
        k2 = self._k2
        b = self._width
        r = self._radius
        n = self._n
        c = self._coh
        K = self._K
        phi = self._phi
        slip = self._skid
        W = self._weight
        # initial guess for theta 1
        th1_initial = degToRad(45.0)
        if( th1_initial/th0 < 1.5):
            th1_initial = th0*2.0
        # iterate until th1 is found
        th1 = sci_opt.fsolve(W_towed_func,th1_initial,args=(th0,th2,W,r,b,n,k1,k2,phi,slip,K,c) )
        w_error = W_towed_func(th1,th0,th2,W,r,b,n,k1,k2,phi,slip,K,c)
        lg.info('theta_1, degrees = ' + str(radToDeg(th1))+', w_error = ' + str(w_error)+'\n')
        if( self._plots):
            # if I want plots, need to find th1 for each value of skid in skid_arr
            i_range = self._slip_arr
            th0_array = self._th0_arr   # th0 = th0(skid)
            th1_out = th1
            th1_array = py.zeros(len(i_range))
            t1_arr = py.zeros(len(i_range))
            t2_arr = py.zeros(len(i_range))
            t3_arr = py.zeros(len(i_range))
            t4_arr = py.zeros(len(i_range))
            for idx in range(0,len(i_range)):
                slip_curr = i_range[idx]
                th0_curr = th0_array[idx]
                th1_out = sci_opt.fsolve(W_towed_func,th1_out,args=(th0_curr,th2,W,r,b,n,k1,k2,phi,slip_curr,K,c))
                th1_array[idx] = th1_out
                # lg.info('slip= '+str(slip_curr) + ', th1= '+str(th1_out) )
                # I want to see how each term of the weight function changes
                [error,t1,t2,t3,t4] = W_towed_func_terms(th1_out,th0_curr,th2,W,r,b,n,k1,k2,phi,slip_curr,K,c)
                # lg.info('slip = '+str(slip_curr) + ', error = ' + str(error))
                t1_arr[idx] = t1
                t2_arr[idx] = t2
                t3_arr[idx] = t3
                t4_arr[idx] = t4               
            fig=plt.figure()
            ax = fig.add_subplot(211,title='Fig. '+str(figNum) + '(d)')
            ax.plot(i_range,radToDeg(th1_array),i_range,radToDeg(self._th0_arr),linewidth=1.5)
            ax.set_xlabel('skid ratio')
            ax.set_ylabel(r'$\theta$ [degrees]' )
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend((r'$\theta_1$',r'$\theta_0$'))
            ax.grid(True)
            
            ax = fig.add_subplot(212)
            ax.plot(i_range,t1_arr,i_range,t2_arr,i_range,t3_arr,i_range,t4_arr,linewidth=1.5)
            ax.set_xlabel('skid ratio')
            if( self._units == 'ips'):
                ax.set_ylabel('weight [lb]')
            else:
                ax.set_ylabel('weight [kg]')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend((r'$\sigma_1$',r'$\sigma_2$',r'$\tau_1$',r'$\tau2$'))
            ax.grid(True)
            # we will need all the values of th1 for further plots
            self._th1_arr = th1_array
            
            # see what we get for weight
#            i_check = py.arange(20.,50.,1.)
#            for i in range(0,len(i_check)):
#                [W_check,ct1,ct2,ct3,ct4] = W_towed_func_terms(degToRad(i_check[i]),degToRad(19.0),th2,W,r,b,n,k1,k2,phi,0.3,K,c)
            # lg.info('weight check, slip = 30%, w_error='+str(W_check))
            # lg.info('term 1='+str(ct1) +', term2= '+str(ct2))
            # lg.info('term3= '+str(ct3) + ', term4= '+str(ct4))
                         
        return th1

    # from Wong/Reece's first 1967 paper, a driven wheel
    # returns theta_1, Appends theta_m
    def eval_W_integral_driven(self,slip):
        def w_func_d1(th,th1,r,b,n,k1,k2):
            outval = sig_1(th,th1,r,b,n,k1,k2)*py.cos(th)
            return outval
        def w_func_d2(th,th_m,th1,th2,r,b,n,k1,k2):
            outval = sig_2(th,th_m,th1,th2,r,b,n,k1,k2)*py.cos(th)
            return outval
        def w_func_d3(th,th1,r,b,n,k1,k2,coh,phi,K,slip):
            outval = tau_d1(th,th1,r,b,n,k1,k2,coh,phi,K,slip)*py.sin(th)
            return outval
        def w_func_d4(th,th_m,th1,th2,r,b,n,k1,k2,coh,K,phi,slip):
            outval = tau_d2(th,th_m,th1,th2,r,b,n,k1,k2,coh,K,phi,slip)*py.sin(th)
            return outval
        # returns error of weight equation, used with fsolve when finding th1
        def W_driven_func(th1,th2,W,r,b,n,k1,k2,phi,slip,K,coh,c1,c2):            
            # Eq 10/11
            th_m = theta_m(th1,c1,c2,slip)
            term1 = sci_int.quad(w_func_d1,th_m,th1,args=(th1,r,b,n,k1,k2) )
            term2 = sci_int.quad(w_func_d2,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2) )
            term3 = sci_int.quad(w_func_d3,th_m,th1,args=(th1,r,b,n,k1,k2,coh,phi,K,slip) )
            term4 = sci_int.quad(w_func_d4,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2,coh,K,phi,slip) )            
            error = r*b*(term1[0] + term2[0] + term3[0] + term4[0]) - W
            return error
        # returns the terms in the weight equation    
        def W_driven_func_terms(th1,th2,W,r,b,n,k1,k2,phi,slip,K,coh,c1,c2):
            # Eq 10/11
            th_m = theta_m(th1,c1,c2,slip)
            term1 = sci_int.quad(w_func_d1,th_m,th1,args=(th1,r,b,n,k1,k2) )
            term2 = sci_int.quad(w_func_d2,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2,coh,K,phi,slip) )
            term3 = sci_int.quad(w_func_d3,th_m,th1,args=(th1,r,b,n,k1,k2,coh,phi,K,slip) )
            term4 = sci_int.quad(w_func_d4,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2,coh,K,phi,slip) )
            
            error = r*b*(term1[0] + term2[0] + term3[0] + term4[0]) - W
            return [error, r*b*term1[0], r*b*term2[0], r*b*term3[0], r*b*term4[0] ]
         
        # end helper functions
        # solve for contact angle, th_1, using equlibrium of vertical forces
        th2 = self._th2
        k1 = self._k1
        k2 = self._k2
        b = self._width
        r = self._radius
        n = self._n
        coh = self._coh
        K = self._K
        phi = self._phi
        W = self._weight
        c1 = self._c1
        c2 = self._c2
        # initial guess for theta 1
        th1_initial = degToRad(35.0)
        if(th1_initial/theta_m(th1_initial,c1,c2,slip) < 1.5):
            th1_initial = 2.0 * theta_m(th1_initial,c1,c2,slip)
        # iterate until th1 is found
        th1 = sci_opt.fsolve(W_driven_func,th1_initial,args=(th2,W,r,b,n,k1,k2,phi,slip,K,coh,c1,c2) )
        thm = theta_m(th1,c1,c2,slip)
        self._thm = thm
        lg.info('slip rate = ' + str(slip))
        lg.info('theta_1, driven [deg] = ' + str(radToDeg(th1)) )
        lg.info('theta_m, driven [deg] = ' + str(radToDeg(thm)) )
        if( self._plots):
            # if I want plots, need to find th1 for each value of skid in skid_arr
            i_range = self._slip_arr
            th_m_array = py.zeros(len(i_range))
            th1_array = py.zeros(len(i_range))
            th_ratio = py.zeros(len(i_range))
            # t1_arr = py.zeros(len(i_range))
            # t2_arr = py.zeros(len(i_range))
            # t3_arr = py.zeros(len(i_range))
            # t4_arr = py.zeros(len(i_range))
            for idx in range(0,len(i_range)):
                slip_curr = i_range[idx]
                th1_out = sci_opt.fsolve(W_driven_func,th1,args=(th2,W,r,b,n,k1,k2,phi,slip_curr,K,coh,c1,c2),xtol=1E-5)
                th1_array[idx] = th1_out
                th_m_curr = theta_m(th1_out,c1,c2,slip_curr)
                th_m_array[idx] = th_m_curr
                th_ratio[idx] = th_m_curr / th1_out
                # I want to see how each term of the weight function changes
                # [error,t1,t2,t3,t4] = W_towed_func_terms(th1_out,th0_curr,th2,W,r,b,n,k1,k2,phi,slip_curr,K,c)
                # t1_arr[idx] = t1
                # t2_arr[idx] = t2
                # t3_arr[idx] = t3
                # t4_arr[idx] = t4
                
            fig=plt.figure()
            ax = fig.add_subplot(211,title='Fig 3(a), weight = '+str(self._weight))
            ax.plot(i_range,radToDeg(th1_array),i_range,radToDeg(th_m_array),linewidth=1.5)
            ax.set_xlabel('slip ratio')
            ax.set_ylabel(r'$\theta$ [degrees]' )
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend(('theta1','theta0'))
            ax.grid(True)
            
            ax = fig.add_subplot(212)
            ax.plot(i_range,th_ratio,linewidth=1.5)
            ax.set_xlabel('slip ratio')
            ax.set_ylabel('theta_m / theta_1')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.set_ylim([0,0.9])
            ax.grid(True)
            # we will need all the values of th1 for further plots
            self._th1_arr = th1_array
            self._thm_arr = th_m_array             
        
        W_check = W_driven_func(th1,0.0,W,r,b,n,k1,k2,phi,slip,K,coh,c1,c2)
        lg.info('weight check = ' + str(W_check) + ' for a th_1 of ' +str(th1))         
        return th1
    
    # this should really be equal to zero for the towed case
    def eval_T_integral_towed(self,th1,slip,figNum=8):
        # torque is only affected by tau terms of the T = integral
        def torque_func_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip)
            return outval
        def torque_func_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip)
            return outval
        # evaluate torque
        def T_towed_func(th1,th0,th2,r,b,n,k1,k2,c,K,phi,slip):                
            term1 = sci_int.quad(torque_func_t1,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term2 = sci_int.quad(torque_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))
            outval = r*r*b*(term1[0] - term2[0])    # same as the error, when T=0
            return outval
        # eval. torque, return the integral terms also    
        def T_towed_func_terms(th1,th0,th2,r,b,n,k1,k2,c,K,phi,slip):
            term1 = sci_int.quad(torque_func_t1,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term2 = sci_int.quad(torque_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))
            outval = r*r*b*(term1[0] - term2[0])
            return [outval, r*r*b*term1[0], r*r*b*term2[0] ]
            
        th0 = self._th0
        th2 = self._th2
        r = self._radius
        b = self._width
        c = self._coh
        k1 = self._k1
        k2 = self._k2
        n = self._n
        phi = self._phi
        K = self._K
        slip = self._skid        
        # should really be solving for th1 w.r.t. W, T and F
        [calc_torque, t1, t2] = T_towed_func_terms(th1,th0,th2,r,b,n,k1,k2,c,K,phi,slip)

        lg.info('towed Torque = ' + str(calc_torque))
        lg.info('torque front = ' + str(t1) + ' , bottom = ' + str(t2) )
        if( self._plots):
            i_range = self._slip_arr
            th0_array = self._th0_arr
            th1_array = self._th1_arr
            torque_arr = py.zeros(len(i_range))
            t1_arr = py.zeros(len(i_range))
            t2_arr = py.zeros(len(i_range))
            for idx in range(0,len(i_range)):
                slip_curr = i_range[idx]
                th0_curr = th0_array[idx]
                th1_curr = th1_array[idx]
                [torque_curr,t1,t2] = T_towed_func_terms(th1_curr,th0_curr,th2,r,b,n,k1,k2,c,K,phi,slip_curr)
                torque_arr[idx] = torque_curr
                t1_arr[idx] = t1
                t2_arr[idx] = t2
                
            fig=plt.figure()
            ax = fig.add_subplot(111,title='Fig. ' + str(figNum) +'(a)')
            ax.plot(i_range,torque_arr,'k--',i_range,t1_arr,i_range,t2_arr,linewidth=1.5)
            ax.set_xlabel('skid ratio')
            if( self._units == 'ips'):
                ax.set_ylabel('Torque [lb-in]' )
            else:
                ax.set_ylabel('Torque [N-m]')
            # ax.set_ylabel('Torque [N-m]')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend(('T',r'$\tau_1$',r'$\tau_2$')) #,loc=2)
            ax.grid(True)
            
        return calc_torque

    # driven wheel should not be zero  
    def eval_T_integral_driven(self,th1,slip,figNum=8):
        # integral terms
        def torque_func_d1(th,th1,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_d1(th,th1,r,b,n,k1,k2,c,phi,K,slip)
            return outval
        def torque_func_d2(th,th_m,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_d2(th,th_m,th1,th2,r,b,n,k1,k2,c,K,phi,slip)
            return outval
        # driven torque function
        def T_driven_func(th1,th_m,th2,r,b,n,k1,k2,c,K,phi,slip):
            term1 = sci_int.quad(torque_func_d1,th_m,th1,args=(th1,r,b,n,k1,k2,c,K,phi,slip))
            term2 = sci_int.quad(torque_func_d2,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2,c,K,phi,slip))
            outval = r*r*b*(term1[0] + term2[0])
            return outval
        # driven torque function, with each term returned
        def T_driven_func_terms(th1,th_m,th2,r,b,n,k1,k2,c,K,phi,slip):
            term1 = sci_int.quad(torque_func_d1,th_m,th1,args=(th1,r,b,n,k1,k2,c,K,phi,slip))
            term2 = sci_int.quad(torque_func_d2,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2,c,K,phi,slip))
            outval = r*r*b*(term1[0] + term2[0])
            return [outval, r*r*b*term1[0], r*r*b*term2[0] ]
            

        # eval_T_integral_towed HELPER FUNCTIONS END HERE
        th_m = self._thm
        th2 = self._th2
        r = self._radius
        b = self._width
        c = self._coh
        k1 = self._k1
        k2 = self._k2
        n = self._n
        phi = self._phi
        K = self._K      
        # should really be solving for th1 w.r.t. W, T and F
        [calc_torque, t1, t2] = T_driven_func_terms(th1,th_m,th2,r,b,n,k1,k2,c,K,phi,slip)

        lg.info('driven Torque = ' + str(calc_torque))
        lg.info('driven terms, front = ' + str(t1) + ' , bottom = ' + str(t2) )
        if( self._plots):
            i_range = self._slip_arr
            thm_array = self._thm_arr
            th1_array = self._th1_arr
            torque_arr = py.zeros(len(i_range))
            t1_arr = py.zeros(len(i_range))
            t2_arr = py.zeros(len(i_range))
            for idx in range(0,len(i_range)):
                slip_curr = i_range[idx]
                thm_curr = thm_array[idx]
                th1_curr = th1_array[idx]
                [torque_curr,t1,t2] = T_driven_func_terms(th1_curr,thm_curr,th2,r,b,n,k1,k2,c,K,phi,slip_curr)
                torque_arr[idx] = torque_curr
                t1_arr[idx] = t1
                t2_arr[idx] = t2
                
            fig=plt.figure()
            ax = fig.add_subplot(111,title='Fig.'+str(figNum)+'(c)')
            ax.plot(i_range,torque_arr,'k--',i_range,t1_arr,i_range,t2_arr,linewidth=1.5)
            ax.set_xlabel('slip ratio')
            if( self._units == 'ips'):
                ax.set_ylabel('Torque [lb-in]' )
            else:
                ax.set_ylabel('Torque [N-m]')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend(('T',r'$\tau_1$',r'$\tau_2$'))
            ax.grid(True)
            
        return calc_torque
    
    # required towing force    
    def eval_F_integral_towed(self,th1,slip,figNum=8):
        # terms in the F = integral
        def F_func_t1(th,th1,r,b,n,k1,k2):
            outval = sig_1(th,th1,r,b,n,k1,k2)*py.sin(th)
            return outval
        def F_func_t2(th,th0,th1,th2,r,b,n,k1,k2):
            outval = sig_2(th,th0,th1,th2,r,b,n,k1,k2)*py.sin(th)
            return outval
        def F_func_t3(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip)*py.cos(th)
            return outval
        def F_func_t4(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip)*py.cos(th)
            return outval
        # total longitudinal force acting on the wheel
        def F_towed_func(th1,th0,th2,r,b,n,k1,k2,c,K,phi,slip):                
            term1 = sci_int.quad(F_func_t1,th0,th1,args=(th1,r,b,n,k1,k2))
            term2 = sci_int.quad(F_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2))
            term3 = sci_int.quad(F_func_t3,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term4 = sci_int.quad(F_func_t4,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))
            outval = r*b*(term1[0]+term2[0]-term3[0]+term4[0])
            return outval
        
        # longitudinal force, and individual terms in the integral
        def F_towed_func_terms(th1,th0,th2,r,b,n,k1,k2,c,K,phi,slip):
            term1 = sci_int.quad(F_func_t1,th0,th1,args=(th1,r,b,n,k1,k2))
            term2 = sci_int.quad(F_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2))
            term3 = sci_int.quad(F_func_t3,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term4 = sci_int.quad(F_func_t4,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))
            outval = r*b*(term1[0]+term2[0]-term3[0]+term4[0])
            return [outval, r*b*term1[0], r*b*term2[0], r*b*term3[0], r*b*term4[0] ]            

        # eval_F_integral_towed() HELPER FUNCTIONS END HERE
        th0 = self._th0
        th2 = self._th2
        r = self._radius
        b = self._width
        k1 = self._k1
        k2 = self._k2
        n = self._n
        phi = self._phi
        K = self._K
        c = self._coh            
        # don't know what this should be apriori
        calc_force = F_towed_func(th1,th0,th2,r,b,n,k1,k2,c,K,phi,slip)

        lg.info('towing Force = ' + str(calc_force))
        if( self._plots):
            i_range = self._slip_arr
            th0_array = self._th0_arr
            th1_array = self._th1_arr
            force_arr = py.zeros(len(i_range))
            t1_arr = py.zeros(len(i_range))
            t2_arr = py.zeros(len(i_range))
            t3_arr = py.zeros(len(i_range))
            t4_arr = py.zeros(len(i_range))
            for idx in range(0,len(i_range)):
                slip_curr = i_range[idx]
                th0_curr = th0_array[idx]
                th1_curr = th1_array[idx]
                [force_curr,t1,t2,t3,t4] = F_towed_func_terms(th1_curr,th0_curr,th2,r,b,n,k1,k2,c,K,phi,slip_curr)
                force_arr[idx] = force_curr
                t1_arr[idx] = t1
                t2_arr[idx] = t2
                t3_arr[idx] = t3
                t4_arr[idx] = t4      
            fig=plt.figure()
            ax = fig.add_subplot(211,title='Fig.' + str(figNum) + '(b) [top] and Fig.'+str(figNum) + '(c) [bottom]')
            ax.plot(i_range,t3_arr-t4_arr,i_range,t1_arr+t2_arr,linewidth=1.5)
            ax.set_xlabel('skid ratio')
            if( self._units == 'ips'):            
                ax.set_ylabel('Force [lb]' )
            else:
                ax.set_ylabel('Force [N]')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend(('H','R'))
            ax.grid(True)
            
            ax = fig.add_subplot(212)
            ax.plot(i_range,force_arr,'k--',i_range,t1_arr,i_range,t2_arr,i_range,t3_arr,i_range,t4_arr,linewidth=1.5)
            # ax.plot(i_range,force_arr,'k--',linewidth=1.5)            
            ax.set_xlabel('skid ratio')
            if( self._units == 'ips'):         
                ax.set_ylabel('Force [lb]')
            else:
                ax.set_ylabel('Force [N]')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend(('F',r'$\sigma_1$',r'$\sigma_2$',r'$\tau_1$',r'$\tau_2$'))
            ax.grid(True)
            
        return calc_force
        
    # driving force    
    def eval_F_integral_driven(self,th1,slip,figNum=8):
        # integral terms
        def F_func_d1(th,th1,r,b,n,k1,k2):
            outval = py.sin(th)*sig_1(th,th1,r,b,n,k1,k2)
            return outval
        def F_func_d2(th,th_m,th1,th2,r,b,n,k1,k2):
            outval = py.sin(th)*sig_2(th,th_m,th1,th2,r,b,n,k1,k2)
            return outval
        def F_func_d3(th,th1,r,b,n,k1,k2,c,K,phi,slip):
            outval = py.cos(th)*tau_d1(th,th1,r,b,n,k1,k2,c,phi,K,slip)
            return outval
        def F_func_d4(th,th_m,th1,th2,r,b,n,c,k1,k2,K,phi,slip):
            outval = py.cos(th)*tau_d2(th,th_m,th1,th2,r,b,n,k1,k2,c,K,phi,slip)
            return outval
        # driven wheel force
        def F_driven_func(th1,th_m,th2,r,b,n,k1,k2,c,K,phi,slip):            
            term1 = sci_int.quad(F_func_d1,th_m,th1,args=(th1,r,b,n,k1,k2))
            term2 = sci_int.quad(F_func_d2,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2))
            term3 = sci_int.quad(F_func_d3,th_m,th1,args=(th1,r,b,n,k1,k2,c,K,phi,slip))
            term4 = sci_int.quad(F_func_d4,th2,th_m,args=(th_m,th1,th2,r,b,n,c,k1,k2,K,phi,slip))
            outval = r*b*(-term1[0]-term2[0]+term3[0]+term4[0])
            return outval
        # driven wheel force, with terms
        def F_driven_func_terms(th1,th_m,th2,r,b,n,k1,k2,c,K,phi,slip):
            term1 = sci_int.quad(F_func_d1,th_m,th1,args=(th1,r,b,n,k1,k2))
            term2 = sci_int.quad(F_func_d2,th2,th_m,args=(th_m,th1,th2,r,b,n,k1,k2))
            term3 = sci_int.quad(F_func_d3,th_m,th1,args=(th1,r,b,n,k1,k2,c,K,phi,slip))
            term4 = sci_int.quad(F_func_d4,th2,th_m,args=(th_m,th1,th2,r,b,n,c,k1,k2,K,phi,slip))
            outval = r*b*(-term1[0]-term2[0]+term3[0]+term4[0])
            return [outval, r*b*term1[0], r*b*term2[0], r*b*term3[0], r*b*term4[0] ]            

        # eval_F_integral_towed() HELPER FUNCTIONS END HERE
        th_m = self._thm
        th2 = self._th2
        r = self._radius
        b = self._width
        k1 = self._k1
        k2 = self._k2
        n = self._n
        phi = self._phi
        K = self._K
        c = self._coh            
        # don't know what this should be apriori
        calc_force = F_driven_func(th1,th_m,th2,r,b,n,k1,k2,c,K,phi,slip)

        lg.info('drawbar pull = ' + str(calc_force))
        if( self._plots):
            i_range = self._slip_arr
            thm_array = self._thm_arr
            th1_array = self._th1_arr
            force_arr = py.zeros(len(i_range))
            t1_arr = py.zeros(len(i_range))
            t2_arr = py.zeros(len(i_range))
            t3_arr = py.zeros(len(i_range))
            t4_arr = py.zeros(len(i_range))
            for idx in range(0,len(i_range)):
                slip_curr = i_range[idx]
                thm_curr = thm_array[idx]
                th1_curr = th1_array[idx]
                [force_curr,t1,t2,t3,t4] = F_driven_func_terms(th1_curr,thm_curr,th2,r,b,n,k1,k2,c,K,phi,slip_curr)
                force_arr[idx] = force_curr
                t1_arr[idx] = t1
                t2_arr[idx] = t2
                t3_arr[idx] = t3
                t4_arr[idx] = t4         
            fig=plt.figure()
            ax = fig.add_subplot(211,title='Fig.' + str(figNum) + '(a) [top] and Fig.'+str(figNum) + '(b) [bottom]')
            ax.plot(i_range,t3_arr+t4_arr,i_range,t1_arr+t2_arr,linewidth=1.5)
            # ax.set_xlabel('slip ratio')
            ax.set_ylabel('Force [lb]')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend(('H','R'))
            ax.grid(True)
            
            ax = fig.add_subplot(212,title='H and R')
            ax.plot(i_range,force_arr,'k--',i_range,t1_arr,i_range,t2_arr,i_range,t3_arr,i_range,t4_arr,linewidth=1.5)
            ax.set_xlabel('slip ratio')
            if( self._units == 'ips'):
                ax.set_ylabel('Force [lb]' )
            else:
                ax.set_ylabel('Force [N]')
            ax.set_xlim([0,i_range[len(i_range)-1]+0.2])
            ax.legend(('D',r'$\sigma_1$',r'$\sigma_2$',r'$\tau_1$',r'$\tau_2$'))
            ax.grid(True)    
            
        return calc_force        
     
     # plot the normal, shear stress distributions along theta, towed wheel 
    def plot_sigTau_towed(self,th1_cs,skid,figNum=11):
        th0 = self._th0
        th2 = self._th2
        r = self._radius
        b = self._width
        k1 = self._k1
        k2 = self._k2
        n = self._n
        phi = self._phi
        K = self._K
        c = self._coh
        
        incr = (th1_cs - th2) / 100.0    # plot increment
        th_arr = py.arange(0,th1_cs + incr, incr) # find sigma, tau at these discrete vals
        sig_arr = py.zeros(len(th_arr))
        tau_arr = py.zeros(len(th_arr))
        slip_arr = py.zeros(len(th_arr))
        
        for idx in range(0,len(th_arr)):
            th = th_arr[idx]
            if(th < th0):
                # we're in the bototm region
                sig_curr = sig_2(th,th0,th1_cs,th2,r,b,n,k1,k2)
                tau_curr = -tau_t2(th,th0,th1_cs,th2,r,b,n,k1,k2,c,K,phi,skid)
                slip_j = j2(th,th0,r,skid)
                sig_arr[idx] = sig_curr
                tau_arr[idx] = tau_curr
                slip_arr[idx] = slip_j
                
            else:
                # we're in the top region ()
                sig_curr = sig_1(th, th1_cs,r,b,n,k1,k2)
                tau_curr = tau_t1(th,th0,th1_cs,r,b,n,k1,k2,c,K,phi,skid)
                slip_j = j1(th,th0,th1_cs,r,skid)
                sig_arr[idx] = sig_curr
                tau_arr[idx] = tau_curr
                slip_arr[idx] = slip_j
                
        if( self._plots):        
            fig = plt.figure()
            ax = fig.add_subplot(211,title='Fig. ' + str(figNum) +' skid=' + str(skid))
            ax.plot(radToDeg(th_arr),sig_arr,radToDeg(th_arr),tau_arr,linewidth=1.5)
            ax.set_xlabel('theta [deg]')
            if( self._units == 'ips'):
                ax.set_ylabel('stress [psi]')
            else:
                ax.set_ylabel('stress [Pa]')
            ax.legend((r'$\sigma$($\theta$)',r'$\tau$($\theta$)'))
            ax.grid(True)
            # take a look at what I"m using for slip displacement also
            ax = fig.add_subplot(212)
            ax.plot(radToDeg(th_arr),slip_arr,linewidth=1.5)
            ax.set_xlabel('theta [deg]')
            if( self._units == 'ips'):
                ax.set_ylabel('slip disp.[in]')
            else:
                ax.set_ylabel('slip disp.[m]')
            ax.grid(True)
            
             # polar plots
            fig=plt.figure()
            ax=fig.add_subplot(111,projection='polar')
            ax.plot(th_arr,sig_arr/1000.,'b',linewidth=1.5)
            ax.plot(th_arr,tau_arr/1000.,'r--',linewidth=1.5)
            # fix the axes
            ax.grid(True)
            if( self._units == 'ips'):
                leg = ax.legend((r'$\sigma$ [kip]',r'$\tau$'))
            else:
                leg = ax.legend((r'$\sigma$ [kPa]',r'$\tau$'))
            leg.draggable()
            ax.set_theta_zero_location('S')
            # also, draw the tire
            polar_r_offset = py.average(sig_arr)/1000.
            theta = py.arange(0.,2.*py.pi+0.05,0.05)
            tire_pos = py.zeros(len(theta))
            ax.plot(theta,tire_pos,color='k',linewidth=1.0)
            ax.set_rmin(-polar_r_offset)
            ax.set_title(r'towed wheel stresses,  $\theta_1$ = %4.3f [rad]' %th1_cs)
            ax.set_thetagrids([-10,0,10,20,30,40,50,60])
            
        return [sig_arr, tau_arr]
    
    # plot the normal, shear stress distributions along theta, driven wheel
    # return the y-vals for [sigma, tau], so I can plot lots of these    
    def plot_sigTau_driven(self,th1_cs,slip,figNum=11,plotBekker=False):
        th_m = self._thm
        th2 = self._th2
        r = self._radius
        b = self._width
        k1 = self._k1
        k2 = self._k2
        n = self._n
        phi = self._phi
        K = self._K
        c = self._coh
        
        incr = (th1_cs - th2) / 100.0    # plot increment
        th_arr = py.arange(0,th1_cs, incr) # find sigma, tau at these discrete vals
        sig_arr = py.zeros(len(th_arr))
        tau_arr = py.zeros(len(th_arr))
        slip_arr = py.zeros(len(th_arr))
        
        for idx in range(0,len(th_arr)):
            th = th_arr[idx]
            if(th <= th_m):
                # we're in the bototm region
                sig_curr = sig_2(th,th_m,th1_cs,th2,r,b,n,k1,k2)
                tau_curr = tau_d2(th,th_m,th1_cs,th2,r,b,n,k1,k2,c,K,phi,slip)
                slip_j = jdriven(th,th1_cs,r,slip)
                sig_arr[idx] = sig_curr
                tau_arr[idx] = tau_curr
                slip_arr[idx] = slip_j
                
            else:
                # we're in the top region ()
                sig_curr = sig_1(th, th1_cs,r,b,n,k1,k2)
                tau_curr = tau_d1(th,th1_cs,r,b,n,k1,k2,c,phi,K,slip)
                slip_j = jdriven(th,th1_cs,r,slip)
                sig_arr[idx] = sig_curr
                tau_arr[idx] = tau_curr
                slip_arr[idx] = slip_j
                
        
        if( self._plots):        
            fig = plt.figure()
            ax = fig.add_subplot(211,title='Fig.'+str(figNum) )
            ax.plot(radToDeg(th_arr),sig_arr,radToDeg(th_arr),tau_arr,linewidth=1.5)
            ax.set_xlabel('theta [deg]')
            ax.set_ylabel('stress [psi]')
            ax.legend((r'$\sigma$',r'$\tau$'))
            ax.grid(True)
            # can also plot Bekker's solution
            if(plotBekker):
                th_bek, sig_bek, tau_bek = self.get_sigTau_Bekker_driven(slip)
                ax.plot(radToDeg(th_bek),sig_bek, radToDeg(th_bek), tau_bek, linewidth=1.5)
                ax.legend((r'$\sigma$',r'$\tau$',r'$\sigma_bek$',r'$\tau_bek$'))
                
            # take a look at what I"m using for slip displacement also
            ax = fig.add_subplot(212)
            ax.plot(radToDeg(th_arr),slip_arr,linewidth=1.5)
            ax.set_xlabel('theta [deg]')
            if( self._units == 'ips'):
                ax.set_ylabel('slip disp.[in]')
            else:
                ax.set_ylabel('slip disp.[m]')
            ax.grid(True)
            
            # polar plots
            fig=plt.figure()
            ax=fig.add_subplot(111,projection='polar')
            ax.plot(th_arr,sig_arr,'b',linewidth=1.5)
            ax.plot(th_arr,tau_arr,'r--',linewidth=1.5)
            # fix the axes
            ax.grid(True)
            if( self._units == 'ips'):
                leg = ax.legend((r'$\sigma$ [psi]',r'$\tau$'))
            else:
                leg = ax.legend((r'$\sigma$ [Pa]',r'$\tau$'))
            leg.draggable()
            ax.set_theta_zero_location('S')
            # also, draw the tire
            polar_r_offset = py.average(sig_arr)
            theta = py.arange(0.,2.*py.pi+0.05,0.05)
            tire_pos = py.zeros(len(theta))
            ax.plot(theta,tire_pos,color='k',linewidth=1.0)
            ax.set_rmin(-polar_r_offset)
            ax.set_title(r'driven wheel stresses, $\theta_1$ = %4.3f [rad]' %th1_cs)
            ax.set_thetagrids([-10,0,10,20,30,40,50,60])
            
    def get_sigTau_Bekker_driven(self,slip):
        '''
        Returns:
            [theta_arr,sigma_array,tau_array]
            for plotting
        '''
        # this is nice, only solve for 1 normal and shear equation
        def w_func_bek_d1(th,th1,r,b,n,k1,k2):
            outval = sig_bek(th,th1,r,b,n,k1,k2)*py.cos(th)
            return outval
        def w_func_bek_d2(th,th1,r,b,n,k1,k2,coh,phi,K,slip):
            outval = tau_bek(th,th1,r,b,n,k1,k2,coh,phi,K,slip)*py.sin(th)
            return outval
        # returns error of weight equation, used with fsolve when finding th1
        def W_bek_driven_func(th1,th2,W,r,b,n,k1,k2,phi,slip,K,coh,c1,c2):
            term1 = sci_int.quad(w_func_bek_d1, th2, th1, args=(th1,r,b,n,k1,k2) )
            term2 = sci_int.quad(w_func_bek_d2, th2, th1, args=(th1,r,b,n,k1,k2,coh,phi,K,slip) )
            error = r*b*(term1[0] + term2[0] ) - W
            return error
         
        # end helper functions
        # solve for contact angle, th_1, using equlibrium of vertical forces
        th2 = self._th2
        k1 = self._k1
        k2 = self._k2
        b = self._width
        r = self._radius
        n = self._n
        coh = self._coh
        K = self._K
        phi = self._phi
        W = self._weight
        c1 = self._c1
        c2 = self._c2
        # initial guess for theta 1
        th1_initial = degToRad(35.0)
        # iterate until th1 is found
        th1 = sci_opt.fsolve(W_bek_driven_func,th1_initial,args=(th2,W,r,b,n,k1,k2,phi,slip,K,coh,c1,c2) )
        lg.info('slip rate = ' + str(slip))
        lg.info('theta_1, Bekker[deg] = ' + str(radToDeg(th1)) )
        incr = (th1 - th2)/100.
        theta_arr = py.arange(th2,th1+incr,incr)
        sig_arr = sig_bek(theta_arr,th1,r,b,n,k1,k2)
        tau_arr = tau_bek(theta_arr,th1,r,b,n,k1,k2,coh,phi,K,slip)
        
        W_check = W_bek_driven_func(th1,0.0,W,r,b,n,k1,k2,phi,slip,K,coh,c1,c2)
        lg.info('weight check, Bekker = ' + str(W_check) + ', ought to be zero')    
        
        # return theta, sigma, tau, for plotting
        return [theta_arr, sig_arr, tau_arr]  

        
        return [sig_arr,tau_arr]

    # plot the sinkage for a driven wheel as a function of slip
    def plot_z0_driven(self,figNum=8):
        """ Only defined for a driven wheel, since a towed wheel will create
        significant build-up in front of the wheel
        """
        fig = plt.figure()
        ax = fig.add_subplot(111,title='weight= ' + str(self._weight) + ', Fig.'+str(figNum)+'(d)' )
        z_array = (1.0-py.cos(self._th1_arr))*self._radius
        ax.plot(self._slip_arr, z_array, linewidth = 1.5 )
        ax.set_xlabel('slip ratio')
        ax.set_ylabel('sinkage [inches]')
        ax.grid(True)
        
    # use this function to directly solve for the unknowns: th0, th1 and torque(skid)=0    
    def eval_towed_vars(self, skid_guess=0.3, th0_guess=py.pi/8.0, th1_guess = py.pi/4.0):
        # individual sigma and tau terms in the weight integral
        def w_func_t1(th,th1,r,b,n,k1,k2):
            outval = sig_1(th,th1,r,b,n,k1,k2)*py.cos(th)
            return outval
        def w_func_t2(th,th0,th1,th2,r,b,n,k1,k2):
            outval = sig_2(th,th0,th1,th2,r,b,n,k1,k2)*py.cos(th)
            return outval
        def w_func_t3(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip)*py.sin(th)
            return outval
        def w_func_t4(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip)*py.sin(th)
            return outval
        # weight function, used for fsolve
        def W_towed_func(th1,th0,th2,W,r,b,n,k1,k2,phi,slip,K,c):                
            term1 = sci_int.quad(w_func_t1,th0,th1,args=(th1,r,b,n,k1,k2) )
            term2 = sci_int.quad(w_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2))
            term3 = sci_int.quad(w_func_t3,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term4 = sci_int.quad(w_func_t4,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))            
            error = r*b*(term1[0] + term2[0] +term3[0] - term4[0]) - W
            return error
        # torque is only affected by tau terms of the T = integral
        def torque_func_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t1(th,th0,th1,r,b,n,k1,k2,c,K,phi,slip)
            return outval
        def torque_func_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip):
            outval = tau_t2(th,th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip)
            return outval
        # evaluate torque
        def T_towed_func(th1,th0,th2,r,b,n,k1,k2,c,K,phi,slip):                
            term1 = sci_int.quad(torque_func_t1,th0,th1,args=(th0,th1,r,b,n,k1,k2,c,K,phi,slip))
            term2 = sci_int.quad(torque_func_t2,th2,th0,args=(th0,th1,th2,r,b,n,k1,k2,c,K,phi,slip))
            outval = r*r*b*(term1[0] - term2[0])    # same as the error, when T=0
            return outval
        # solve unkowns using equlibrium of vertical forces
        th0 = self._th0
        th2 = self._th2
        k1 = self._k1
        k2 = self._k2
        b = self._width
        r = self._radius
        n = self._n
        c = self._coh
        K = self._K
        phi = self._phi
        slip = self._skid
        W = self._weight
        
        return [2]
        
        
        
if __name__ == '__main__':
    lg.basicConfig(fileName = 'logFile.log',
                   level=lg.INFO,
                   format='%(message)s')
    # default font size
    font = {'size' : 14}
    plt.rc('font', **font)

    # Now, implement the procedure to find the steady-state sinkage, slip, and
    #   motion resistance according to Wong/Reece [1967]
    
    # use the constants from Table 1 to validate procedure
    skid = 0.216
    skid_ls = 0.445
    skid_ls_wide = .52
#    slip = 0.221    # so I can match Fig 5. for pressure(theta), paper 1
    coh = 0.1   # cohesion
    phi = degToRad(33.3)  # degrees
    k1 = 20.0
    k2 = 2.5
    n = 0.47
    K = 1.5 # shear modulu=s
    rad = 24.7  # wheel radius, inches
    wid = 6.0   # wheel width, inches
    Weight = 2006.0 # wheel weight, lbs
    # note: must input a value for skid rate, from [0 - x], but should be at least 0.1
    Compact_sand = WongReece(coh,phi,k1,k2,n,K,rad,wid,Weight,skid,0.43,0.32,True)
    Loose_sand = WongReece(0.12,degToRad(31.1),0.0,2.0,1.15,1.5,24.7,6.0,1981,skid,0.18,0.32,True)
    Loose_sand_wide = WongReece(0.12,degToRad(31.1),0.0,2.0,1.15,1.5,24.7,12.0,2080,skid,0.18,0.32,True)
    
    # neglegct rut recovery, th2 = 0
    # solve for th1 in W = f(th1)
    th1_cs = Compact_sand.eval_W_integral_towed(skid,8)
    Torque = Compact_sand.eval_T_integral_towed(th1_cs, skid,8)
    Force = Compact_sand.eval_F_integral_towed(th1_cs, skid,8)
    Compact_sand.plot_sigTau_towed(th1_cs,skid,11)
#    [th0_check,th1_check,T_check,F_check] = Compact_sand.eval_towed_vars(skid)
    '''
    #on loose sand
    th1_ls = Loose_sand.eval_W_integral_towed(skid_ls,9)
    Loose_sand.eval_T_integral_towed(th1_ls, skid_ls,9)
    Loose_sand.eval_F_integral_towed(th1_ls,skid_ls,9)
    Loose_sand.plot_sigTau_towed(th1_ls, skid_ls,12)
    
    # loose sand, width = 12"
    th1_ls = Loose_sand_wide.eval_W_integral_towed(skid_ls_wide,10)
    Loose_sand_wide.eval_T_integral_towed(th1_ls, skid_ls_wide,10)
    Loose_sand_wide.eval_F_integral_towed(th1_ls,skid_ls_wide,10)
    Loose_sand_wide.plot_sigTau_towed(th1_ls, skid_ls_wide,13)
   '''    
    py.show()