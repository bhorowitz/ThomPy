import sys
from MDSplus import *
import MDSplus
import re
import numpy as numpy
from scipy import integrate
from scipy import optimize
tree = 0

array = set(['p4','p6','t1','t2','L2times','L1times'])


def gtefitP(wv_all,nshelf,nkt,Nsc,dNsc,cg,Rguess,QT,fb,theta,kt,ishelf):
     
     i = 1
     #wave = conf.wv_all[i,:]
     pos = [0,1,3,4,5]
     a_tot = numpy.zeros([2,nkt])
     a_tot_above = numpy.zeros([2,nkt])
     print "ishelf", ishelf
     wave = wv_all[ishelf,:]
     wave = wave[wave>0]
     for i in range(0,nkt):
          ii = kt[i]
          #not robust?!
          amp = numpy.sum(Nsc[ishelf,ii,pos])/1e3
          a = [amp , numpy.polyval(cg[ishelf,::-1],Rguess[ishelf,ii])] #parameters [ne,Te]
          x = numpy.arange(5)
          y=Nsc[ishelf,ii,pos]
          dy=dNsc[ishelf,ii,pos]
          iter1=8

          Tfb = fb[ishelf]
    
          w=1./dy**2
          theta_ = theta[ishelf]
          def gte(a,QT=QT,Tfb=Tfb,theta_=theta_,wave=wave):
               Sc = ts_spectra([a[1]],[theta_],wave)
               p = numpy.zeros([5])
               for k in range(0,5):
                    l = Sc*QT[k,wave>0]*Tfb[wave>0]
                    #l = l[0,0]
                    #l = l[wave>0]
                    int = integrate.simps(l,wave)
                    p[k]=int
               f = a[0]*p
               return f

          def f(a,y=y):
               return gte(a)-y
          def g(y,a):
               return gte(a)-y
     #a_maybe = optimize.curve_fit(g,y,numpy.zeros([5]))[0]
          #a_new_below = optimize.leastsq(f,[0.01,0.01])[0]
          a_new_above = optimize.leastsq(f,a,maxfev = 6000)[0]
          #a_tot[:,ishelf,i]=a_new_above[:]
          a_tot_above[:,i]= a_new_above[:]
     return a_tot_above

def ts_spectra(Te,theta,ls,li=10640):
     #Functional form of the Thompson Scattering Spectrum
     
     #Te,electron temperature [keV]
     #Theta, scattering angle [Degrees]
     #ls, wavelength [angstroms]
     #li, laser intensity (default is 10640) 

     Te = numpy.array(Te)
     theta = numpy.array(theta)
     ls = numpy.array(ls)

     nTe = Te.size
     ntheta = theta.size
     nls = ls.size
     S = numpy.zeros((ntheta,nTe,nls))
     thr = numpy.radians(theta)
     x = ls/li - 1
     alpha = 255.5/Te
     C = numpy.sqrt(alpha/numpy.pi)*(1 - 15/(16*alpha)+345/(512*alpha*alpha))

     A1 =(1+x)*(1+x)*numpy.sqrt(2*(1-numpy.cos(thr))*(1+x)+x*x) 
     B = numpy.sqrt(1+x*x/(2*(1-numpy.cos(thr))*(1+x)))-1
     for i in range(0,nTe):
          S[i]=C[i]/A1*numpy.exp(-2*alpha[i]*B)
     #returns number of scattered photons
     return S

def gte(a,QT,Tfb,theta_,wave):
     # Corrected number of scatter photons when including gain of filter
     # and quantum effeciency of electronics.
     
     Sc = ts_spectra([a[1]],[theta_],wave)

     p = numpy.zeros([5]) 
     for k in range(0,5):
          l = Sc*QT[k,:]*Tfb
          int = integrate.simps(l,wave)
          p[k]=int
     f = a[0]*p

     dTe=a[1]/1e3
     if dTe < 0: dTe=1e-4
     Sc1=ts_spectra([a[1]+dTe],[theta_],[wave])
     pder = numpy.zeros([2,5])
     pder[0,:]=p
     for i in range(0,5):
          pder[1,i]=a[0]*(integrate.simps(Sc1*QT[i,:]*Tfb,wave)-p[i])/dTe
     return f,pder

def Rguess(Nsc, nchannels,nkl):
     #Initial guess for Thomspon fiting routine... could be improved to speed
     #up computation time.
     Rguess=numpy.zeros((nchannels,nkl))
     Nsc1=Nsc[0:20,:,3:6]
     Nsc1=(Nsc1[0:20,:,0] + Nsc1[0:20,:,1] + Nsc1[0:20,:,2])
     Rguess[0:20,:]=Nsc[0:20,:,1]/Nsc1
     Nsc1 = Nsc
     Nsc1=(Nsc1[20:30,:,0] + Nsc1[20:30,:,1] + Nsc1[20:30,:,2])
     Rguess[20:30,:]=Nsc[20:30,:,2]/Nsc1
     return Rguess

     
def dblarr(x,y):
     #Compatibility function with MatLab
     r = numpy.zeros((y,x))
     return r

def  photoelectrons(it_las01, it_las02, data, nchannels,nmux,gain,Times_all,cfb, Gfast, cF):
     #physics constants loaded
     e=1.6022e-19  #coulombs
     Rfb=5e4       # ohms
     nmux=12
     tau=50e-9

     #ksort is array of indicies of data where there is data
     ksort = numpy.append(it_las01,it_las02)
     ksort = numpy.sort(ksort)
     ksort = numpy.unique(ksort)

     #roundabout way of initializing the mpts_phase array
     a = numpy.zeros(10)
     r = numpy.append([],a+1)
     r = numpy.append(r,a+2)
     mpts_phase = numpy.append(r,a+3)

     
     dl=data[0:nchannels,ksort,:]
     dl0 = numpy.zeros((nchannels,nmux))
     nkl = ksort.size

     indices = numpy.arange(Times_all.size)
     conditions = (Times_all < -0.25)
     kx = indices[conditions]

     #initializing arrays for number of scattered photoelectrons and their error
     Nsc = numpy.zeros((nchannels,nkl,12))
     dNsc = numpy.zeros((nchannels,nkl,12))
     M = gain #silly... fix eventually

  

     for j in range(0,nchannels):
          #casework for gslow
          if mpts_phase[j] == 1 or mpts_phase[j] == 2:
               gslow=[20.,20.,40.,40.,40.,40.0]
          if mpts_phase[j] ==  3:
               gslow=[40.,40.,40.,40.,1.,1.] # last two are dummies             
          for i in range(0,6):
               dl0[j,i]  = numpy.std(data[j,kx,i])  # fast
               dl0[j,i+6]= numpy.std(data[j,kx,i+6])  # slow
               
               dl[j,0:ksort.size,i]  =data[j,ksort,i]- numpy.mean(data[j,kx,i]) #subtract off DC offset
               dl[j,0:ksort.size,i+6]=data[j,ksort,i+6]-numpy.mean(data[j,kx,i+6]) #subtract off DC offset
               
               Nsc[j,:,i]=dl[j,:,i]*cfb[j,i]/M[j,i]/e/Gfast[j]  # fast
               constants = cF[j,i,:]
               Fex=numpy.polyval(constants[::-1],M[j,i])

               #error estimation
               dNsc[j,:,i]= numpy.sqrt((dl0[j,i]*cfb[j,i]/M[j,i]/e/Gfast[j])**2 +Fex*dl[j,:,i]*cfb[j,i]/M[j,i]/e/Gfast[j]  +Fex*tau*dl[j,:,i+6]/M[j,i]/Rfb/e/gslow[i] )
               Nsc[j,:,i+6]=dl[j,:,i+6]/Rfb/M[j,i]/e/gslow[i]
     return dl0, dl, Nsc, dNsc, nkl, ksort, mpts_phase

def get_raw(R, type):
     #short routine for raw data.
     Raw = numpy.empty((30,598,12))
     for i in range(1,31,1):
         x = 'p'+str(i)
         K = R[i][:7176]
         L = numpy.reshape(K,(-1,12))
         Raw[i-1]=L
     return Raw
     

def  sph_int_laser(R, type, Times_all, it_Dark, it_las01, it_las02):

     #integrating spherical laser
    edata = R[type['LaserEnergy']]
    mx = max(edata)
    nn = Times_all.size

    ed = numpy.reshape(edata,(nn,12))
    erebin = rebin(ed, (nn,1))
    ed = numpy.reshape(erebin,(1,nn))[0,:]

    k0 = ed[it_Dark]
    k0 = k0[k0<-3]
    eavg0=numpy.average(k0)
    
    
    k1 = ed[it_las01]
    k1 = k1[k1>-3]
    eavg1=numpy.average(k1)-eavg0

    k2 = ed[it_las02]
    k2 = k2[k2>-3]
    eavg2=numpy.average(k2)-eavg0

    ev= ed-eavg0

    return ev, eavg0

def rebin(a, shape):
     #Compatibility function with MatLab, resizes arrays with averaging
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def gain(R,type, cgv):
    #Get Vdac Values of APDs
    VDac = numpy.empty((1,8*8))
    VT = numpy.empty((1,8))
    for i in range(1, 3+1, 1):
        VT = []
        for j in range(1, 8+1, 1):
            x = 'VDAC'+str(i)+'_'+str(j)
            L = R[type[x]]
            #K = numpy.array((L[0],L[3],L[1],L[4],L[2],L[5]))
            VT = numpy.append(VT,L)
        VDac = numpy.vstack([VDac,VT])
    Vdac = numpy.delete(VDac,0,0)
    Vdac = Vdac[:,0:60]
    Vdac = numpy.reshape(Vdac,(-1,6))
    #reorder=[0,3,1,4,2,5]
    vdac = numpy.zeros((30,6))
    vdac[:,0] = Vdac[:,0]
    vdac[:,1] = Vdac[:,3]
    vdac[:,2] = Vdac[:,1]
    vdac[:,3] = Vdac[:,4]
    vdac[:,4] = Vdac[:,2]
    vdac[:,5] = Vdac[:,5]

    #Gain Calculation
  
    a0 = numpy.reshape(cgv[:,:,0],(30,6))
    a1 = numpy.reshape(cgv[:,:,1],(30,6))
    a2 = numpy.reshape(cgv[:,:,2],(30,6))
    a3 = numpy.reshape(cgv[:,:,3],(30,6))

    a0 = a0 - vdac
    #weave THIS!
    gain = numpy.zeros((30,6))
    for i in range(0,6):
        for j in range(0,30):
            #gain[j,i]=numpy.roots((a0[j,i],a1[j,i],a2[j,i],a3[j,i]))
            r = numpy.roots((a3[j,i],a2[j,i],a1[j,i],a0[j,i]))
            r =  numpy.exp(abs(r[numpy.isreal(r)]))
            if r.size == 0:
                gain[j,i]= 1
            else:
                #multiple roots are bad! deal with them by finding value closest to 60 (usually around 35 or 70)
                gain[j,i] = r[min(range(len(r)), key=lambda z: abs(r[z]-60))]

    return gain

def F(x,a0,a1,a2,a3):
     #Compatability function with MatLab
    return a0 + a1*x+a2*x*x+a3*x*x*x

def times(R, type):

    nall_times_rack1 = R[type['nall_times_rack1']]
    nall_times_rack2 = R[type['nall_times_rack2']]
    nlaser1 = R[type['nlaser1']]
    nlaser2 = R[type['nlaser2']]

    time_rack1 = numpy.empty((1,1024,))
    time_rack1[:] = numpy.NAN
    time_rack2 = time_rack1
    Laser_times1 = time_rack1
    Laser_times2 = time_rack1

    if nall_times_rack1.size is not 0:
        time_rack1 = R[type['t1']][0:nall_times_rack1]

    if nall_times_rack2.size is not 0:
        time_rack2 = R[type['t2']][0:nall_times_rack2]

    nlaser = min(nlaser1,nlaser2)

    if nlaser1 > 2:
        Laser_times1 = R[type['L1times']][0:nlaser]

    if nlaser2 > 2:
        Laser_times2 = R[type['L2times']][0:nlaser]

        time_rack = time_rack1
    if max(time_rack2.size, time_rack1.size) is time_rack2.size:
        time_rack = time_rack2


    cond_laser1 = (Laser_times1 < max(time_rack)) & (Laser_times1 > min(time_rack))

    cond_laser2 = (Laser_times2 < max(time_rack)) & (Laser_times2 > min(time_rack))

    it_las01 = []
    for x in xrange(0, Laser_times1.size):
        #it_las01.append(min(abs(time_rack - Laser_times1[x])))
        TEMP1= abs(time_rack - Laser_times1[x])
        it_las01.append(numpy.argmin(TEMP1))

    it_las02 = []
    for x in xrange(0, Laser_times2.size):
        TEMP2= abs(time_rack - Laser_times2[x])
        it_las02.append(numpy.argmin(TEMP2))

        
    it_las01 = numpy.unique(it_las01)
    it_las02 = numpy.unique(it_las02)
    
    it_Dark = numpy.setdiff1d(numpy.arange(time_rack.size),numpy.union1d(it_las01[cond_laser1[:it_las01.size]],it_las02[cond_laser2[:it_las02.size]]))

    Times_Laser1 = time_rack[it_las01[cond_laser1[1:it_las01.size]]]
    Times_Laser2 = time_rack[it_las02[cond_laser2[0:it_las02.size]]]
    Times_dark = time_rack[it_Dark]
    Times_all = time_rack

    return Times_Laser1, Times_Laser2, Times_dark, Times_all, it_las01, it_las02, it_Dark, cond_laser1, cond_laser2

def shut_config(s1,s2):
    shut_config = 99
    if s1 == 'OPEN' and s2 == 'OPEN':
	    shut_config = 0
    if s1 == 'OPEN' and s2 == 'SHUT':
	    shut_config= 1
    if s1 == 'SHUT' and s2 == 'OPEN':
	    shut_config= 2
    if s1 == 'SHUT' and s2 == 'SHUT':
	    shut_config= 3
    return shut_config
    

def i_m(map):
     #Inverse dictionary mapping utility
	inv_map = {}
	for k, v in map.iteritems():
		inv_map[v] = inv_map.get(v, [])
	        inv_map[v].append(k)
	return inv_map

def f(regexStr,target):
    mo = re.search(regexStr,target)
    if not mo:
        print "No Match"
    else:
        print "Match!",mo.group()

def goto_shutters(shot):
    global tree1
    tree1 = MDSplus.Tree('nstx', shot, 'Readonly')   # MDSplus.tree.Tree
    n = tree1.getNode('ACTIVESPEC')          # MDSplus.treenode.TreeNode
    tree1.setDefault(n)
    return n
					    
def goto_rawdata(shot):
    global tree
    tree = MDSplus.Tree('nstx', shot, 'Readonly')   # MDSplus.tree.Tree
    n = tree.getNode('ACTIVESPEC')          # MDSplus.treenode.TreeNode
    tree.setDefault(n)
    p2 = tree.getNode('.mpts.shutter217')
    shutter1 = p2.getData().data()
    p2 = tree.getNode('.mpts.shutter218')
    shutter2 = p2.getData().data()
    n2 = tree.getNode('ACTIVESP_RAW')
    tree.setDefault(n2)
    n3 = tree.getNode('MPTS')
    tree.setDefault(n3)
    n = tree.getNode('RAWDATA')
    tree.setDefault(n)
    return shutter1, shutter2, n
					    
def goto_rawdata1(shot):
    global tree
    tree = MDSplus.Tree('nstx', shot, 'Readonly')   # MDSplus.tree.Tree
    n = tree.getNode('ACTIVESPEC')          # MDSplus.treenode.TreeNode
    tree.setDefault(n)
    n2 = tree.getNode('ACTIVESP_RAW')
    tree.setDefault(n2)
    n3 = tree.getNode('MPTS')
    tree.setDefault(n3)
    n = tree.getNode('RAWDATA')
    tree.setDefault(n)
    return n

def goto_sampledata(shot):
    global tree
    tree = MDSplus.Tree('nstx', shot, 'Readonly')   # MDSplus.tree.Tree
    n = tree.getNode('ACTIVESPEC')          # MDSplus.treenode.TreeNode
    tree.setDefault(n)
    n = tree.getNode('ACTIVESP_RAW')
    tree.setDefault(n)
    n = tree.getNode('MPTS')
    tree.setDefault(n)
    n = tree.getNode('RAWDATA')
    tree.setDefault(n)                             
    n = tree.getNode('TS_H908_01')
    tree.setDefault(n)
  #  n = tree.getNode('INPUT_01')
    return n

def goto_sampledata1(shot):
    global tree
    tree = MDSplus.Tree('nstx', shot, 'Readonly')   # MDSplus.tree.Tree
    n = tree.getNode('ACTIVESPEC')          # MDSplus.treenode.TreeNode
    tree.setDefault(n)
    n = tree.getNode('ACTIVESP_RAW')
    tree.setDefault(n)
    n = tree.getNode('MPTS')
    tree.setDefault(n)
    n = tree.getNode('RAWDATA')
    tree.setDefault(n)
    n = tree.getNode('TS_H908_01')
    tree.setDefault(n)
    n = tree.getNode('INPUT_01')
    return n
                                                  

def goto_p(n):
    global tree
    n = n.getParent()
    tree.setDefault(n)
    return n

def goto_c(obj,n):
    global tree
    n = tree.getNode(obj)
    tree.setDefault(n)
    return n

def d_shot():

    return 139047

def Print(obj):
    ans = str(obj)
    print ans
