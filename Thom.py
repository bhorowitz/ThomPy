from ts import *
import numpy
from MDSplus import *
import csv
import os
from scipy.io import netcdf
import time

"""
NAME
     PyThom: NSTX Thompson Scattering Analysis Code

PURPOSE
     Analyze Output of NSTX Thomson Scattering Diagnostic

STATUS
     Under development and minimally operational. Modules are in ts.py.
     
RUN
     module load nstx/python
     Python PyThom.py [SHOT NUMBER]

USER INPUTS
     Shot Number

OTHER INPUTS
    - setup.netcdf (see utils, netcdf_gen.py)
    - NSTX MDSplus structure
    - ./ss/ folder (a crutch till new datastream)

OTHER PARAMETERS
    - # of processors on system (nproc)

OUTPUTS
    [shot number].npz file which contains timings and \
    profile information
    - timing information is a (X) length array
    - profile information is a (30,2,X) array
    
LOG
    20/06/13: Research Begins
    03/07/13: Prototyping Begins
    22/07/13: Importing Functional(?)
    24/07/13: Scattered Photo-electrons Functional
    03/08/13: Fitting Functional
    03/13/13: MOC Operational

CONTACT
    Ben Horowitz (Yale University)
    - Benjamin.a.Horowitz@yale.edu
    - Horowitz.Ben@gmail.com
"""
tic = time.time()
shot = int(sys.argv[1])
good = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
nproc = 8 #Number of system processors

#configuration data
#generated with netcdf_gen.py
conf = netcdf.netcdf_file('setup.netcdf', 'r')

#read in list of paths

ins = open("paths", "r") # file that holds the name of all MDS routes of interest and their natural name
path = [] # path list that will be formed
nb = {} # numerical name to path name
type = {} # key from "natural name" to numerical name in data array

for x,line in enumerate(ins):
    columns = line.split()
    path.append(str(columns[0]))
    nb[x] = str(columns[0])
    type[str(columns[1])] = x
ins.close()

# navigate to data
tree = MDSplus.Tree('nstx', shot, 'Readonly') 

# get plasma current

p = tree.getNode('\engineering::ip1')
ip1 = p.getData().data()
p = tree.getNode('\engineering::ip2')
ip2 = p.getData().data()

ip = ip2/(1e6)
#shutters, s1 and s2 are either OPEN or CLOSED
s1, s2, r = ts.goto_rawdata(shot)

shut_config = ts.shut_config(s1,s2) # assigns numerical name (0,1,2,3,99) for shutter configuration, used in nfactor querry

R = [[]]
for x in path:
    #reads in all data to numpy array structure
    p = tree.getNode(x)
    ll = p.getData().data()
    R.append(ll)

Times_Laser1, Times_Laser2, Times_dark, Times_all, it_las01, it_las02, it_Dark, cond_laser1, cond_laser2  = ts.times(R,type) #gets timing information for shot

data = ts.get_raw(R,type) #loads raw data into array

#anomolout entries in it_las01/02...
it_las01 = numpy.delete(it_las01,-1)
it_las01 = numpy.delete(it_las01,0)
it_las02 = numpy.delete(it_las02,-1)
it_las02 = numpy.delete(it_las02,0)

cgv = numpy.reshape(conf.cgv,(30,6,4))
gain = ts.gain(R,type,cgv) #gets gain for shot

#integrating spherical laser
ev_sph,eavg_sph = ts.sph_int_laser(R, type, Times_all,it_Dark, it_las01, it_las02)

nchannels=30
nshelf=30

e=1.6022e-19  #coulombs
Rfb=5e4       # ohms
nmux=12
tau=50e-9

Gfast = numpy.zeros((1,30))+10
Gfast = Gfast[0]

#anomolous entries at begining of it_las01/02 (2,5)
cfb = conf.cfb #to be passed into photoelectrons

cF = conf.cF
cF = numpy.reshape(cF,(30,6,2))
cfb = numpy.reshape(cfb,(30,6))
dl0, dl, Nsc, dNsc, nkl, ksort, mpts_phase = ts.photoelectrons(it_las01, it_las02, data, nchannels,nmux,gain,Times_all,cfb,Gfast, cF)

#select which wavelengths for initial guess

Rguess = ts.Rguess(Nsc,nchannels,nkl)

tip = numpy.loadtxt('./ss/tip')
dip = numpy.loadtxt('./ss/dip')

tl = Times_all[ksort]

kdrop = ((dip > -40) & (dip > 0))

t_maxip = tip[numpy.argmax(ip)]

tmin=0.005
tmax=10.0
ip_min=0.025

i1 = (tip > tmin)
i2 = (ip > ip_min)
kip = numpy.all([i1,i2], axis=0)
kt = numpy.all([(tl > tmin),(tl < tip[kip].max()),(tl < tmax)],axis=0)

kt = numpy.arange(kt.size)[kt]
kt  = numpy.delete(kt, -1)
kt = numpy.delete(kt, -1)


#brief sanatization
time = tl[kt]

nkt = time.size
valid = numpy.zeros((nchannels,nkt))
sat_ok = numpy.zeros((nchannels,nkt))

for i in range(0,nkt):
    for j in range(0,nchannels):
        ksat_=(dl[j,kt[i],:] > 9.7)
        if ksat_.max() == False:
            sat_ok[j,i]=1

#nfactor_30, requires ramray_select, shutter select, qt_cal

nfactor_20 = numpy.loadtxt('./ss/nfactor_20_'+str(shut_config),delimiter=',')
nfactor_30 = numpy.loadtxt('./ss/nfactor_30_'+str(shut_config),delimiter=',')

    #nfactor can be pre-fabricated!!

Elaser_R_st = numpy.loadtxt('./ss/Elaser_R_st')
Elaser_R_st_20 = numpy.loadtxt('./ss/Elaser_R_st_20')
Elaser_R_sph = numpy.loadtxt('./ss/Elaser_R_sph')
Elaser_R_sph_20 = numpy.loadtxt('./ss/Elaser_R_sph_20')

ev=ev_sph
iElaser_R=Elaser_R_sph
ksort1 = numpy.loadtxt('./ss/ksort', delimiter=',')

#qt = numpy.loadtxt('./QT', delimiter=',')
#QT = numpy.loadtxt('./setup/QT1', delimiter=',')

wv_all =numpy.reshape(conf.wv_all,[30,330])
cg = numpy.reshape(conf.cg,[30,5])
fb = numpy.reshape(conf.fb,[30,330])
theta = conf.theta

#a_tot_above = gtefit(wv_all,nshelf,nkt,Nsc,dNsc,cg,Rguess,QT,fb,theta,kt,ishelf)
QT = numpy.reshape(conf.QT,[5,330])
#gte def

kr = time.time()
print kr-tic

import sys
sys.path.append("/u/bhorowit/pprocess-0.5")
import pprocess

results = pprocess.Map(limit=nproc, reuse=1)
parallel_function = results.manage(pprocess.MakeReusable(gtefitP))
tic=time.time()
[parallel_function(wv_all,nshelf,nkt,Nsc,dNsc,cg,Rguess,QT,fb,theta,kt,ishelf) for ishelf in range(0,nshelf-2)]

a = numpy.zeros([nshelf,2,nkt])
for i in range(0,nshelf-2):
    a[i,:,:]=results[i]

print time.time()-tic

timings = tl[kt]

numpy.savez(str(shot),a,timings,Nsc)
