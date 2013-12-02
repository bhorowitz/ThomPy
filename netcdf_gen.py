from config import *
from scipy.io import netcdf
from datetime import date
import os

version = '0.1a'

today = date.today()
t = today.timetuple()

su = netcdf.netcdf_file('setup.netcdf', 'w')

su.history = 'Version: ' + version + ", Made on " + str(t[2]) + '/' + str(t[1]) + '/' + str(t[0])

su.author = 'Benjamin Horowitz at horowitz.ben@gmail.com'

su.cF = cF

su.qt12_lambda = qt12_lambda
su.qt12d_lambda = qt12d_lambda

su.qt12_fb = qt12_fb
su.qt12_dfb = qt12_dfb


su.qt3_lambda = qt3_lambda
su.qt3d_lambda = qt3d_lambda

su.qt3_fb = qt3_fb
su.qt3_dfb = qt3_dfb

su.qt12_wv_all = qt12_wv_all
su.qt3_wv_all = qt3_wv_all

su.nwv_ph12 = nwv_ph12
su.nwv_ph3 = nwv_ph3

su.nwvl5 = nwvl5
su.nrad5 =nrad5
su.qt_lambda = qt_lambda
su.dqt_lambda = dqt_lambda
su.wv_all = wv_all
su.fb = fb
su.dfb = dfb

su.cg = cg

su.cf2 = cf2
su.cf = cf
su.cfb = cfb
su.dcfb = dcfb
su.theta = theta

cF = numpy.loadtxt('./setup/cF')
cF = numpy.reshape(cF,(30,6,2))

su.cF = cF

cgv12 = numpy.loadtxt('./setup/cgv12')
cgv12 = cgv12.reshape(20,6,4)
cgv3 = numpy.loadtxt('./setup/cgv3')
cgv3 = cgv3.reshape(10,4,4)
cgv = numpy.zeros([30,6,4])
cgv[0:20,:,:]=cgv12
cgv[20:30,0:4,:]=cgv3

su.cgv = cgv

su.close()
