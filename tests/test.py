import os
import scipy
from scipy import constants
from numpy import testing
from matplotlib import pyplot
from pywindow import *

nthreads = 8
ells = [0,2,4]

def test_2pcf_multi(typewin='global',los='midpoint',losn=0):
	"""
	s = scipy.linspace(2,100.,1)
	costheta = scipy.linspace(0.8,1.,201)
	distance = scipy.linspace(10.,11.,2)
	"""
	s = scipy.linspace(0.1,100.,100)
	costheta = scipy.linspace(-0.5,1.,201)
	distance = scipy.linspace(1000.,2000.,101)
	angular = 2.*scipy.ones(len(costheta))
	radial = scipy.ones(len(distance))
	
	pywindow = PyWindow()
	pywindow.set_2pcf_multi(s,costheta,angular,distance,radial,ells=ells,los=los,losn=losn,typewin=typewin,nthreads=nthreads)
	print pywindow.y

def test_3pcf_multi(typewin='global',los='midpoint',losn=0):
	
	if 'angular' in typewin:
		s = [scipy.linspace(0.1,100.,100)]*2
		costheta = scipy.linspace(0.8,1.,1000)
		distance = scipy.linspace(1000.,2000.,1000)
		angular = 2.*scipy.ones(len(costheta))
		radial = scipy.ones(len(distance))
	else:
		s = [scipy.linspace(0.1,100.,100)]*2
		costheta = [scipy.linspace(0.8,1.,1000)]*2
		distance = scipy.linspace(1000.,2000.,1000)
		angular = 2.*scipy.ones(map(len,costheta))
		radial = scipy.ones(len(distance))
	
	pywindow = PyWindow()
	pywindow.set_3pcf_multi(s,costheta,angular,distance,radial,ells=ells,los=los,losn=losn,typewin=typewin,nthreads=nthreads)
	print pywindow.y

def test_4pcf_multi(typewin='global-global',los='midpoint',losn=0):
	
	if 'angular' in typewin and typewin != 'angular-angular':
		s = [scipy.linspace(0.1,100.,100)]*2
		costheta = [scipy.linspace(0.8,1.,1000)]*2
		distance = scipy.linspace(1000.,2000.,1000)
		"""
		s = [scipy.linspace(0.1,100.,2)]*2
		costheta = [scipy.linspace(0.8,1.,10)]*2
		distance = scipy.linspace(100.,200.,10)
		"""
		angular = 2.*scipy.ones(map(len,costheta))
		radial = scipy.ones(len(distance))
	else:
		s = [scipy.linspace(0.1,100.,100)]*2
		#s = [scipy.linspace(100.,100.,1)]*2
		costheta = scipy.linspace(0.8,1.,1000)
		distance = scipy.linspace(1000.,2000.,1000)
		angular = 2.*scipy.ones(len(costheta))
		radial = scipy.ones(len(distance))
	
	pywindow = PyWindow()
	pywindow.set_4pcf_multi(s,costheta,angular,distance,radial,ells=ells,los=los,losn=losn,typewin=typewin,nthreads=nthreads)
	print pywindow.y

def test_interpol_bilin():
	x,y = 1.5,2.5
	fx = [1.,2.]
	fy = [2.,3.]
	f = [1.,2.,2.,3.]
	pywindow = PyWindow()
	print pywindow.interpol_bilin(x,y,fx,fy,f)

def plot_angular():
	costheta = scipy.linspace(-1.,1.,10)
	angular = 1. + scipy.cos(costheta)
	pywindow = PyWindow()
	pywindow.set_angular_selection(costheta,angular)
	costheta_ = costheta[1::2] + 0.05
	angular_ = pywindow.find_angular_selection(costheta_,interpol='poly')
	pyplot.plot(costheta,angular,color='k')
	pyplot.plot(costheta_,angular_,color='r')
	pyplot.show()

def test_verbosity():
	print('No output in between <<')
	pywindow = PyWindow()
	pywindow.set_verbosity('quiet')
	print('>>')
	
#test_2pcf_multi(typewin='global',los='endpoint')
#test_2pcf_multi(typewin='radial')
#test_2pcf_multi(typewin='angular')
test_3pcf_multi(typewin='global-dlos')
#test_3pcf_multi(typewin='angular')
#test_4pcf_multi(typewin='global-global')
#test_4pcf_multi(typewin='radial-radial')
#test_4pcf_multi(typewin='radial-global')
#test_4pcf_multi(typewin='angular-angular')
#test_4pcf_multi(typewin='angular-global')
#test_4pcf_multi(typewin='angular-radial')
#test_interpol_bilin()
#plot_angular()
#test_verbosity()

