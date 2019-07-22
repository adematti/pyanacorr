import os
import ctypes
import scipy
from scipy import constants
from numpy import ctypeslib

def new(tab,dtype=ctypes.c_double,copy=False):
	if copy: return scipy.asarray(tab,dtype=dtype).flatten()
	return tab.astype(dtype).flatten()

class PyWindow(object):

	C_TYPE = ctypes.c_double
	PATH_WINDOW = os.path.join(os.path.dirname(os.path.realpath(__file__)),'window.so')

	def __init__(self):

		self.window = ctypes.CDLL(self.PATH_WINDOW,mode=ctypes.RTLD_LOCAL)
		self.clear()

	def set_verbosity(self,mode='info'):
		self.window.set_verbosity.argtypes = (ctypes.c_char_p,)
		self.window.set_verbosity(mode)
		self._verbose = mode
	
	def set_pole(self,num=1,ells=[0,2,4,6,8,10,12]):

		nells = len(ells)
		if (scipy.mod(ells,2)==0).all():
			ells = scipy.arange(nells)*2
			multitype = 'even'
		elif (scipy.mod(ells,2)==1).all():
			ells = 1+scipy.arange(nells)*2
			multitype = 'odd'
		else:
			ells = scipy.arange(nells)
			multitype = 'all'
		
		self.window.set_pole.argtypes = (ctypes.c_size_t,ctypes.c_char_p,ctypes.c_size_t)
		self.window.set_pole(num,multitype,nells)
		self._ells[num] = ells.tolist()

		return list(self._ells[num])
	
	def set_los(self,num=1,los='midpoint',n=0):
		
		self.window.set_los.argtypes = (ctypes.c_size_t,ctypes.c_char_p,ctypes.c_size_t)
		self.window.set_los(num,los,n)

		self._los[num] = (los,n)

		return tuple(self._los[num])
	
	def set_n_poles(self,ells,n=1):
		if scipy.isscalar(ells[0]): ells = [ells]*n
		ells = [self.set_pole(ill+1,ells=ells[ill]) for ill in range(n)]
		return ells
	
	def set_n_los(self,los,losn=0,n=1):
		if isinstance(los,(str,unicode)): los = [los]*n
		if scipy.isscalar(losn): losn = [losn]*n
		los = [self.set_los(ilos+1,los=los[ilos],n=losn[ilos]) for ilos in range(n)]
		return los

	def find_angular_selection(self,x,y=None,interpol='lin'):

		x = scipy.asarray(x)
		if len(self._costheta) == 2:
			y = scipy.asarray(y)
			self.window.find_angular_selection_2d.argtypes = (self.C_TYPE,self.C_TYPE)
			self.window.find_angular_selection_2d.restype = self.C_TYPE
			toret = [[self.window.find_angular_selection_2d(x_,y_) for y_ in y] for x_ in x]
		else:
			interpol = 1 if interpol=='poly' else 0
			self.window.find_angular_selection_1d.argtypes = (self.C_TYPE,ctypes.c_size_t)
			self.window.find_angular_selection_1d.restype = self.C_TYPE
			toret = [self.window.find_angular_selection_1d(x_,interpol) for x_ in x]

		return scipy.array(toret)

	def set_2pcf_multi(self,s,costheta,angular,distance,radial,ells=[0,2,4,6,8,10,12],los='midpoint',losn=0,typewin='global',nthreads=8):

		self.ells = self.set_pole(1,ells=ells)
		self.los = self.set_los(1,los=los,n=losn)
		self.set_angular_selection(costheta,angular)
		self.set_n_radial_selection(distance,radial,n=2)
		self.set_window(s)
		assert typewin in ['global','radial','angular']
		
		self.window.run_2pcf_multi.argtypes = (ctypes.c_char_p,ctypes.c_size_t,)
		self.window.run_2pcf_multi(typewin,nthreads)
	
	def set_3pcf_multi(self,s,costheta,angular,distance,radial,ells=[0,2,4,6,8,10,12],los='midpoint',losn=0,typewin='global',nthreads=8):

		self.ells = self.set_n_poles(ells,n=2)
		self.los = self.set_n_los(los,losn,n=2)
		self.set_angular_selection(costheta,angular)
		self.set_n_radial_selection(distance,radial,n=3)
		self.set_window(s)
		assert typewin in ['global','radial','angular']
		
		self.window.run_3pcf_multi.argtypes = (ctypes.c_char_p,ctypes.c_size_t,)
		self.window.run_3pcf_multi(typewin,nthreads)

	def set_4pcf_multi(self,s,costheta,angular,distance,radial,ells=[0,2,4,6,8,10,12],los='midpoint',losn=0,typewin='global-global',nthreads=8):

		self.ells = self.set_n_poles(ells,n=2)
		self.los = self.set_n_los(los,losn,n=2)
		self.set_angular_selection(costheta,angular)
		self.set_n_radial_selection(distance,radial,n=4)
		self.set_window(s)
		assert typewin in ['global-global','radial-radial','radial-global','angular-angular','angular-global','angular-radial']
		
		self.window.run_4pcf_multi.argtypes = (ctypes.c_char_p,ctypes.c_size_t,)
		self.window.run_4pcf_multi(typewin,nthreads)
		
	def set_n_radial_selection(self,distances,radials,n=1):

		if not isinstance(distances,list): distances = [distances]
		if not isinstance(radials,list): radials = [radials]
		distances += [distances[-1]]*(n-len(distances))
		radials += [radials[-1]]*(n-len(radials))

		for num,(distance,radial) in enumerate(zip(distances,radials)):
			self.set_radial_selection(num+1,distance,radial)

	def set_radial_selection(self,num,distance,radial,copy=False):

		size = len(distance)
		
		self._distance[num] = new(distance,dtype=self.C_TYPE,copy=copy)
		self._radial[num] = new(radial,dtype=self.C_TYPE,copy=copy)

		typedistance = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=size)
		typeradial = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=size)
		
		self.window.set_radial_selection.argtypes = (ctypes.c_size_t,typedistance,typeradial,ctypes.c_size_t)
		self.window.set_radial_selection(num,self._distance[num],self._radial[num],size)

	def set_angular_selection(self,costheta,angular,copy=False):
		
		if scipy.isscalar(costheta[0]): costheta = [costheta]
		shapeangular = angular.shape
		self._shape_angular = scipy.array(shapeangular,dtype=ctypes.c_size_t)
		
		self._costheta = [new(x_,dtype=self.C_TYPE,copy=copy) for x_ in costheta]
		self.__costheta = scipy.concatenate(self._costheta)
		shape = tuple(map(len,self._costheta))
		self._angular = new(angular,dtype=self.C_TYPE,copy=copy)

		typecostheta = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=scipy.sum(shape))
		typeangular = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=scipy.prod(shape))
		
		self.window.set_angular_selection.argtypes = (typecostheta,typeangular,ctypeslib.ndpointer(dtype=ctypes.c_size_t,shape=len(shape)),ctypes.c_size_t)
		self.window.set_angular_selection(self.__costheta,self._angular,self._shape_angular,len(shapeangular))
		
		self._angular.shape = shapeangular
		
	def set_window(self,x,nells=None):

		if scipy.isscalar(x[0]): x = [x]
		self.x = [new(x_,dtype=self.C_TYPE) for x_ in x]
		self._x = scipy.concatenate(self.x)
		shapex = tuple(map(len,self.x))
		self._shape_window = scipy.array(shapex,dtype=ctypes.c_size_t)
		
		if nells is None: nells = tuple(len(self._ells[num+1]) for num in range(len(shapex)))
		shapey = shapex + nells
		self.y = scipy.zeros(shapey,dtype=self.C_TYPE).flatten()
		
		typex = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=scipy.sum(shapex))
		typey = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=scipy.prod(shapey))
		
		self.window.set_window.argtypes = (typex,typey,ctypeslib.ndpointer(dtype=ctypes.c_size_t,shape=len(shapex)),ctypes.c_size_t)
		self.window.set_window(self._x,self.y,self._shape_window,len(shapex))

		self.y.shape = shapey
	
	def interpol_bilin(self,x,y,fx,fy,f):
		
		fx = new(fx,dtype=self.C_TYPE,copy=True)
		fy = new(fy,dtype=self.C_TYPE,copy=True)
		f = new(f,dtype=self.C_TYPE,copy=True)
		typep = ctypeslib.ndpointer(dtype=self.C_TYPE,shape=2)
		self.window.interpol_bilin.argtypes = (self.C_TYPE,self.C_TYPE,typep,typep,ctypeslib.ndpointer(dtype=self.C_TYPE,shape=4))
		self.window.interpol_bilin.restype = self.C_TYPE
		
		return self.window.interpol_bilin(x,y,fx,fy,f)

	def clear(self):
		self._distance = {}
		self._radial = {}
		self._ells = {}
		self._los = {}
		self.window.clear_radial_selections()
		self.window.clear_poles()
		self.set_verbosity()
		
