import numpy as np
from scipy import fftpack
import logging
logger = logging.getLogger('ramses_pp.analysis.halos')

try:
	import numexpr as ne
	numexpr_available = True
except:
	numexpr_available = False

def power_spectrum_nd(input_array, box_dims):
	''' 
	Calculate the power spectrum of input_array and return it as an n-dimensional array,
	where n is the number of dimensions in input_array
	box_side is the size of the box in comoving Mpc. If this is set to None (default),
	the internal box size is used
	
	Parameters:
		* input_array (numpy array): the array to calculate the 
			power spectrum of. Can be of any dimensions.
		* box_dims = None (float or array-like): the dimensions of the 
			box. If this is None, the current box volume is used along all
			dimensions. If it is a float, this is taken as the box length
			along all dimensions. If it is an array-like, the elements are
			taken as the box length along each axis.
			[boxsize]*ndim
	
	Returns:
		The power spectrum in the same dimensions as the input array.		
	'''

	logger.info( 'Calculating power spectrum...')
	ft = fftpack.fftshift(fftpack.fftn(input_array.astype('float64')))
	power_spectrum = np.abs(ft)**2
	logger.info( '...done')

	# scale
	boxvol = np.product(map(float,box_dims))
	pixelsize = boxvol/(np.product(input_array.shape))
	power_spectrum *= pixelsize**2/boxvol
	
	return power_spectrum


def cross_power_spectrum_nd(input_array1, input_array2, box_dims):
	''' 
	Calculate the cross power spectrum two arrays and return it as an n-dimensional array,
	where n is the number of dimensions in input_array
	box_side is the size of the box in comoving Mpc. If this is set to None (default),
	the internal box size is used
	
	Parameters:
		* input_array1 (numpy array): the first array to calculate the 
			power spectrum of. Can be of any dimensions.
		* input_array2 (numpy array): the second array. Must have same 
			dimensions as input_array1.
		* box_dims = None (float or array-like): the dimensions of the 
			box. If this is None, the current box volume is used along all
			dimensions. If it is a float, this is taken as the box length
			along all dimensions. If it is an array-like, the elements are
			taken as the box length along each axis.
	
	Returns:
		The cross power spectrum in the same dimensions as the input arrays.
		
	TODO:
		Also return k values.
	'''

	assert(input_array1.shape == input_array2.shape)

	logger.info( 'Calculating power spectrum...')
	ft1 = fftpack.fftshift(fftpack.fftn(input_array1.astype('float64')))
	ft2 = fftpack.fftshift(fftpack.fftn(input_array2.astype('float64')))
	power_spectrum = np.real(ft1)*np.real(ft2)+np.imag(ft1)*np.imag(ft2)
	logger.info( '...done')

	# scale
	#boxvol = float(box_side)**len(input_array1.shape)
	boxvol = np.product(map(float,box_dims))
	pixelsize = boxvol/(np.product(map(float,input_array1.shape)))
	power_spectrum *= pixelsize**2/boxvol

	return power_spectrum


def radial_average(input_array, box_dims, kbins=10):
	'''
	Radially average data. Mostly for internal use.
	
	Parameters: 
		* input_array (numpy array): the data array
		* box_dims = None (float or array-like): the dimensions of the 
			box. If this is None, the current box volume is used along all
			dimensions. If it is a float, this is taken as the box length
			along all dimensions. If it is an array-like, the elements are
			taken as the box length along each axis.
		* kbins = 10 (integer or array-like): The number of bins,
			or a list containing the bin edges. If an integer is given, the bins
			are logarithmically spaced.
			
	Returns:
		A tuple with (data, bins, n_modes), where data is an array with the 
		averaged data, bins is an array with the bin centers and n_modes is the 
		number of modes in each bin

	'''

	k_comp, k = _get_k(input_array, box_dims)

	kbins = _get_kbins(kbins, box_dims, k)
	
	#Bin the data
	logger.info('Binning data...')
	dk = (kbins[1:]-kbins[:-1])/2.
	#Total power in each bin
	outdata = np.histogram(k.flatten(), bins=kbins,
						weights = input_array.flatten())[0]
	#Number of modes in each bin
	n_modes = np.histogram(k.flatten(), bins=kbins)[0].astype('float')
	outdata /= n_modes
	
	return outdata, kbins[:-1]+dk, n_modes
	

def power_spectrum_1d(input_array_nd, box_dims, kbins=100, return_n_modes=False):
	''' Calculate the spherically averaged power spectrum of an array 
	and return it as a one-dimensional array.
	
	Parameters: 
		* input_array_nd (numpy array): the data array
		* kbins = 100 (integer or array-like): The number of bins,
			or a list containing the bin edges. If an integer is given, the bins
			are logarithmically spaced.
		* box_dims = None (float or array-like): the dimensions of the 
			box. If this is None, the current box volume is used along all
			dimensions. If it is a float, this is taken as the box length
			along all dimensions. If it is an array-like, the elements are
			taken as the box length along each axis.
		* return_n_modes = False (bool): if true, also return the
			number of modes in each bin
			
	Returns: 
		A tuple with (Pk, bins), where Pk is an array with the 
		power spectrum and bins is an array with the k bin centers.
	'''

	input_array = power_spectrum_nd(input_array_nd, box_dims=box_dims)	

	ps, bins, n_modes = radial_average(input_array, kbins=kbins, box_dims=box_dims)
	if return_n_modes:
		return ps, bins, n_modes
	return ps, bins


def cross_power_spectrum_1d(input_array1_nd, input_array2_nd, box_dims, kbins=100, return_n_modes=False):
	''' Calculate the spherically averaged cross power spectrum of two arrays 
	and return it as a one-dimensional array.
	
	Parameters: 
		* input_array1_nd (numpy array): the first data array
		* input_array2_nd (numpy array): the second data array
		* kbins = 100 (integer or array-like): The number of bins,
			or a list containing the bin edges. If an integer is given, the bins
			are logarithmically spaced.
		* box_dims = None (float or array-like): the dimensions of the 
			box. If this is None, the current box volume is used along all
			dimensions. If it is a float, this is taken as the box length
			along all dimensions. If it is an array-like, the elements are
			taken as the box length along each axis.
		* return_n_modes = False (bool): if true, also return the
			number of modes in each bin
			
	Returns: 
		A tuple with (Pk, bins), where Pk is an array with the 
		cross power spectrum and bins is an array with the k bin centers.
	'''

	input_array = cross_power_spectrum_nd(input_array1_nd, input_array2_nd, box_dims=box_dims)	

	ps, bins, n_modes = radial_average(input_array, kbins=kbins, box_dims = box_dims)
	if return_n_modes:
		return ps, bins, n_modes
	return ps, bins

#Some methods for internal use

def _get_k(input_array, box_dims):
	'''
	Get the k values for input array with given dimensions.
	Return k components and magnitudes.
	For internal use.
	'''
	dim = len(input_array.shape)
	if dim == 1:
		x = np.arange(len(input_array))
		center = x.max()/2.
		kx = 2.*np.pi*(x-center)/box_dims[0]
		return [kx], kx
	elif dim == 2:
		x,y = np.indices(input_array.shape, dtype='int32')
		center = np.array([(x.max()-x.min())/2, (y.max()-y.min())/2])
		kx = 2.*np.pi * (x-center[0])/box_dims[0]
		ky = 2.*np.pi * (y-center[1])/box_dims[1]
		k = np.sqrt(kx**2 + ky**2)
		return [kx, ky], k
	elif dim == 3:
		x,y,z = np.indices(input_array.shape, dtype='int32')
		center = np.array([(x.max()-x.min())/2, (y.max()-y.min())/2, \
						(z.max()-z.min())/2])
		kx = 2.*np.pi * (x-center[0])/box_dims[0]
		ky = 2.*np.pi * (y-center[1])/box_dims[1]
		kz = 2.*np.pi * (z-center[2])/box_dims[2]

		k = get_eval()('(kx**2 + ky**2 + kz**2 )**(1./2.)') 		
		return [kx,ky,kz], k


def _get_kbins(kbins, box_dims, k):
	'''
	Make a list of bin edges if kbins is an integer,
	otherwise return it as it is.
	'''
	if isinstance(kbins,int):
		kmin = 2.*np.pi/min(box_dims)
		kbins = 10**np.linspace(np.log10(kmin), np.log10(k.max()), kbins+1)
	return kbins


def _get_nonzero_idx(ps_shape, los_axis):
	'''
	Get the indices where k_perp != 0
	'''
	x,y,z = np.indices(ps_shape)
	if los_axis == 0:
		zero_idx = (y == ps_shape[1]/2)*(z == ps_shape[2]/2)
	elif los_axis == 1:
		zero_idx = (x == ps_shape[0]/2)*(z == ps_shape[2]/2)
	else:
		zero_idx = (x == ps_shape[0]/2)*(y == ps_shape[1]/2)
	good_idx = np.invert(zero_idx)
	return good_idx

def get_eval():
	'''
	Evaluate an expression using numexpr if
	available. For internal use.
	''' 
	if numexpr_available:
		return ne.evaluate
	return eval
