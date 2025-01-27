#######################################################################
#
# Functions used by plotting figures
# Written by Ho Sang (Leon) Chan at 2023
# The Chinese University of Hong Kong and University of Colorado
#
#######################################################################

#import packages#
import math
import scipy
import numpy as np
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from matplotlib.colors import TwoSlopeNorm

#######################################################################

# define domain #
domain = [-4700, 4700, -4700, 4700]
smalldomain = [-4700, 4700, -4700, 4700]
stream = [0, 4700, -4700, 4700]
smallstream = [0, 4700, -4700, 4700]

bkey = ['bx','by','bz','bfield']


#######################################################################
# define function #

# plot contour #
def contour(x1_in,y1_in,z1_in,vx1_in=None,vy1_in=None,z2_in=None,rbh=None,fname=None,time_in=None):

	# set figure #
	fig, ax = plt.subplots(1, 1)
	plt.title('Time = '+str(time_in), size=15)
	cbar_ax = fig.add_axes([0.05, 0.05, 0.9, 0.01])

	# set label #
	ax.set_xlabel('x1', size=15)
	ax.set_ylabel('x2', size=15)
	
	# log plot #
	zmin = np.nanmin(z1_in[z1_in != -np.inf])
	zmax = np.nanmax(z1_in[z1_in != np.inf])
	norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
	tick = np.linspace(zmin, zmax, 10,  endpoint=True)
	cmap = 'plasma'
		
  #plot#
	ax.pcolormesh(x1_in, y1_in, z1_in, 100, norm=norm, cmap=cmap, extend='min')

	# ax.contourf(x1_in, y1_in, z1_in, 100, norm=norm, cmap=cmap, extend='min')
	if(z2_in is not None):
		ax.contourf(-x1_in, y1_in, z2_in, 100, norm=norm, cmap=cmap, extend='min')

	# streamline #
	if((vx1_in is not None) and (vy1_in is not None)):
		ax.streamplot(x1_in,y1_in,vx1_in,vy1_in,density=1,linewidth=1,arrowsize=1,color='black',broken_streamlines=True)

	#domain#
	time_in = 8000

	if any(keyword in fname for keyword in bkey):
		ax.set_xlim(smalldomain[0], smalldomain[1])
		ax.set_ylim(smalldomain[2], smalldomain[3])
	else:
		ax.set_xlim(domain[0], domain[1])
		ax.set_ylim(domain[2], domain[3])

  #colorbar#
	cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, ticks=tick,orientation='horizontal')

	# further plot parameters # 
	ax.tick_params(axis='both', labelsize=15)
	ax.set_aspect('equal')
	plt.subplots_adjust(bottom=0.2)
	plt.savefig(str(fname)+'.png',dpi=500)
	plt.clf()
	plt.close()

#######################################################################
# define function #

# plot contour #
def contour_interpolate_stream(x1_in,y1_in,z1_in,xaxis_in=None,yaxis_in=None,vx1_in=None,vy1_in=None,z2_in=None,rbh=None,fname=None,time_in=None):

	# set figure #
	fig, ax = plt.subplots(1, 1)
	plt.title('Time = '+str(time_in), size=15)
	cbar_ax = fig.add_axes([0.05, 0.05, 0.9, 0.01])

	# set label #
	ax.set_xlabel('x1', size=15)
	ax.set_ylabel('x2', size=15)
	
	# log plot #
	zmin = np.nanmin(z1_in[z1_in != -np.inf])
	zmax = np.nanmax(z1_in[z1_in != np.inf])
	norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
	tick = np.linspace(zmin, zmax, 10,  endpoint=True)
	cmap = 'plasma'
		
  #plot#
	ax.pcolormesh(x1_in, y1_in, z1_in, norm=norm, cmap=cmap)

	# ax.contourf(x1_in, y1_in, z1_in, 100, norm=norm, cmap=cmap, extend='min')
	if(z2_in is not None):
		ax.contourf(-x1_in, y1_in, z2_in, 100, norm=norm, cmap=cmap, extend='min')

	# streamline #
	if((xaxis_in is not None) and (yaxis_in is not None) and (vx1_in is not None) and (vy1_in is not None)):
		x_cart, y_cart, bx_stream, by_stream = interpolate_stream_rtheta(xaxis_in, yaxis_in, vx1_in, vy1_in, time_in)
		ax.streamplot(x_cart, y_cart, bx_stream, by_stream, density=3, color='black', linewidth=0.5, arrowsize=0.5)

	#domain#
	time_in = 8000
	if any(keyword in fname for keyword in bkey):
		ax.set_xlim(smalldomain[0], smalldomain[1])
		ax.set_ylim(smalldomain[2], smalldomain[3])
	else:
		ax.set_xlim(domain[0], domain[1])
		ax.set_ylim(domain[2], domain[3])

  	#colorbar#
	cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, ticks=tick,orientation='horizontal')

	# further plot parameters # 
	ax.tick_params(axis='both', labelsize=15)
	ax.set_aspect('equal')
	plt.subplots_adjust(bottom=0.2)
	plt.savefig(str(fname)+'.png',dpi=500)
	plt.clf()
	plt.close()

#######################################################################
# define fuction #

# find vector potential #
def interpolate_stream_rtheta(xaxis, yaxis, bx, by, time_in): 

	# interpolation #
	time_in = 8000
	if time_in > 27000:
		x_cart = np.linspace(smallstream[0], smallstream[1], 100)
		y_cart = np.linspace(smallstream[2], smallstream[3], 100)
	else:
		x_cart = np.linspace(stream[0], stream[1], 100)
		y_cart = np.linspace(stream[2], stream[3], 100)

	# sine and cosine #
	x_cart, y_cart = np.meshgrid(x_cart, y_cart)
	radius = np.sqrt(x_cart**2 + y_cart**2)
	theta = np.arccos(y_cart/radius) 
	sine = x_cart/radius
	cosine = y_cart/radius

	# interpolation #
	fx = scipy.interpolate.RegularGridInterpolator((xaxis, yaxis), bx, method='pchip', bounds_error=False, fill_value=None)
	fy = scipy.interpolate.RegularGridInterpolator((xaxis, yaxis), by, method='pchip', bounds_error=False, fill_value=None)

	# coordinate transformation #
	br_stream = fx((radius, theta))
	bth_stream = fy((radius, theta))
	bx_stream = br_stream*sine + bth_stream*cosine
	by_stream = br_stream*cosine - bth_stream*sine

	# return #
	return x_cart, y_cart, bx_stream, by_stream

#######################################################################
