#!/usr/bin/python
import subprocess
import numpy as np
import micro_module
from micro_module import *

def make_stars():
	#arguments to subprocess must be strings; they are converted to doubles in the c program
	#compile simu.c script from Sebastiano (makefile must be in the same folder as simu.c)
	subprocess.call('make')

	#length of each side of image crop
	trim_size = 100

	#position of total image crop
	imcrop_x = 800
	imcrop_y = 3000

	#convert input values to strings for c program call
	trim_size_s = str(trim_size)
	imcrop_x_s = str(imcrop_x)
	imcrop_y_s = str(imcrop_y)

	#define ranges for the star's x-y position and flux (check flux range)
	x_pos_range = np.linspace(0, trim_size, trim_size+1)
	y_pos_range = np.linspace(0, trim_size, trim_size+1)
	flux_range = np.linspace(1000, 10000, 1001)

	#i denotes the number of stars in images that will be created
	#there is currently one star produced per image
	for i in range(10):
		print "\n------------------Iteration Number: ", i+1, " ------------------"
		#injected star position
		x_star_position =  str(np.random.choice(x_pos_range))
		y_star_position = str(np.random.choice(y_pos_range))
		print "Injected Star Parameters:"
		print 'x = ', x_star_position
		print 'y = ', y_star_position

		#injected star flux; flux conversion needed (not sure what value 1000 actually converts to)
		injs_flux = str(np.random.choice(flux_range))
		print "flux = ", injs_flux
		star_output_imname = 'injection_output_{}.fits'.format(i+1)

		#call injection simulation c programs
		#run simu script from Sebastiano
		subprocess.call(['./simu', 'images/UKIRT_21_1_0001.fits', 'images/snap_pre_UKIRT_21_1_0001.fits.fits', \
			x_star_position, y_star_position, injs_flux, trim_size_s, star_output_imname, imcrop_x_s, imcrop_y_s])
		print "\n-----------------------------------------------------------"


times = np.linspace(0., 100., 101)
u0 = 1.
tE = 50.
t0 = 0.
Ftot = 1.
fb = 1.
values = fluxCurve(times,u0,tE,t0,Ftot,fb)
print len(values)
plt.plot(times, values)
plt.plot((times*-1.), values)
plt.show()
