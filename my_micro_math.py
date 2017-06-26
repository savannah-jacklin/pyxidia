import numpy as np
import matplotlib.pyplot as plt

#change these values to match ukirt inputs, but mimcs sumi's injection strategy
#parameters for artificial event injection
#does the sample size matter in the linspace for statistics beyond what we got going here

def random_injection_parameters():
	#impact parameter
	u_0_range = np.linspace(0.0, 1.5, 1000)

	#time of peak magnification
	t_0_range = np.linspace(2453824, 2454420, 10000)

	#einstein radius crossing time
	t_E_range = np.linspace(0.1, 250, 1000)

	#source magnitude
	I_s_range = np.linspace(14.25, 21.0, 250)

	count = 0
	while (count <=10):
		#draw random input values from ranges, to make 'count' number of light curve sets; adjust for UKIRT?
		u_0 = np.random.choice(u_0_range)
		t_0 = np.random.choice(t_0_range)
		t_E = np.random.choice(t_E_range)
		I_s = np.random.choice(I_s_range)
	
		print 'Iteration ', count
		print 't_0: ', t_0
		print 'u_0: ', u_0
		print 't_E: ', t_E
		print 'I_s: ', I_s
		print '\n'
		count+=1


#input parameters
G = 6.67259e-8
#source mass?
M = 1.99e33
c = 2.99792458e10
D_s = 2.4685e22
D_l = 1.2343e22

u_min = 2.5#1.1
t_0 = 0.
t_E = 45.
theta = 0.001#angular separation of images of the source and lens on the sky
kappa = (4.*G)/c**2.

def einstein_radius(D_s, D_l, G, M, c):
	theta_E = np.sqrt(((4.*G*M)/c**2.) * ((D_s-D_l)/(D_s*D_l)))
	return theta_E


#this should be the angular seapartion between the lens and source, not sure if this is correct
def total_magnification(D_s, D_l, G, M, c):
	theta_E = einstein_radius(D_s, D_l, G, M, c)
	u = (D_s - D_l)/theta_E
	A = ((u**2.)+2.)/(u*np.sqrt((u**2.)+4))
	return a

def radians_to_arcseconds(radians):
	arcseconds = radians*206265.
	return arcseconds

#this function never turns around so im not sure that this is correct
def magnification_over_time(u_min, t_0, t_E):
	t= np.linspace(0, 200, 1000)
	u = np.zeros(len(t))
	for i in range(len(t)):
		u[i] = np.sqrt((u_min**2.) + ((t[i] - t_0)/t_E)**2.)
		#t[i]=t[i-1]+1
		#print t[i]
		i+=1
	return u, t

def projected_einstein_radius(G, M, c, D_s, D_l):
	r_tilde_e = np.sqrt(((4.*G*M)/c**2.)*((D_s*D_l)/(D_s-D_l)))
	return r_tilde_e

def einstein_parallax(G, M, c, D_s, D_l):
	r_tilde_e = projected_einstein_radius(G, M, c, D_s, D_l)
	return 1./r_tilde

def lens_source_angular_sep(theta, G, M, c, D_s, D_l):
	beta = theta - ((4.*G*M)/((c**2.)*theta))*((D_s - D_l)/(D_s*D_l))
	return beta


theta_E = einstein_radius(D_s, D_l, G, M, c)
#print radians_to_arcseconds(theta_E)

y, x = magnification_over_time(u_min, t_0, t_E)
y2, x2 = magnification_over_time(u_min, t_0, t_E)
a, b = magnification_over_time(1.1, 0., 19.)
a2, b2 = magnification_over_time(1.1, 0., 19.)
#plt.plot(x, y)
#plt.plot((x2*-1.)+(np.max(x)*2.), y2)
plt.plot(a, b)
plt.plot((a2*-1.)+(np.max(a)*2.), b2)
plt.show()