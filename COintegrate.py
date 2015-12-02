import numpy as np
import pylab as plt
import physical as ph
import reaction as re
import parameters as par
import abundances as ab
import os


# Create the run directory:
if not os.path.exists('runs/%s'%par.directory): os.makedirs('runs/%s'%par.directory)

#----------- Initial fluid parameters ------------
T_cs = par.T_cs_init
c_s = np.sqrt(ph.kB * par.T_cs_init/(ph.mp))
R_cs = par.R_cs_init
delta = par.delta_init
n = par.n_init

#------------ Initial reaction rates ------------
# In the future K_nonthermal will depend on time, but 
# for now it is constant.
K_nontherm = par.epsilon*par.L_gamma*ph.mp/par.M_cs/par.W_d
K_ra = re.formation(T_cs)
K_th = re.thermal(T_cs)
K_nth = K_nontherm
# K_rd is the sum of both destruction rates.
K_rd = K_th + K_nth

#------------ Initial abundances ------------
Y_C_free = par.Y_C_tot
Y_O_free = par.Y_O_tot
Y_CO = ab.Y_CO_equil(par.Y_C_tot, par.Y_O_tot, n, K_ra, K_rd)
Y_C_free = par.Y_C_tot - Y_CO
Y_O_free = par.Y_O_tot - Y_CO
integrate_flag = 0 # flag = 1 when integrating
adaptive_dt_flag = 0 # flag = 1 when dt is smaller than default

# Initial time
t = par.tmin

#------------- OUTPUT FILES -----------------------#
fractionfile = open("runs/%s/fractions.txt"%par.directory, 'w')
ratesfile = open("runs/%s/rates.txt"%par.directory, 'w')
thermofile = open("runs/%s/thermo.txt"%par.directory, 'w')
#Then open for appending
fractionfile = open("runs/%s/fractions.txt"%par.directory, 'a')
ratesfile = open("runs/%s/rates.txt"%par.directory, 'a')
thermofile = open("runs/%s/thermo.txt"%par.directory, 'a')
# headers: 
fractionfile.write("# t     Y_CO      Y_C_free     Y_O_free     flag\n")
ratesfile.write("# K_ra      K_therm     K_nontherm \n")
thermofile.write("# T_cs     n      delta      R_cs       c_s  \n")
# first line: 
fractionfile.write("%.5f %.5f %.5f %.5f %i %i \n"%(t, Y_CO, Y_C_free, Y_O_free, integrate_flag, adaptive_dt_flag))
ratesfile.write("%.5f %.5e %.5e \n"%(K_ra, K_th, K_nth))
thermofile.write("%.5f %.5f %.5f %.5f %.5f \n"%(T_cs, n, delta, R_cs, c_s))


#----------- TIME EVOLUTION ----------------------# 
while t < par.tmax:
	# What to do when Y_C_free goes negative?? For now, break
	#if Y_C_free < 0.: break

	#---------------- First choose dt ------------------#
	dt = par.dt_init
	t = t + dt
	Y_CO_equilibrium = ab.Y_CO_equil(par.Y_C_tot, par.Y_O_tot, n, K_ra, K_rd)
	#if t > par.t_integrate:
	if Y_CO_equilibrium > 1e-5:
		integrate_flag = 1
		adaptive_dt_flag = 0
		delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		# Reduce dt if delta_YCO is too large
		while (np.absolute(delta_YCO) > Y_C_free \
			or np.absolute(delta_YCO) > Y_O_free):
			adaptive_dt_flag = 1
			dt = dt/2.
			delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		# Also reduce dt so we don't form more than some percent of CO
		#while ((Y_CO > par.CO_amount  or Y_CO + delta_YCO > par.CO_amount) \
		if (Y_CO > par.CO_amount  or Y_CO + delta_YCO > par.CO_amount):
			while (np.absolute(delta_YCO) > par.percentCO*Y_CO):
				adaptive_dt_flag = 1
				dt = dt/2.
				delta_YCO = ab.dYCO_dt(Y_CO, Y_C_free, Y_O_free, n, K_ra, K_rd) * dt
		# Evolve Y_CO
		Y_CO = np.maximum(delta_YCO + Y_CO,0.)
	else:
		# Otherwise solve for the equilibrium density 
		integrate_flag = 0
		Y_CO = Y_CO_equilibrium
	#--------------------------------------------------#

	# Expand the shell
	R_cs = par.R_cs_init + par.v_ej*t
	delta = par.delta_init + c_s*t

	if t < (R_cs**2 * delta/par.v_ej**3)**(1./3.): 
		n = par.M_cs/(4.*np.pi*ph.mp*R_cs**2 * delta)
		t_exp = t/2.
	else:
		n = par.M_cs/(4.*np.pi*ph.mp*(par.v_ej*t)**3)
		t_exp = t/3.

	# Temperature and sound speed decrease
	T_cs = par.T_cs_init * (n / par.n_init)**(par.gamma-1.)
	c_s = np.sqrt(ph.kB*T_cs/(ph.mp))
	#---------------------------------------------#

	# evolve the rates
	K_ra = re.formation(T_cs)
	K_th = re.thermal(T_cs)
	K_nth = K_nontherm
	K_rd = K_th + K_nth

	n_CO = Y_CO * n
	Y_C_free = par.Y_C_tot - Y_CO
	Y_O_free = par.Y_O_tot - Y_CO
	n_C_free = Y_C_free * n
	n_O_free = Y_O_free * n

	# Append to text file
	fractionfile.write("%.5f %.5f %.5f %.5f %i %i \n"%(t, Y_CO, Y_C_free, Y_O_free, integrate_flag, adaptive_dt_flag))
	ratesfile.write("%.5e %.5e %.5e \n"%(K_ra, K_th, K_nth))
	thermofile.write("%.5f %.5f %.5f %.5f %.5f \n"%(T_cs, n, delta, R_cs, c_s))

	#if Y_C_free[i] + Y_O_free[i] + 2.*Y_CO[i] != 1.: 
	#	print 'ERROR: Sum of number fractions is ne 1!'
	#	break
	print t
#---------------------------------------------------------#

fractionfile.close()
ratesfile.close()
thermofile.close()


