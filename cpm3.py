import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import sys, os
import re
import random

from calib import convert_to_pressure, import_HCL_data
from Table import Tableread, import_PSI
global RHO

K_1 = 1.3
A_PLUS = 26
KIN_VIS = 14.9e-6
RHO = 1.2
CHORD_LENGTH = 0.2

TIMEAXIS = 0
LIMIT_K = True

xc_pos = [0.21, 0.32, 0.41]
pos_table = Tableread("pos.txt", format = ['f', 'f'], separator = '\t', skiplines=2)

diameters_mm = [0.5,0.8,1]
diameters = [i / 1000 for i in diameters_mm]

d_eff = [1.3 * i / 2 for i in diameters]

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "data_UTF8"
HCL_DATA_RE = re.compile("HCL_(\\d+)_ms_(\\d+)_deg_Pos(\\d).dat")

measurements = []

for file_entry in os.scandir(DATA_DIR):

	if file_entry.is_file():
		#print(file_entry.name)
		match = HCL_DATA_RE.match(file_entry.name)
		if match is not None:
			groups = match.groups()
			assert len(groups) == 3

			measurements.append(groups)

Vel_set = list(set(([i[0] for i in measurements])))
AoA_set = list(set(([i[1] for i in measurements])))
Pos_set = list(set(([i[2] for i in measurements])))

Vel_set.sort()
AoA_set.sort()
Pos_set.sort()

# def stoch_grad_descent(f, initial_value):
#	 x = random.randint()

def van_dries(y_max, K, p):
	y = np.linspace(0,y_max,1000)
	pressure_corr_term = 1 + p * y
	exp_term = 1 - np.e**(- y * np.sqrt(pressure_corr_term) / A_PLUS)

	numerator = 2 * pressure_corr_term
	denominator = 1 + np.sqrt(1 + 4 * (K * y)**2 * pressure_corr_term * exp_term**2)
	
	u_plus = np.trapz(numerator/denominator, y)
	return u_plus


def cpm_estimation(K,q_target, p_grad, d):
	#print(q_target)
	q_plus = q_target * d**2 / ( 4 * RHO * KIN_VIS**2) 
	y_eff = K_1 * d/2

	tau_w = 1 # initial tau_guess
	u_tau = np.sqrt(tau_w/RHO)
	tau_old = -np.inf
	n = 0 
	while (abs(tau_w - tau_old) > 0.001) and (n < 1000):
		if n % 100 == 0:
			print(n)
		y_plus = u_tau * y_eff / KIN_VIS
		p = KIN_VIS/RHO / u_tau**3 * p_grad
		res = van_dries(y_plus, K, p)
		tau_plus = 2 * q_plus / res**2

		u_tau = np.sqrt(tau_plus) * 2 * KIN_VIS / d
		tau_old = tau_w
		tau_w = u_tau**2 * RHO
		n += 1
		#print(tau_w)

	return tau_w

def check(tab, k_i, i, p_grad):
	q = convert_to_pressure(i, tab['avg'][i])
	print("Q: ", q)
	return cpm_estimation(k_i, q, p_grad, diameters[i])


def cpm3_iter(tab,k, p_grad):
	taus = [check(tab, k, i, p_grad) for i in range(3)]
	print(taus)
	tau = np.average(taus)
	curr_error = np.sqrt(sum([abs(i - tau)**2 for i in taus]))
	return tau, k, curr_error

def calc_cpm3(tab, p_grad):
	global RHO
	tab['avg'] = np.mean(tab['data'], axis=TIMEAXIS)
	tab['std'] = np.std(tab['data'], axis=TIMEAXIS)
	RHO = tab['rho']

	k = 0.4
	n = 3
	
	tau, k, error = cpm3_iter(tab, k, p_grad)

	while n < 1000:
		new_tau_p, _, new_error_p = cpm3_iter(tab, k + 1/n, p_grad)
		new_tau_n, _, new_error_n = cpm3_iter(tab, k - 1/n, p_grad)

		if new_error_p < error:
			k = k + 1/n
			error = new_error_p
			tau = new_tau_p

		elif new_error_n < error:
			k = k - 1/n
			error = new_error_n
			tau = new_tau_n
			
		else:
			n *= 2

		n += 1

	if LIMIT_K:
		if k > 1:
			tau, k, error = cpm3_iter(tab, 1, p_grad)
		
		if k < 0:
			tau, k, error = cpm3_iter(tab, 0, p_grad)

	return (tau, k, error)

print(Vel_set)
print(Pos_set)
print(AoA_set)

print(measurements)

try:
	f = open("res.txt", "x")
except:
	f = open("res.txt", "w")

f.write("vel\tAoA\tpos\ttau\tk\terror\n\n")
for vel in Vel_set: 
	for pos in Pos_set:
		for AoA in AoA_set:
			if (vel, AoA, pos) not in measurements:
				# if a given combination doesn't exist, skip it
				continue
			
			p_filename = DATA_DIR + "/" + f"PSI_{vel}_ms_{AoA}_deg_Pos{pos}.dat"
			filename = DATA_DIR + "/" + f"HCL_{vel}_ms_{AoA}_deg_Pos{pos}.dat"

			
			print(filename)
			hcl_data = import_HCL_data(filename)
			p_data = import_PSI(p_filename)
			i = int(pos) + 8
			p_grad = p_data['p'][i+1] - p_data['p'][i-1] / (pos_table[i+1][0]-pos_table[i-1][0]) / CHORD_LENGTH
			tau, k, error = calc_cpm3(hcl_data, p_grad)
			f.write(f"{vel}\t{AoA}\t{pos}\t{tau}\t{k}\t{error}\n")

		

