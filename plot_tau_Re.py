import sys

import numpy as np
import matplotlib.pyplot as plt

from Table import Tableread, import_PSI
from calib import import_HCL_data

RESULT_FILE = sys.argv[1] if len(sys.argv) > 1 else "res.txt"
DATA_DIR = sys.argv[2] if len(sys.argv) > 2 else "data_UTF8"

KIN_VIS = 14.9e-6
CHORD_LENGTH = 0.2

ADD = "LIMITED"

res = Tableread(RESULT_FILE, format=['i'] * 3 + ['f'] * 3, separator="\t", skiplines=2)
xc_pos = [0.21, 0.32, 0.41]

plot_data = {}

for i in [0,6,10]:
    for j in [10,15]:
        plot_data[(i,j)] = [[],[],[],[],[]]


for line in res:
    vel, AoA, pos, tau, k, error = tuple(line)
    p_data = import_PSI(DATA_DIR + "/" + f"PSI_{vel}_ms_{AoA}_deg_Pos{pos}.dat")
    hcl_data = import_HCL_data(DATA_DIR + "/" + f"HCL_{vel}_ms_{AoA}_deg_Pos{pos}.dat")
    
    q = p_data['p'][-1]
    u = np.sqrt(q / hcl_data['rho'])
    
    Re = u * xc_pos[pos-1] / KIN_VIS / CHORD_LENGTH
    
    plot_data[(AoA,vel)][0].append(Re)
    plot_data[(AoA,vel)][1].append(tau)
    plot_data[(AoA,vel)][2].append(tau/q)
    plot_data[(AoA,vel)][3].append(k)
    plot_data[(AoA,vel)][4].append(error)

colorscheme = [
    ([0.5, 0.5, 0.5], "o-"),
    ([0.4, 0.4, 0.4], "--x"),
    ([0.2, 0.2, 0.2], "-.^"),
    ([0, 0, 0],       ":h")
]

global style_idx
global label 
style_idx = 0

def cplot(x,y):
    global style_idx, label
    if label != None:
        plt.semilogx(x,y,colorscheme[style_idx][1], color=colorscheme[style_idx][0], label=label)
    plt.semilogx(x,y,colorscheme[style_idx][1], color=colorscheme[style_idx][0])
    style_idx += 1


for vel in [10,15]:
    for AoA in [6,10]:
        label = f"${vel} m/s$ mit {AoA}° Anstellwinkel"
        cplot(plot_data[(AoA,vel)][0], plot_data[(AoA,vel)][1])

plt.ylabel("Wandschubspannung $\\tau_w$ [$N/m^2$]")
plt.xlabel("Reynoldszahl [Re]")
plt.legend()
plt.savefig(f"images/tau_over_re{ADD}.eps", format="eps")
plt.clf()
plt.close()

style_idx = 0
for vel in [10,15]:
    for AoA in [6,10]:
        label = f"${vel} m/s$ mit {AoA}° Anstellwinkel"
        cplot(plot_data[(AoA,vel)][0], plot_data[(AoA,vel)][2])

plt.ylabel("Reibungskoeffizient $C_f$ [-]")
plt.xlabel("Reynoldszahl $Re$ [-]")
plt.legend()
plt.savefig(f"images/cf_over_re{ADD}.eps", format="eps")
plt.clf()
plt.close()

style_idx = 0
for vel in [10,15]:
    for AoA in [6,10]:
        label = f"${vel} m/s$ mit {AoA}° Anstellwinkel"
        cplot(plot_data[(AoA,vel)][0], plot_data[(AoA,vel)][-1])

plt.xlabel("Reynoldszahl [Re]")
plt.xlabel("Residual of the CPM3 Method [-]")
plt.legend()
plt.savefig(f"images/error_over_re{ADD}.eps", format="eps")
plt.clf()
plt.close()


    
