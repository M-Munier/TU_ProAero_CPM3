import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import sys, os
import re

from calib import convert_to_pressure, import_HCL_data


DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "data_UTF8"
HCL_DATA_RE = re.compile("HCL_(\\d+)_ms_(\\d+)_deg_Pos(\\d)_new.dat")

measurements = []

for file_entry in os.scandir(DATA_DIR):

    if file_entry.is_file():
        print(file_entry.name)
        match = HCL_DATA_RE.match(file_entry.name)
        if match is not None:
            groups = match.groups()
            assert len(groups) == 3

            measurements.append(groups)

Vel_set = set(([i[0] for i in measurements]))
AoA_set = set(([i[1] for i in measurements]))
Pos_set = set(([i[2] for i in measurements]))

def calc_cpm3(tab):
    
    pass

for vel in Vel_set: 
    for pos in Pos_set:
        for AoA in AoA_set:
            filename = f"HCL_{vel}_{AoA}_deg_Pos{pos}_new.dat"
            hcl_data = import_HCL_data(filename)
            calc_cpm3(hcl_data)
