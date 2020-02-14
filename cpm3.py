import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import sys, os
import re

from calib import convert_to_pressure


DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "data_UTF8"
HCL_DATA_RE = re.compile("HCL_(\d+)_ms_(\d+)_deg_Pos(\d)_new.dat")

for file_entry in os.scandir(DATA_DIR):

    if file_entry.is_file():
        print(file_entry.name)
        match = HCL_DATA_RE.match(file_entry.name)
        if match is not None:
            groups = match.groups()
            assert len(groups) == 3

            vel, AoA, pos = *groups
            

        
