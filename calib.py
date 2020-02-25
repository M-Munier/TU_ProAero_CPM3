import os, sys
import json
import numpy as np
import matplotlib.pyplot as plt

NUM_SENSORS = 3
CALIB_DATA_PATH = "data_UTF8/calib"

CALIBRATION_RES_FILENAME = "calibration_result.json"
CALIBRATION_TO_P_FILENAME = "calibration_inverse.json"

voltage_to_p_curve = None
try:
    voltage_to_p_curve = json.load(open(CALIBRATION_TO_P_FILENAME))
except IOError:
    pass

def convert_to_pressure(sensor_idx: int, voltage_data : np.ndarray):
    if voltage_to_p_curve is None:
        raise FileNotFoundError("Calibration Data was not Found")
    if sensor_idx < NUM_SENSORS:
        coeffs = voltage_to_p_curve[f'sensor_{sensor_idx}']
        return np.polyval(coeffs, voltage_data)
    else:
        raise ValueError(f"Sensor IDX {sensor_idx} is out of bounds for a total of {NUM_SENSORS}")

def import_HCL_data(filename):
    """
Projekt Aerodynamik Übung CPM3
P_amb in [Pa]:  102415 T_amb in [K]: 297
Channel 1-3 HCL Druckdaten [V] 
2,744	2,781	2,759
"""

    def parseheader(hdrline):
        hdrline = hdrline.split("P_amb in [Pa]:  ")[1]
        hdrline = hdrline.split(" T_amb in [K]: ")

        p_amb = float(hdrline[0])
        T = float(hdrline[1])

        rho = p_amb/(T * 287.06)

        return {
            'p_amp' : p_amb,
            'T'     : T,
            'rho'   : rho,
        }

    with open(filename) as f:
        title = f.readline()
        header = f.readline()
        tabhead = f.readline()
        
        data = [[float(num) for num in line.replace(",",".").split("\t")] for line in f]

    table = parseheader(header)
    table['data'] = np.array(data)

    return table


if __name__ == "__main__":
    calib_data = {}
    delta_p_list = []

    calib_dir = os.scandir(CALIB_DATA_PATH)
    for i in calib_dir:
        delta_p = float(i.name[4:].split("_mbar")[0])
        delta_p_list.append(delta_p)

    delta_p_list.sort()

    for delta_p in delta_p_list:
        print(delta_p)
        print(int(delta_p))
        filename = "KAL_" + str(int(delta_p)).zfill(3) + "_mbar.dat"
        path = os.path.join(CALIB_DATA_PATH, filename)
        tab = import_HCL_data(path)
        print(f"Kalibrationspunkt {delta_p} Pa")
        print(tab)
        print()
        # plt.plot(tab['data']), axis=1
        # plt.show()

        TIMEAXIS = 0 
        tab['avg'] = np.mean(tab['data'], axis=TIMEAXIS)
        tab['std'] = np.std(tab['data'], axis=TIMEAXIS)

        print(f"Average: {tab['avg']}")
        print(f"Standard deviation: {tab['std']}")
        
        calib_data[delta_p] = tab


    calib_result = {}
    calib_inverse = {}

    for sensor_idx in range(NUM_SENSORS):

        averages =      [calib_data[dp]['avg'][sensor_idx] for dp in delta_p_list]
        deviations_3 =    [calib_data[dp]['std'][sensor_idx] * 3 for dp in delta_p_list]

        calib_coefficients, res = np.polyfit(delta_p_list, averages, deg=1, full=True)
        print(res)
        print(f"Sensor {sensor_idx}:")
        print(f"Coefficients: {calib_coefficients}")

        x = np.linspace(0,25000,2)
        plt.plot(x, np.polyval(calib_coefficients, x))
        plt.errorbar(delta_p_list, averages, yerr=deviations_3, fmt='none')

        calib_result[f"sensor_{sensor_idx}"]  = list(calib_coefficients)
        calib_inverse[f"sensor_{sensor_idx}"] = [1/calib_coefficients[0], -calib_coefficients[1]/calib_coefficients[0]]


    json.dump(calib_result, open(CALIBRATION_RES_FILENAME, "w"))
    json.dump(calib_inverse, open(CALIBRATION_TO_P_FILENAME, "w"))

    plt.show()