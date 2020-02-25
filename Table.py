import numpy
from binascii import hexlify
def Tableread(filename, format, separator = ' ', emptyvalues = False, skiplines=0):
    def convertfloat(s):
        return float(s.replace(",","."))
    validformats = {
    'f': convertfloat,
    'i' : int,
    's' : str,
    }
    table = []
    with open(filename) as f:
        for i in range(skiplines):
            f.readline()
        for line in f:
            row = []
            entries = line.split(separator)
            if not emptyvalues:
                entries = list(filter(None, entries))
            for i in range(len(format)):
                #print(entries)
                #print(i)
                row.append(validformats[format[i]](entries[i]))
            table.append(row)
    return table


def convertfloat(s):
    return float(s.replace(",","."))

def importdat(filename):
    d = {}
    with open(filename) as f:
        header = f.readline()
        #13:07:27	Temperatur:	21,25	Luftdruck: 	99969,42	Dynamischer Druck: 239,92
        d['T'] = float(header.split("Temperatur:")[1][1:].split("\t")[0].replace(",",".").strip())
        d['p_amb'] = float(header.split("Luftdruck: ")[1][1:].split("\t")[0].replace(",",".").strip())
        d['d'] = float(header.split("Dynamischer Druck:")[1][1:].replace(",",".").strip())
        d['table_fields'] = f.readline().strip().split("\t")
        #d[''] = float(lines[1].split("T_amb in [K]:")[1].strip())

        lines = f.readlines()
        d['values'] = numpy.array([list(filter(None, map(convertfloat, line.split("\t")))) for line in lines if line]).transpose()
    return d

def import_PSI(filename):
    d = {}

    with open(filename) as f:
        lines = f.readlines()
        d['p_amb'] = float(lines[1].split("P_amb in [Pa]:")[1].split(" T_amb")[0].strip())
        d['T'] = float(lines[1].split("T_amb in [K]:")[1].strip())
        d['p'] = list(filter(None, map(convertfloat, lines[-1].split("\t"))))

    return d