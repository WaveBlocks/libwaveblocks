import csv

from WaveBlocksND import *

def WriteHaWpCoefficientsAsCSV(coefficients, shape_enum, filename):
    with open(filename, "w") as filehandle:
        csvfile = csv.writer(filehandle, delimiter=' ')
        for mindex, coeff in zip(shape_enum, coefficients):
            csvfile.writerow(mindex + tuple(coeff))