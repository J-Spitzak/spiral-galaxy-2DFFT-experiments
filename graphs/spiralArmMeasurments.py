import csv
import datetime
import math
import sys
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.optimize as optim

rownum = 0

manualPA = []
predictionPA = []
SFRavg = []
name = []

surveyfile = open("spiralArmMeasurments.csv",'r')	#open csv file with UG survey data

for row in csv.reader(surveyfile): 		#Read pairs file row by row
    if (rownum > 2):
        if len(row[13].strip()) != 0 and len(row[6].strip()) != 0:
            name.append(str(row[1]))
            manualPA.append(int(row[6]))
            predictionPA.append(int(row[10]))
            SFRavg.append(float(row[13]))
    rownum = rownum + 1


plt.scatter(manualPA,SFRavg,color='blue')
plt.scatter(predictionPA,SFRavg,color='red')

plt.ylabel("SFR")
plt.xlabel("Pitch angle")

plt.show()
