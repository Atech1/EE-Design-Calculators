
"""
Created by: Alec Ross
Date: 9/02/2025

This script is intended to be able to generate a set of graphs showing Resistance vs. Temperature for a Thermistor, given the Stein-Hart Equation parameters of the Thermistor

because some Thermistor manufacturers do not give all 3 A, B, C coiefficents,
I also implemented a way to take 3 points along the Resistance vs. Temperature to calculate the coefficents
It can be handy to estimate about how accurate the provided values in a datasheet are in being able to compare a calculated coefficent and the manufacturer provided ones. 

for the table of 3 points along a thermistors curve, at least 2 should be provided by the manufacturer. Thermistors are rated as X Resistance @ 25C, and Y Resistance at @ maximum rated temperature. 

To try to not be too confusing, the variable names in math functions are trying to stick to the common scienfic names derived in the classic formulas. 

the Matplotlib library is imported for handling the plotting
unum is a library to handle unit conversions, to try to make that part of the code simpler. 

                             THERMISTOR                                      
  +-----------------------------------------------------------------------+
R |                                                                       |
E |                                                                       |
S |                                                                       |
I +------+                                                                |
S |      |                                                                |
T |      |                                                                |
A |      |                                                                |
N |      |                                                                |
C |      |                                                                |
E |      +---------+                                                      |
  |                |                                                      |
  |                |                                                      |
  |                |                                                      |
  |                |                                                      |
  |                |                                                      |
  |                +----------+                                           |
  |                           |                                           |
  |                           |                                           |
  |                           |                                           |
  |                           +-----------------+                         |
  |                                             |                         |
  |                                             |                         |
  |                                             +-------------------+     |
  |                                                                 +---- |
  |                                                                       |
  +-----------------------------------------------------------------------+
                                                            TEMPERATURE            

"""


import math
import unum
import unum.units
from unum import Unum
import matplotlib.pyplot as plt

KOHM = Unum.unit('KOHM', 1000* unum.units.OHM)

# Thermistor Beta Equation Parameters

T_0 = 298.15 # Temperature 0K, for Kelvin Conversion
R25 = 10_000
β25 = 3435 # B25 and B85 are parameters common in the Beta thermistor function, some manufacturer's prefer this instead.
β85 = 3984
Temp_Range = [-55, 150] # basic range of temperatures handled for graphing, -55 - 150 C handles almost all deployable enivornments. 

# SteinHart Thermistor Equation Parameters
RTable = [[233.15, 334274], [298.15, 9_900], [398.15, 336.7]]    # for a different thermistor, a new list of list(2) should be made with the temperature in K, and resistance in Ohms of at least 3 points. 
NTCALUG01T103FATable = [[-30 + 273.15, 176132.54], [298.15, 10_000], [85 + 273.15, 1066.11]] 
A = 0.001141940
B = 0.000232183
C = 0.000000094

def Calculate_SteinHart_Hart_Coeiff(Points):
    """
    This is a crude implmentation derived directly from https://en.wikipedia.org/wiki/Steinhart%E2%80%93Hart_equation
    """
    L = [math.log(Points[0][1]), math.log(Points[1][1]), math.log(Points[2][1])]
    Y = [1/Points[0][0], 1/Points[1][0], 1/Points[2][0]]
    gamma = [(Y[1] - Y[0]) / (L[1] - L[0]), (Y[2] - Y[0]) / (L[2] - L[0])]
    c = ((gamma[1] - gamma[0]) / (L[2] - L[1])) * math.pow((L[0] + L[1] + L[2]), -1)
    b = gamma[0] - c * ((L[0]**2 + L[0] * L[1] + L[2]**2))
    a = Y[0] - ((b + L[0]**2 * c) * L[0])
    return (a, b, c)  

def Resistance_From_Temp_Beta(T, βerror = 0.00, Rerror = 0.00):
    if math.hypot(T - 25) <= math.hypot(T - 85) or β85 is None:
        return round((R25 + (Rerror * R25)) * math.exp((β25 + (β25 * βerror)) * ((1/(T + 273.15) - (1/T_0)))), 0) / 1000
    else:
        return round(R25 * math.exp(β85 * ((1/(T + 273.15) - (1/T_0)))), 0) / 1000

def Temp_From_Resistance_Beta(R):
     return (((1 / T_0 + (1 / β25) * math.log(R / R25)) ** -1) - 273.15) * unum.units.C

def Resistance_From_Temp_SH(T):
    Tk = T + 273.15
    x  = 1 / (2 * C) * (A - 1 / Tk)
    y  = math.sqrt(math.pow(B / (3 * C), 3) + pow(x, 2))
    R  = math.exp(pow(y - x, 1.0/3) - pow(y + x, 1.0/3))
    return round(R, 0) / 1000


def Temp_From_Resistance_SH(R):
    lnR = math.log(R)
    invT = A + B * lnR + C * math.pow(lnR, 3)
    return 1 / invT - 273.15

def Create_Beta_Plot(step = 5):
    BRValues = [Resistance_From_Temp_Beta(X) for X in range(Temp_Range[0], Temp_Range[1], step)]
    fig, ax = plt.subplots()
    ax.plot( [X for X in range(Temp_Range[0], Temp_Range[1], step)], BRValues, "o")
    ax.set_title("Thermistor Beta Curve")
    ax.set_ylabel("Resistance in KΩ")
    ax.set_xlabel("Temperature in C")
    plt.show()

def Create_SH_Plot(step = 5):
    #A, B, C = Calculate_SteinHart_Hart_Coeiff(DatasheetTable)
    RValues = [Resistance_From_Temp_SH(X) for X in range(Temp_Range[0], Temp_Range[1], step)]
    fig, ax = plt.subplots()
    ax.plot( [X for X in range(Temp_Range[0], Temp_Range[1], step)], RValues, "o")
    ax.set_title("Thermistor Steinhart - Hart Curve")
    ax.set_ylabel("Resistance in KΩ")
    ax.set_xlabel("Temperature in C")
    ax.text(50, 500, f"a = {A},\nb = {B},\nc= {C}")
    plt.show()

def Compare_SH_Beta_Plot(step = 5):
    SHRValues = [Resistance_From_Temp_SH(X) for X in range(Temp_Range[0], Temp_Range[1], step)]
    BRValues = [Resistance_From_Temp_Beta(X) for X in range(Temp_Range[0], Temp_Range[1], step)]
    fig, ax = plt.subplots()
    ax.plot( [X for X in range(Temp_Range[0], Temp_Range[1], step)], SHRValues, "o", label="Steinhart Eq")
    ax.plot([X for X in range(Temp_Range[0], Temp_Range[1], step)], BRValues, label="Beta Eq")
    ax.set_title("Steinhart vs Beta Curve")
    ax.set_ylabel("Resistance in KΩ")
    ax.set_xlabel("Temperature in C")
    ax.legend()
    plt.show()


print(Resistance_From_Temp_Beta(25))
print(Temp_From_Resistance_Beta(10000))
print(Resistance_From_Temp_SH(25))
print(Calculate_SteinHart_Hart_Coeiff(NTCALUG01T103FATable) )
print(Calculate_SteinHart_Hart_Coeiff(RTable) )

A, B, C = Calculate_SteinHart_Hart_Coeiff(NTCALUG01T103FATable)
print(Resistance_From_Temp_SH(25))
print(Temp_From_Resistance_SH(10_000))

A, B, C = Calculate_SteinHart_Hart_Coeiff(RTable)
print(Resistance_From_Temp_SH(25))
print(Temp_From_Resistance_SH(10_000))
Create_SH_Plot(5)
#Create_Beta_Plot(5)
#Compare_SH_Beta_Plot()