import math
import numpy as np
from statistics import mean 
import matplotlib.pyplot as plt
import pandas as pd

# Constants (assumed for air at standard conditions)
RHO_AIR = 1.225  # density of air (kg/m^3)
MU_AIR = 1.81e-5  # dynamic viscosity of air (Pa.s)
KINEMATIC_VISCOSITY = MU_AIR / RHO_AIR  # kinematic viscosity (m^2/s)

# Function to calculate Reynolds number
def reynolds_number(velocity, hydraulic_diameter):
    return velocity * hydraulic_diameter / KINEMATIC_VISCOSITY

def hydraulic_diameter(width, height):
    return 2 * (width * height) / (width + height)

# Colebrook-White equation for friction factor in turbulent flow
def colebrook_white(d, Re, e_d):
    return 0.25 / (math.log10((e_d / (3.7 * d)) + (5.74 / Re**0.9)))**2

# Darcy-Weisbach equation for head loss
def darcy_weisbach(friction_factor, length, velocity, hydraulic_diameter):
    return friction_factor * (length / hydraulic_diameter) * (velocity**2 / (2 * 9.81))

def minor_loss(k,v):
    return k * (v ** 2) / ( 2 * 9.81 )

def oriface_k(d0,d2,beta): # Computes the loss through a multi hole sharp oriface with transiton


    l = 1 + 0.622 *(1 - 0.215 * beta ** 2 - 0.785 * beta ** 5) # Velocity ratio

    k = 0.0696 * (1 - beta ** 5) * (l ** 2) + (l - (d0/d2)**2 ) ** 2 # Compute k factor 

    return k

def Calculate_Loss(airgap,Q,car_effects=False):
    '''
    |------------------------------------------------------------------|
    |~~~~~~~~~~~~~~~~~~~~~~~~~~PACK Variables~~~~~~~~~~~~~~~~~~~~~~~~~~|
    |------------------------------------------------------------------|
    '''
    #Input everything in SI for the love of christ

    N_Slots =  8  # Number of airflow slots 
    N_Hole_Per_Slot_In = 8
    N_Hole_Per_Slot_Out = 6
    D_Hole_In = 0.0098
    D_Hole_Out = 0.01
    t_hole = 0.0043
    L_Cell = 0.12
    W_Cell = 0.0965
    Airgap = airgap
    A_wall_in = 0.0313 # Area of the Inlet Wall
    A_wall_out = 0.000217 * N_Slots # Area of the Outlet Wall of Channel; think of this as the wall space in the bottom chanel
    A_naca_in = 0.0048 # Naca in


    '''
    |------------------------------------------------------------------|
    |~~~~~~~~~~~~~~~~~~~~~~~~~Operating Points~~~~~~~~~~~~~~~~~~~~~~~~~|
    |------------------------------------------------------------------|
    '''


    Q_Sweep = np.linspace(0,4/60,500)

    '''
    |------------------------------------------------------------------|
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~Loss Factors~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    |------------------------------------------------------------------|
    '''
    k_in = 0.57
    k_90 = 1.5
    k_contraction = 0.4 # Estimate from literature table

    d1 = hydraulic_diameter(Airgap,L_Cell) # Hydralic Diameter of Channel
    beta1 = np.sqrt ((((np.pi * (D_Hole_In) / 2) ** 2 ) * N_Hole_Per_Slot_In * N_Slots) / A_wall_in) # Porosity Of Oriface 
    d2 = hydraulic_diameter(Airgap,W_Cell) # Hydralic Diameter of Channel
    beta2 = np.sqrt (A_wall_out / ((np.pi * (D_Hole_Out / 2) ** 2 ) * N_Hole_Per_Slot_Out * N_Slots)) # Porosity Of Oriface 
    d1_5 = mean([d1,d2])

    k_hole_in = oriface_k(D_Hole_In,d1,beta1)
    k_hole_out = oriface_k(d2,D_Hole_Out,beta2)
    k_out = 1
    k_mesh = 0.254 # Online German Calculator http://www.pressure-drop.online
    k_tune = 12
    k_tune_car = 15


    '''
    |------------------------------------------------------------------|
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~Velocities~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    |------------------------------------------------------------------|
    '''

    Q_Slot = Q / N_Slots
    v_0 = Q_Slot / (N_Hole_Per_Slot_In * (np.pi * (D_Hole_In / 2) ** 2))
    v_1 = Q_Slot / ( Airgap * L_Cell)
    v_2 = Q_Slot / ( Airgap * W_Cell)
    v_1_5 = mean([v_1,v_2])
    v_3 = Q_Slot / (N_Hole_Per_Slot_Out * (np.pi * (D_Hole_Out / 2) ** 2))
    v_naca = Q * 3 / A_naca_in

    '''
    |------------------------------------------------------------------|
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~Minor Losses~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    |------------------------------------------------------------------|
    '''
    hin = minor_loss(k_in,v_0)
    h0 = minor_loss(k_hole_in,v_0)
    h1 = minor_loss(k_90,v_2) # Average of velocities 
    h2 = minor_loss(k_hole_out,v_3)
    hout = minor_loss(k_out,v_3)
    hcar = (Q ** 2) * 39026.733 # Inertia Loss Accociated with the intake of the car to the pack from CFD using System Curve Formula
    hfanshroud = (Q ** 2) * 22619.56461
    hmesh_in = minor_loss(k_mesh,v_naca)
    hmesh_out = minor_loss(k_mesh,v_3)
    hcontraction = minor_loss(k_contraction,v_2)


    '''
    |------------------------------------------------------------------|
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~Major Losses~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    |------------------------------------------------------------------|
    '''

    Re_Channel = reynolds_number(v_1_5,d1_5)
    e_cell = 0.005 * 10 ** (-3)

    ed = e_cell / d1_5

    f = colebrook_white(d1_5,Re_Channel,ed)

    H1 = darcy_weisbach(f,L_Cell,v_1_5,d1_5)



    '''
    |------------------------------------------------------------------|
    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~Sum Losses~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
    |------------------------------------------------------------------|
    '''
    if car_effects == True:
        H_Total = (sum([H1,hin,h0,h1,h2,hout,hcar,hmesh_in,hfanshroud])) * k_tune_car
    elif car_effects == False:
        H_Total = (sum([H1,hin,h0,h1,h2,hout,hfanshroud,hmesh_out,hcontraction])) * k_tune

    print(H_Total) 
    
    return(H_Total)


'''
|------------------------------------------------------------------|
|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Sweep~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|------------------------------------------------------------------|
'''


plt.figure(figsize=(10, 6))

dp_sweep = []

Q_Sweep = np.linspace(0,0.05,100)

aigap_range = np.linspace(0.003,0.0055,17)
'''
for g in aigap_range:
    dp_res = []
    for i in Q_Sweep:
        dp_res.append(Calculate_Loss(g,i,car_effects=False))
    f = round(g * 1000,2)
    plt.plot(Q_Sweep,dp_res, label=f'{f}mm Airgap')
'''
#'''
g = 0.004

dp_res = []
for i in Q_Sweep:
    dp_res.append(Calculate_Loss(g,i,car_effects=False))
f = round(g * 1000,2)
dp1 = dp_res[0]
plt.plot(Q_Sweep,dp_res, label=f'{f}mm Airgap')

dp_res = []

g = 0.0055
Q_Sweep = np.linspace(0,0.05,100)

for i in Q_Sweep:
    dp_res.append(Calculate_Loss(g,i,car_effects=False))
f = round(g * 1000,2)
dp1 = dp_res[0]
plt.plot(Q_Sweep,dp_res, label=f'{f}mm Airgap')
'''
for i in Q_Sweep:
    dp_res.append(Calculate_Loss(g,i,car_effects=True))
f = round(g * 1000,2)
plt.plot(Q_Sweep,dp_res, label=f'{f}mm Airgap Car effects')
dp2 = dp_res[0]

print(dp1/dp2)
'''
#'''
'''
|------------------------------------------------------------------|
|~~~~~~~~~~~~~~~~~~~~~~~~~Fan Curve~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|------------------------------------------------------------------|
'''

Q_fan = [0.040852, 0.449512, 1.048852, 1.430268, 2.056852, 2.601732, 2.983148, 3.092096, 3.419024]
Q_fan_m3_s = [element / 60 for element in Q_fan]
dP_fan =[359.71, 319.55, 267.33, 227.17, 178.95, 150.79, 108.63,90.57, 46.41]


MFE_24_Fan = pd.read_csv("MFE24_Fan.csv")
Q_fan = MFE_24_Fan.values[:,0] / 60
dP_fan = MFE_24_Fan.values[:,1]

J01 = pd.read_csv("J04_Fan.csv")
Q_fan_J01 = J01.values[:,0] / 60
dP_fan_J01 = J01.values[:,1]

Test_Fan = pd.read_csv("Beefy_Fan.csv")
Q_fan_Test = Test_Fan.values[:,0]
dP_fan_Test = Test_Fan.values[:,1]


 # Plot the first dataset
plt.plot(Q_fan, dP_fan, label="MFE24 Fan",color='black', linewidth=2, linestyle='--')
plt.plot(Q_fan_J01, dP_fan_J01, label="J04 Fan",color='red', linewidth=2, linestyle='--')
plt.plot(Q_fan_Test, dP_fan_Test, label="Test Fan",color='Purple', linewidth=2, linestyle='--')


'''
|------------------------------------------------------------------|
|~~~~~~~~~~~~~~~~~~~~~~~~~Fan Curve~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
|------------------------------------------------------------------|
'''

# Add labels and title
plt.xlabel(r"Volumetric Flow Rate (m$^3$/s)")
plt.ylabel("Pressure Drop (Pa)")

# Add a legend
plt.legend()
plt.title("Module System Curve With Respect to Airgap")
# Show the plot
plt.grid(True)
plt.savefig(fname="Fan_Sweep.png",dpi=1000)
plt.show()


