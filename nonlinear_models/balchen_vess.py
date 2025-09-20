# ===============================
# VESSEL CONSTANTS
# ===============================
# Wind coefficients
C_X = 0.6
C_Y = 0.8
C_N = 0.1

# Windage areas [m^2]
A_F = 500 # frontal windage area [m^2]
A_L = 1100 # lateral windage area [m^2]

# Vessel dimensions
LENGTH = 73.2  # overall length [m]

# Mass and inertia
M_SURGE = 4e6   # [kg]
M_SWAY = 4e7    # [kg]
I_YAW = 4.7e10  # [kg.m^2]

# Damping coefficients  # taken from page 151 and 152 in Balchen paper

D_SURGE = 5e-5
D_SWAY = 21e-5
D_YAW1 = 1.1e-10
D_YAW2 = 201e-15

# Air density
RHO_AIR = 1.23  # [kg/m^3]


import math
import numpy as np

def smooth(a, WSZ):
    """
    SMOOTH smoothing func

    Input:
        a - ...
        wSZ - ...
    Output:
       
    """
    WSZ = int(WSZ)
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))
    
def slowly(N):
    """
    SLOWLY slowly warying 

    Input:
        N - Length samples
    Output:
        A -...
    """
    A = np.zeros([N,1])
    a = np.random.rand(4,1)
    for i in range(0, N):
        if i <= N/4:
            A[i,0] = a[0,0]
        elif i > N/4 and i <= N/2:
            A[i,0] = a[1,0]
        elif i > N/2 and i <= 3*N/4:
            A[i,0] = a[2,0]
        else:
            A[i,0] = a[3,0]
    return A   
    
def wind_force(V_w, betta, psi, v_su, v_sw):
    """
    WIND_FORCES Calculate the wind forces acting on a vessel.

    Input:
        V_w - Wind speed [m/s] in the NED frame.
        betta - Wind direction [deg] in the NED frame.
        psi - Vessel heading [rad].
        v_su - Wind relative speed in surge.
        v_sw - Wind relative speed in sway.

    Output:
        F_w - wind forces in surge, sway, and yaw moment.
    """

    betta = betta * math.pi/180 # convert to radian
    u_w = V_w * math.cos(betta - psi) # wind speed in surge
    v_w = V_w * math.sin(betta - psi)# wind speed in sway
    u_rw = v_su - u_w # wind relative speed in surge
    v_rw = v_sw - v_w # wind relative speed in sway
    V_r = math.sqrt(u_rw**2 + v_rw**2) # wind relative speed
    F_w_su = 0.5*RHO_AIR*C_X*A_F*math.cos(betta)* V_r**2
    F_w_sw = 0.5*RHO_AIR*C_Y*A_L*math.sin(betta)* V_r**2
    N_w = 0.5*RHO_AIR*C_N*A_L*LENGTH*math.sin(2*betta)* V_r**2
    F_w = np.zeros([3,1])
    F_w[0] = F_w_su
    F_w[1] = F_w_su
    F_w[2] = N_w
    return F_w
 
def current_force(V_c, x):
    """
    CURRENT_FORCE Calculate the wind forces acting on a vessel.

    Input:
        V_w - Wind speed [m/s] in the NED frame.
        betta - Wind direction [deg] in the NED frame.
        psi - Vessel heading [rad].
        v_su - Wind relative speed in surge.
        v_sw - Wind relative speed in sway.

    Output:
        F_w - wind forces in surge, sway, and yaw moment.
    """
    R = np.matrix([
    [math.cos(x[2]), -math.sin(x[2]), 0],
    [math.sin(x[2]),  math.cos(x[2]) ,0],
    [0, 0 ,1]])
    F_c = R.T*V_c.T #tranform from NED to body coordinate system
    return F_c
    
def nonlinear(x, u, F_c, F_w):
    """
    NONLINEAR "NONLINEAR" state space model of a vessel

    Input:
        x - states body
        u - input
        F_c - Current forces
        F_W - Wind Forces 

    Output:
       f - dot x
    """
	f = np.zeros([6,1])
	f[0] = x[3]
	f[1] = x[4]
	f[2] = x[5]
	f[3] = -D_SURGE/M_SURGE*np.abs(x[3]-F_c[0])@(x[3]-F_c[0]) + 1/M_SURGE*(u[0]+F_w[0])
	f[4] = -D_SWAY/M_SWAY*np.abs(x[4]-F_c[1])@(x[4]-F_c[1]) + 1/M_SWAY*(u[1]+F_w[1])
	f[5] = -D_YAW1/I_YAW*np.abs(x[5])@x[5] - D_YAW2/I_YAW*np.abs(x[4]-F_c[1])@(x[4]-F_c[1]) + 1/I_YAW*(u[2]+F_w[2]+F_c[2])
	return f    
 
