import numpy as np
import matplotlib.pyplot as plt
import math

# ---------------- Initializing variables ----------------
T_0 = 298.0          # K, inlet temperature
G_0 = 1.0e5          # K/m, inlet gradient dT/dx
L = 0.0015           # m
ATM = 101325.0       # Pa
N = 10000            # number of steps
dx = L / N

phi = 0.46           # equivalence ratio
p_in = 1.0 * ATM     # inlet pressure

# ---------------- Table 1 parameters (Case 1) ----------------
R = 8.3145           # universal gas constant
A = 1.4e8            # pre-exponential factor
NU = 1.6             # reaction order
EA = 121417.2        # activation energy
W = 0.02897          # molecular weight
LAM = 0.026          # thermal conductivity
CP = 1000.0          # heat capacity
QF = 5.0e7           # fuel calorific value

# ---------------- Inlet state ----------------
rho_0 = p_in * W / (R * T_0)       # inlet density (ideal gas)
YF_0 = phi / (4.0 + phi)           # inlet fuel mass fraction
T_ad = T_0 + QF * YF_0 / CP        # adiabatic flame temperature (all fuel burned)

# ---------------- Bisection bracket for the flame speed ----------------
u_L = 0.0            # m/s, lowest guess
u_R = 1.0            # m/s, highest guess
epsilon = 1e-16
converged = False

T_all = np.zeros(N + 1)            # stores T at every x (for plotting later)

# ================= Outer loop: bisection on u_0 (= s_L) =================
while (u_R - u_L) > epsilon and not converged:

    u_0 = (u_L + u_R) / 2          # current guess for the flame speed

    # restart the march from the inlet for this guess
    T = T_0
    G = G_0
    T_all[:] = 0.0
    T_all[0] = T

    # omega at the inlet (i = 0)
    c = u_0 + R * T_0 / (W * u_0)
    u = (c - math.sqrt(c**2 - 4 * R * T / W)) / 2
    rho = rho_0 * u_0 / u
    YF = YF_0 + CP * (T_0 - T) / QF
    omega = A * math.exp(-EA / (R * T)) * (rho * YF)**NU

    converged = True               # stays True unless the march fails below

    # ============ Inner loop: march from x = 0 to x = L ============
    i = 0
    while (i < N):

        G = G + (dx / LAM) * (rho_0 * u_0 * CP * G - QF * omega)
        T = T + dx * G

        if T > T_ad:               # temperature overshoots: u_0 too large (blow out)
            u_R = u_0
            converged = False
            break
        if G < 0:                  # temperature starts falling: u_0 too small (flashback)
            u_L = u_0
            converged = False
            break

        # update everything that depends on the new T (closure relations)
        c = u_0 + R * T_0 / (W * u_0)
        u = (c - math.sqrt(c**2 - 4 * R * T / W)) / 2
        rho = rho_0 * u_0 / u
        YF = YF_0 + CP * (T_0 - T) / QF
        YF = max(YF, 0.0)          # guard: tiny negative rounding error would make (rho*YF)**NU complex
        omega = A * math.exp(-EA / (R * T)) * (rho * YF)**NU

        T_all[i + 1] = T
        i = i + 1

x_reference = np.linspace(0, L, N + 1)
T_all[i + 1:] = T_all[i]        # i = last filled index when the loop stopped

if __name__ == "__main__":
    print("Flame speed s_L =", u_0, "m/s")
    print("Adiabatic flame temperature =", T_ad, "K")
    
    plt.plot(x_reference * 1e3, T_all, label = 'Analytical Solution')
    plt.xlabel('x')
    plt.ylabel('Temperature (K)')
    plt.title('BVP Solution w/ Shooting Method of 1D Flame')
    plt.legend()
    plt.show()