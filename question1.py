import star_stuff as star
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns

#useful constants
sig = 5.6704*10**(-5)
c = 2.98*10**10
G = 6.67*10**-8


r_sun = 69.634*10**9
m_sun = 1.989*10**33
l_sun = 3.846*10**33

x_init = 0  #zeta at center of star (radius)
y_init = 1  #theta at center of star (temperature)
z_init = 0 # d_theta/d_zeta at center of star
dz_dx_init = -1/3  #d^2_theta/d^2_zeta and center of star from expansion and limit
init = (x_init, y_init, z_init, dz_dx_init)



#FOR QUESTION 1
h = 0.1
n = 1.5
polytrope = star.polytropic_solver(n, h, init)


#STAR 1a
m = 1*1.989*10**33
rho_c = 160 # cgs
mu = 0.61
gas = "yes"
radiation = "no"

center_1a, profile_1a, beta_1a = star.stellar_estimates(polytrope, n, m, rho_c, mu, gas, radiation)
print("QUESTION 1A")
print("{:e}".format(profile_1a[0][-1])+" cm", profile_1a[0][-1]/r_sun, " Rsun")
print("{:e}".format(center_1a[1])+" g cm^-1 s^-2")
print("{:e}".format(center_1a[2])+" K")
print("beta:", beta_1a)






#STAR 1c
m = 15*1.989*10**33
rho_c = 6 # cgs
mu = (2*(0.715/1)+3*(0.271/4)+8*(0.014/14))**-1
gas = "yes"
radiation = "yes"

center_1c, profile_1c, beta_1c = star.stellar_estimates(polytrope, n, m, rho_c, mu, gas, radiation)
print("QUESTION 1C")
print("{:e}".format(profile_1c[0][-1])+" cm", profile_1c[0][-1]/r_sun, " Rsun")
print("{:e}".format(center_1c[1])+" g cm^-1 s^-2")
print("{:e}".format(center_1c[2])+" K")
print("beta:", beta_1c)
