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

h = 0.1
n = 1.5
polytrope = star.polytropic_solver(n, h, init)


#FOR QUESTION 2A
m = 1*1.989*10**33
rho_c = 49.4 #how to get this
mu = (2*(0.715/1)+3*(0.271/4)+8*(0.014/14))**-1
gas = "yes"
radiation = "no"


center_2a, profile_2a, beta_2a = star.stellar_estimates(polytrope, n, m, rho_c, mu, gas, radiation)

#contribution for total energy
pp = 1
cno = 0
he = 0

shells_2a = star.stellar_burning(pp, cno, he, polytrope, center_2a, profile_2a)
print("QUESTION 2A")
print("Radius:", "{:e}".format(profile_2a[0][-1])+" cm", profile_2a[0][-1]/r_sun, " Rsun")
print("Central Density:", "{:e}".format(center_2a[0])+" g cm^-3")
print("{:e}".format(center_2a[1])+" g cm^-1 s^-2")
print("Central Temperature:", "{:e}".format(center_2a[2])+" K")
print("Luminosty:", "{:e}".format(sum(shells_2a[3]))+" ergs s^-1", sum(shells_2a[3])/l_sun, " Lsun")
print("beta:", beta_2a)

M_coord = []
M_internal = shells_2a[0][0]
L = []
L_internal = shells_2a[3][0]

for i, j in zip(shells_2a[0], shells_2a[3]):

    M_internal = M_internal + i
    L_internal = L_internal + j
    M_coord.append(M_internal)
    L.append(L_internal)

fig, ax = plt.subplots()
ax.scatter((M_coord/M_internal), L, s=0.2, c="red")
ax.set_xlabel("Mass Coordinate")
ax.set_ylabel("Luminosity erg $s^{-1}$")


fig1, ax1 = plt.subplots()
ax1.scatter(M_coord/M_internal, shells_2a[2], s=0.2, c="red")
ax1.set_xlabel("Mass Coordinate")
ax1.set_ylabel("Specific Energy Generation Rate erg $s^{-1}$ $g^{-1}$")








#FOR QUESTIONS 2B
print("QUESTION 2B")
h = 0.1
n = 3
polytrope = star.polytropic_solver(n, h, init)
xi = polytrope[0][-1]
zi = polytrope[2][-1]
fig2, ax2 = plt.subplots()

ax2.set_xlabel("Mass Coordinate")
ax2.set_ylabel("Temperature K")

fig3, ax3 = plt.subplots()

ax3.set_xlabel("Mass Coordinate")
ax3.set_ylabel("Density g$cm^{-3}$")

m = 100000*1.989*10**33
temps = (2*10**7, 2.5*10**7, 3*10**7, 3.5*10**7, 4*10**7)

for t in temps:
    print("Temperature:", t, " K")
    #extra bit for the question
    temp_c = t
    pressure_c = ((4*sig)/(c))*(temp_c**4)/3
    rho_c = (((pressure_c*(n+1))/(4*np.pi*G))*((-m)/(4*np.pi*(xi**2)*zi))**(-2/3))**(3/4)
    mu = (2*(0.715/1)+3*(0.271/4)+8*(0.014/14))**-1
    gas = "no"
    radiation = "yes"
    center_2b, profile_2b, beta_2b = star.stellar_estimates(polytrope, n, m, rho_c, mu, gas, radiation)
    pp = 0
    cno = 1
    he = 0
    shells_2b = star.stellar_burning(pp, cno, he, polytrope, center_2b, profile_2b)
    print("Radius:", "{:e}".format(profile_2b[0][-1])+" cm", profile_2b[0][-1]/r_sun, " Rsun")
    print("Central Density:", "{:e}".format(center_2b[0])+" g cm^-3")
    print("Central Temperature:", "{:e}".format(center_2b[2])+" K")
    print("Luminosty:", "{:e}".format(sum(shells_2b[3]))+" ergs s^-1", sum(shells_2b[3])/l_sun, " Lsun")
    print("beta:", beta_2b)

    M_coord = []
    M_internal = shells_2b[0][0]

    for i in shells_2b[0]:

        M_internal = M_internal + i
        M_coord.append(M_internal)

    ax2.plot((M_coord/M_internal), profile_2b[3], lw=1, label="{:e}".format(t)+"K")
    ax3.plot(M_coord/M_internal, profile_2b[1], lw=1, label="{:e}".format(t)+"K")


ax2.legend()
ax3.legend()
