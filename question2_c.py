import star_stuff as star
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns

print("QUESTION 2C")

pp = []
cno = []
he3 = []
temp = []

rho = 49.4
fig, ax = plt.subplots()

first_temps = np.arange(0.01*10**7, 5*10**8, 0.01*10**7)
more_temps = np.arange(5*10**8, 5*10**9, 1*10**8)

for t in np.hstack((first_temps, more_temps)):

    pp.append(star.e_pp(t, rho))
    cno.append(star.e_cno(t, rho))
    he3.append(star.e_he3(t, rho))
    temp.append(t)






ax.plot(temp, pp, lw=0.8, c="green", label="PP $\\rho_{c} = 49.4$", alpha=0.7)
ax.plot(temp, cno, lw=0.8, c="blue", label="CNO $\\rho_{c} = 49.4$", alpha=0.7)
ax.plot(temp, he3, lw=0.8, c="red", label="Triple Alpha $\\rho_{c} = 49.4$", alpha=0.7)

pp = []
cno = []
he3 = []
temp = []

rho = 3.29*10**-3

first_temps = np.arange(0.01*10**7, 5*10**8, 0.01*10**7)
more_temps = np.arange(5*10**8, 5*10**9, 1*10**8)

for t in np.hstack((first_temps, more_temps)):

    pp.append(star.e_pp(t, rho))
    cno.append(star.e_cno(t, rho))
    he3.append(star.e_he3(t, rho))
    temp.append(t)

ax.plot(temp, pp, lw=0.8, c="darkgreen", label="PP $\\rho_{c} = 0.003$", alpha=0.7, linestyle=":")
ax.plot(temp, cno, lw=0.8, c="darkblue", label="CNO $\\rho_{c} = 0.003$", alpha=0.7, linestyle=":")
ax.plot(temp, he3, lw=0.8, c="darkred", label="Triple Alpha $\\rho_{c} = 0.003$", alpha=0.7, linestyle=":")

ax.axvline(1.8522e7, ymin=0, ymax=100000, c="black", label="T = 1.8522e7", lw=0.5)

ax.set_title("Reaction strengths")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(0.0001,1000)
ax.set_xlim(3*10**6,11**8.4)
ax.set_xlabel("Log Tc")
ax.set_ylabel("Log e nuc")
ax.legend(loc="upper right", prop={'size': 6})
