
import numpy as np
import math
import sympy as sy
from scipy.optimize import minimize



def polytropic_solver(n, h, init):

    #SOLVE POLYTROPE WITH LARGE STEP SIZE
    data = runge_kutta(func_1, func_2, func_3, init, n, h)
    #REFINE ESTIMATE FOR EDGE VALUES
    edge = star_edge(data, n, h)

    results = interpolate_solution(data, edge)

    return results



def constants(n, m, rho_c, xi, zi):
    G = 6.67*10**-8
    alpha = (((-4*np.pi*(rho_c)*(xi**2)*(zi))**-1)*m)**(1/3)
    k = ((alpha**2)*(4*G*np.pi)*(rho_c**(((n-1)/n))))/(n+1)
    return alpha, k



def stellar_profile(polytrope, center, m, n):
    alpha, k = constants(n, m, center[0], polytrope[0][-1], polytrope[2][-1])
    rho = center[0]*polytrope[1]**(n)
    pressure = center[1]*polytrope[1]**(1+n)
    temp = center[2]*polytrope[1]
    radius = alpha*polytrope[0]
    return (radius, rho, pressure, temp)


def stellar_estimates(polytrope, n, m, rho_c, mu, gas, radiation):

    r = 8.314462*(10**7)
    sig = 5.6704*10**(-5)
    c = 2.98*10**10

    xi = polytrope[0][-1]
    zi = polytrope[2][-1]
    alpha, k = constants(n, m, rho_c, xi, zi)
    pressure_c = k*rho_c**((n+1)/n)
    if gas == "yes":
        if radiation == "yes":

            T = sy.Symbol('T')
            eqn = (r*T*rho_c)/mu + ((4*sig)/(c))*(T**4)/3 - pressure_c
            temp_c = sy.solvers.solve(eqn)[1]
            center = (rho_c, pressure_c, temp_c)

        else:

            temp_c = (pressure_c*mu)/(rho_c*r)
            center = (rho_c, pressure_c, temp_c)

    if gas == "no":
        if radiation == "yes":

            T = sy.Symbol('T')
            eqn = ((4*sig)/(c))*(T**4)/3 - pressure_c
            temp_c = sy.solvers.solve(eqn)[1]
            center = (rho_c, pressure_c, temp_c)

        else:
            print("gas and radiaition havnt been specified corrrectly")
            center = (rho_c, pressure_c, np.nan)


    profile = stellar_profile(polytrope, center, m, n)
    beta = (r*center[2]*center[0]/mu)/(pressure_c)

    return center, profile, beta


def interpolate_solution(data, edge):

        x_temp = np.hstack((data[0], edge[0]))
        y_temp = np.hstack((data[1], edge[1]))
        z_temp = np.hstack((data[2], edge[2]))
        dz_dx_temp = np.hstack((data[3], edge[3]))

        x = np.arange(0, edge[0], 0.001)
        y = np.interp(x, x_temp, y_temp)
        z = np.interp(x, x_temp, z_temp)
        dz = np.interp(x, x_temp, dz_dx_temp)

        return (x, y, z, dz)


def func_1(x, y, z, n):
    dy_dx = z
    return dy_dx


def func_2(x, y, z, n):

    dz_dx = (-1/(x**2))*(2*z*x + (y**n)*(x**2))

    return dz_dx

def func_3(z_init):
    dz_dx = z_init
    return dz_dx


def runge_kutta(func_1, func_2, func_3, init, n, h):

    x, y, z, dz_dx = [], [], [], []
    x_n, y_n, z_n, dz_dx_init = init[0], init[1], init[2], init[3]


    x.append(x_n), y.append(y_n), z.append(z_n), dz_dx.append(dz_dx_init)

    counter, threshhold = 0, 0
    while threshhold < 1:

        if counter == 0:
            k1 = h*func_1(x_n, y_n, z_n, n)
            l1 = h*func_3(dz_dx_init)
            counter = counter + 1

        else:
            k1 = h*func_1(x_n, y_n, z_n, n)
            l1 = h*func_2(x_n, y_n, z_n, n)

        k2 = h*func_1(x_n + h/2, y_n + k1/2, z_n + l1/2, n)
        l2 = h*func_2(x_n + h/2, y_n + k1/2, z_n + l1/2, n)

        k3 = h*func_1(x_n + h/2, y_n + k2/2, z_n + l2/2, n)
        l3 = h*func_2(x_n + h/2, y_n + k2/2, z_n + l2/2, n)

        k4 = h*func_1(x_n + h, y_n + k3, z_n + l3, n)
        l4 = h*func_2(x_n + h, y_n + k3, z_n + l3, n)

        x_1 = x_n + h
        y_1 = y_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        z_1 = z_n + (1/6)*(l1 + 2*l2 + 2*l3 + l4)

        if isinstance(z_1, complex):
            threshhold=1

        if isinstance(z_1, complex):
            threshhold=1
        elif  y_1 < 0 or math.isnan(y_1):
            threshhold=1
        else:
            x.append(x_1), y.append(y_1), z.append(z_1)
            dz_dx.append(func_2(x_1, y_1, z_1, n))
            x_n, y_n, z_n = x_1, y_1, z_1

    return (np.asarray(x)[:-1], np.asarray(y)[:-1], np.asarray(z)[:-1], np.asarray(dz_dx)[:-1])


def star_edge(data, n, h):

    h_new = h/100
    x_init = data[0][-2]
    y_init = data[1][-2]
    z_init = data[2][-2]
    dz_dx_init = data[3][-2]

    init = (x_init, y_init, z_init, dz_dx_init)

    i = 0
    while i < 5:
        i = i + 1

        data_new = runge_kutta(func_1, func_2, func_3, init, n, h_new)

        h_new = h_new/10

        try:
            x_init = data_new[0][-2]
            y_init = data_new[1][-2]
            z_init = data_new[2][-2]
            dz_dx_init = data_new[3][-2]
            init = (x_init, y_init, z_init, dz_dx_init)


        except:
            break

    final = (data_new[0][-1], data_new[1][-1], data_new[2][-1], data_new[3][-1])

    return final

def e_pp(t, rho):

    t = float(t)
    rho = float(rho)

    x1 = 0.715

    t9 = t/(10**9)

    g11 = (1+3.82*t9 + 1.51*t9**2 + 0.144*t9**3 - 0.0114*t9**4)

    e = (2.57*10**4)*(1)*(2)*g11*(rho)*(x1**2)*(t9**(-2/3))*(np.exp(-3.381/((t9)**(1/3))))

    return e


def e_cno(t, rho):

    t = float(t)
    rho = float(rho)

    x1 = 0.715
    x_cno = 0.014

    g14 = 1-2*(t/(10**9)) + 3.41*(t/10**9)**2 - 2.43*(t/10**9)**3

    e = (8.24*10**25)*g14*(x_cno)*(x1)*(rho)*((t/(10**9))**(-2/3))*np.exp(-15.231*((t/(10**9))**(-1/3))-((t/(10**9))/0.8)**2)

    return e

def e_he3(t, rho):

    t = float(t)
    rho = float(rho)

    x4 = 0.271

    t9 = t/(10**9)

    term_1 = (2.43*(10**9)*(t9**(-2/3))*np.exp(-13.490*(t9**(-1/3)) - (t9/0.15)**2)*(1+74.5*t9)+6.09*(10**5)*(t9**(-3/2))*np.exp(-1.054/t9))
    term_2 = (2.76*(10**7)*(t9**-(2/3))*np.exp(-23.570*(t9**(-1/3))-(t9/0.4)**2)*(1+5.47*t9 + 326*t9**2)+130.7*(t9**(-3/2))*np.exp(-3.338/t9)+2.51*(10**4)*(t9**(-3/2))*np.exp(-20.307/t9))

    e = 6.272*(rho**2)*(x4**3)*(1+0.01589*t9**(-0.65))*term_1*term_2

    return e

def stellar_burning(pp, cno, he, polytrope, center, profile):

    radius, rho, pressure, temp = profile
    rho_c, pressure_c, temp_c = center

    x, y, z, dz_dx = polytrope

    de = np.empty([np.shape(x)[0]])
    dv = np.empty([np.shape(x)[0]])
    dm = np.empty([np.shape(x)[0]])
    dl = np.empty([np.shape(x)[0]])

    dpp = np.empty([np.shape(x)[0]])
    dcno = np.empty([np.shape(x)[0]])
    dhe = np.empty([np.shape(x)[0]])

    dpp[0] = e_pp(temp_c, rho_c)
    dcno[0] = e_cno(temp_c, rho_c)
    dhe[0] = e_he3(temp_c, rho_c)

    de[0] = pp*e_pp(temp_c, rho_c) + cno*e_cno(temp_c, rho_c) + he*e_he3(temp_c, rho_c)
    dv[0] = 0
    dm[0] = 0
    dl[0] = 0

    for i in range(1, np.shape(x)[0]):

        j = i-1

        dv[i] = ((4*np.pi)/3)*(radius[i]**3 - radius[j]**3)

        dm[i] = dv[i]*((rho[j] + rho[i])/2)

        de[i] = pp*e_pp(temp[i], rho[i]) + cno*e_cno(temp[i], rho[i]) + he*e_he3(temp[i], rho[i])
        dpp[i] = e_pp(temp[i], rho[i])
        dcno[i] = e_cno(temp[i], rho[i])
        dhe[i] = e_he3(temp[i], rho[i])

        dl[i] = dm[i]*(de[i]+de[j])/2


    return (dm, dv, de, dl, dpp, dcno, dhe)
