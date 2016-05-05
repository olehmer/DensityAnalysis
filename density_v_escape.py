"""
This file wil calculate the ocean lifetime in Ga using the analytic model
derived from the isothermal case. It will be over a range of masses and radii
to see if there is a realtion between the escape rate and the observed 1.6 Earth
radius limit on rocky exoplanets
"""

from math import exp, pi
import numpy as np
import matplotlib.pyplot as plt

#####CONSTANTS#######
k = 1.380662E-23        #Boltzmann's constant [J K^-1]
G = 6.672E-11           #Gravitational constant [N m^2 kg^-2]
sigma = 5.6704E-8       #Stefan-Boltzmann constant [W m^-2 K^-4]
m_H = 1.66E-27          #Mass of H atom [kg]
R_w = 8.314/0.018       #Specific gas constant for water vapor [J K^-1 kg^-1] = R/(0.018 kg/mol for water)
mu = 2.0                #H2 gas in a.m.u
m_gas = mu*m_H          #Molecular mass of the gas


def calculateLifetime(masses, radii, flux, albedo):
    """
    Calculates the ocean lifetime in Ga.

    Returns an array of size len(masses)xlen(radii)
    """

    life_array = np.zeros((len(masses), len(radii)))

    for i in range(0,len(masses)):
        M = masses[i]
        for j in range(0,len(radii)):
            r_s = radii[j]
            T_guess = ((0.25/sigma)*(1.0-albedo)*flux)**0.25
            #print("for M=%0.3e, r_s=%0.3e, flux=%f, T_guess=%f" % (M,r_s,flux,T_guess))

            rho_s = 1.13E11*exp(-5200.0/T_guess)/(R_w*T_guess)
            a = (k*T_guess/m_gas)**0.5
            r_c = G*M/(2.0*a**2.0)

            if r_c < r_s:
                life_array[j][i] = 0.0
            else:
                u_s = a*(r_c/r_s)**2.0*exp( -0.5+(1.0-(r_s/r_c)**2.0)*(G*M/(2.0*r_c*a**2.0)-1.0)+G*M/(r_c*a**2.0)*(1.0-r_c/r_s) )
                mass_loss = 4.0*pi*rho_s*u_s*r_s**2.0
                life_array[j][i] = M*(0.4)/mass_loss/(3.1536E7*1.0E9)


    return life_array


def plotRadVLife():
    M_earth = 5.972E24      #Mass of the Earth [kg]
    R_earth = 6.371E6       #radius of Earth [m]
    masses = np.linspace(0.0,9.0*M_earth,100)
    radii = np.linspace(R_earth,4.0*R_earth,100)

    life = calculateLifetime(masses, radii, 100000.0, 0.1)

    plt.plot(np.linspace(1,4,100), life[50][:])
    plt.yscale('log')
    plt.show()


def plotData():
    M_earth = 5.972E24      #Mass of the Earth [kg]
    R_earth = 6.371E6       #radius of Earth [m]
    masses = np.linspace(0.0,9.0*M_earth,100)
    radii = np.linspace(R_earth,4.0*R_earth,100)

    flux = 10000.0
    albedo = 0.1
    T_guess = ((0.25/sigma)*(1.0-albedo)*flux)**0.25
    print("T_guess is: %0.3f"%(T_guess))
    life = calculateLifetime(masses, radii, flux, albedo)

    density_data = np.zeros((len(masses), len(radii)))
    for i in range(0,len(masses)):
        for j in range(0,len(radii)):
            density = masses[i]/(4./3.*pi*radii[j]**3.)/1000.0 #in g/cm^3
            #print("density_data[%d][%d]=%0.3f"%(i,j,density))
            density_data[j][i] = density

    plt.xlabel("Mass [Earth Masses]")
    plt.ylabel("Radius [Earth Radius]")
    plt.plot((0, 9), (1.6, 1.6), "k:")
    CP = plt.contour(np.linspace(0,9,100), np.linspace(1,4,100), life, [0.01, 1, 100, 10000], colors='blue')
    CP2 = plt.contour(np.linspace(0,9,100), np.linspace(1,4,100), density_data, [1, 6, 10], colors='red')
    plt.clabel(CP, inline=1, fontsize=12, colors='blue', fmt="%0.2f [Ga]")
    plt.clabel(CP2, inline=1, fontsize=12, colors='red', fmt="%0.0f [g/cm3]")
    plt.show()

plotData()
