import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import exp, pi, log
from moviepy.editor import *
import numpy as np

SAVE_TO_FILE = True

#####CONSTANTS#######
k = 1.380662E-23        #Boltzmann's constant [J K^-1]
G = 6.672E-11           #Gravitational constant [N m^2 kg^-2]
sigma = 5.6704E-8       #Stefan-Boltzmann constant [W m^-2 K^-4]
m_H = 1.66E-27          #Mass of H atom [kg]
R_w = 8.314/0.018       #Specific gas constant for water vapor [J K^-1 kg^-1] = R/(0.018 kg/mol for water)
mu = 1.0                #H2 gas in a.m.u
m_gas = mu*m_H          #Molecular mass of the gas
M_earth = 5.972E24      #Mass of the Earth [kg]
R_earth = 6.371E6       #radius of Earth [m]

def calculate_step(masses, radii, time_step, albedo=0.1, flux=2000.0):

    for i in range(0,len(masses)):
        M = masses[i]
        R = radii[i]

        core_m = M*2.0/3.0
        core_r = R/2.0

        T_guess = ((0.25/sigma)*(1.0-albedo)*flux)**0.25

        g = core_m*G/core_r**2.0
        r_s = core_r
        atmos_m = M - core_m
        p_s = g*atmos_m/(4.0*pi*r_s**2.0)
        rho_s = p_s*m_gas/k/T_guess
        a = (k*T_guess/m_gas)**0.5
        r_c = G*M/(2.0*a**2.0)

        if r_c < r_s or p_s <= 0.0:
            #the whole planet is an exosphere or airless, assume atmosphere is lost instantly
            masses[i] = 0.0
            #radii[i] = core_r
        else:
            u_s = a*(r_c/r_s)**2.0*exp( -0.5+(1.0-(r_s/r_c)**2.0)*(G*M/(2.0*r_c*a**2.0)-1.0)+G*M/(r_c*a**2.0)*(1.0-r_c/r_s) )
            mass_loss = 4.0*pi*rho_s*u_s*r_s**2.0 #in [kg/s]

            new_M = M - mass_loss*time_step
            if new_M < 0:
                masses[i] = 0.0
            else:
                masses[i] = new_M



def plotMassOverTime():
    #set the flux we will test at
    dist = 0.2
    flux = 1366.0/dist**2.0

    fig = plt.figure()
    ax = fig.add_subplot(111, xlim=(0.5, 9), ylim=(1, 4))

    mass_range = np.linspace(0.5*M_earth,9.0*M_earth,10)
    radius_range = np.linspace(R_earth,4.0*R_earth,10)
    masses = np.zeros(len(mass_range)*len(radius_range))
    radii = np.zeros(len(mass_range)*len(radius_range))
    for i in range(0, len(mass_range)):
        for j in range(0, len(radius_range)):
            masses[i*len(mass_range)+j] = mass_range[i]
            radii[i*len(mass_range)+j] = radius_range[j]


    density_data = np.zeros((100, 100))
    hd_mass = np.linspace(0.5*M_earth,9.0*M_earth,100)
    hd_radius = np.linspace(R_earth,4.0*R_earth,100)
    for i in range(0,100):
        for j in range(0,100):
            density = hd_mass[i]/(4./3.*pi*hd_radius[j]**3.)/1000.0 #in g/cm^3
            density_data[j][i] = density

    plt.xlabel("Mass [Earth Masses]")
    plt.ylabel("Radius [Earth Radius]")
    plt.title("Mass Loss at %0.2f [AU]"%dist)
    CP = plt.contour(np.linspace(0.5,9,100), np.linspace(1,4,100), density_data, [1, 6, 10], colors='red')
    plt.clabel(CP, inline=1, fontsize=12, colors='red', fmt="%0.0f [g/cm3]")
    plt.plot((0.5, 9), (1.6, 1.6), "k:")

    time_step = 3.154E13 #10 Ma in seconds

    scat, = ax.plot(masses/M_earth, radii/R_earth, linestyle='', marker='o', color='b')

    time_text = ax.text(0.6,4.05,"Time: %4.0d [Ma]"%(0))
    def animate(i):
        time_val = 0
        if i > 24 and i < 125:
            time_val = (i-24)*10
            calculate_step(masses, radii, time_step, flux=flux)
            scat.set_data(masses/M_earth, radii/R_earth)
        else:
            scat.set_data(masses/M_earth, radii/R_earth)

        if i > 124:
            time_val = 1000

        time_text.set_text("Time: %4.0d [Ma]"%(time_val))

        return scat

    ani = animation.FuncAnimation(fig, animate, frames=148, interval=100, blit=False, repeat_delay=1000, repeat=False)
    plt.show()

    if SAVE_TO_FILE:
        print("Creating video file...")
        ani.save('mass_loss.mp4', fps=24, extra_args=['-vcodec', 'libx264'])
        clip = VideoFileClip("mass_loss.mp4")
        clip.write_gif("mass_loss.gif")
        print("Cleaning up...")
        os.remove('mass_loss.mp4')

plotMassOverTime()
