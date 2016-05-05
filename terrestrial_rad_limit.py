import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import exp, pi, log
from moviepy.editor import *
import os

SAVE_TO_FILE = False

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

def calculate_step(masses, radii, cores, time_step, albedo=0.1, flux=2000.0):

    for i in range(0,len(masses)):
        M = masses[i]
        R = radii[i]

        core_m, core_r = cores[i]
        #core_m = M_earth
        #core_r = R_earth

        atmos_m = M - core_m
        core_sa = 4.0*pi*core_r**2.0 #surface area of the core

        T_guess = ((0.25/sigma)*(1.0-albedo)*flux)**0.25

        g = core_m*G/core_r**2.0
        r_s = core_r
        p_s = g*atmos_m/core_sa
        rho_s = p_s*m_gas/k/T_guess
        a = (k*T_guess/m_gas)**0.5
        r_c = G*M/(2.0*a**2.0)

        if r_c < r_s or p_s <= 0.0:
            #the whole planet is an exosphere or airless, assume atmosphere is lost instantly
            masses[i] = core_m
            radii[i] = core_r
        else:
            H = (k*T_guess/m_gas/g)
            z_H = log(10.0*p_s) #pressure at 0.1 bar
            est_R = z_H*H + core_r
            u_s = a*(r_c/r_s)**2.0*exp( -0.5+(1.0-(r_s/r_c)**2.0)*(G*M/(2.0*r_c*a**2.0)-1.0)+G*M/(r_c*a**2.0)*(1.0-r_c/r_s) )
            mass_loss = 4.0*pi*rho_s*u_s*r_s**2.0 #in [kg/s]

            new_M = M - mass_loss*time_step
            if new_M < core_m:
                masses[i] = core_m
            else:
                #print("M=%0.03f, R=%0.03f, updating mass!"%(new_M/M_earth,R/R_earth))
                masses[i] = new_M

            if est_R < core_r:
                radii[i] = core_r
            else:
                #print("New R=%0.3f"%(est_R/R_earth))
                radii[i] = est_R
    #return (result_m, result_r)


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


    cores = [] #tracks mass of the core
    for i in range(0, len(masses)):
        #estimate the core size, assume the core is 6 g/cm3
        #assume the core % scales linearly with the initial radius (from 100% core at R=1, to 50% core at R=4)
        M = masses[i]
        R = radii[i]
        core_percent = exp(-5.0*(R-radii[0])/(radii[-1]-radii[0]))
        if core_percent == 0.0:
            core_percent = 0.01
        print("M=%0.3f has core_percent=%3.2f"%(M/M_earth, core_percent*100.0))
        #core_rho = 6000.0 #density in kg/m3 = 6 g/cm3
        core_r = R_earth
        core_m = core_percent*M
        #core_r = (3.0/4.0/pi*(core_m/core_rho))**(1.0/3.0)
        cores.append( (core_m,core_r) )


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

    time_step = 3.154E14 #10 Ma in seconds

    scat, = ax.plot(masses/M_earth, radii/R_earth, linestyle='', marker='o', color='b')

    time_text = ax.text(0.6,4.05,"Time: %4.0d [Ma]"%(0))
    def animate(i):
        time_text.set_text("Time: %4.0d [Ma]"%(i*10+10))
        calculate_step(masses, radii, cores, time_step, flux=flux)
        scat.set_data(masses/M_earth, radii/R_earth)
        return scat

    ani = animation.FuncAnimation(fig, animate, frames=100, interval=100, blit=False, repeat_delay=1000, repeat=False)
    plt.show()

    if SAVE_TO_FILE:
        print("Creating video file...")
        ani.save('mass_loss.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
        clip = VideoFileClip("mass_loss.mp4")
        clip.write_gif("mass_loss.gif")
        print("Cleaning up...")
        os.remove('mass_loss.mp4')

plotMassOverTime()
