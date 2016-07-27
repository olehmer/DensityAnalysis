#=
Data taken from: http://www.astro.keele.ac.uk/jkt/tepcat/tepcat.html
=#

using PyPlot
using PyCall
@pyimport matplotlib.colors as COL

G = 6.67E-11 #m3/kg/s2
M_jup = 1.898E27 #mass of Jupiter in kg
R_jup = 69911000 #radius of Jupiter in m
rho_jup = 1.33 #density of Jupiter [g/cm3]
AU = 149597870700 #AU in m
M_earth = 5.972E24      #Mass of the Earth [kg]
R_earth = 6371393 #earth radius [m]
M_venus = 4.867E24 #mass of Venus [kg]
M_mars = 6.39E23 #mass of mars [kg]
M_saturn = 5.683E26 #mass of saturn [kg]
M_neptune = 1.024E26 #mass of neptune [kg]
M_uranus = 8.681E25 #mass of Uranus [kg]

function init_solar_system_planets()
  solar_system = []

  venus = Planet()
  venus.name = "Venus"
  venus.mass = M_venus / M_jup
  venus.distance = 0.72
  venus.radius = 6051133 / R_jup
  venus.period = 225
  venus.density = 5.24 / rho_jup
  venus.flux = 1366/venus.distance^2
  push!(solar_system, venus)

  earth = Planet()
  earth.name = "Earth"
  earth.mass = M_earth / M_jup
  earth.distance = 1
  earth.radius = R_earth / R_jup
  earth.period = 365
  earth.density = 5.51 / rho_jup
  earth.flux = 1366/earth.distance^2
  push!(solar_system, earth)

  mars = Planet()
  mars.name = "Mars"
  mars.mass = M_mars / M_jup
  mars.distance = 1.524
  mars.radius = 3389278 / R_jup
  mars.period = 687
  mars.density = 3.39 / rho_jup
  mars.flux = 1366/mars.distance^2
  push!(solar_system, mars)

  jupiter = Planet()
  jupiter.name = "Jupiter"
  jupiter.mass = 1
  jupiter.distance = 5.2
  jupiter.radius = 1
  jupiter.period = 4332
  jupiter.density = 1
  jupiter.flux = 1366/jupiter.distance^2
  push!(solar_system, jupiter)

  saturn = Planet()
  saturn.name = "Saturn"
  saturn.mass = M_saturn / M_jup
  saturn.distance = 9.5
  saturn.radius = 58232000 / R_jup
  saturn.period = 10585
  saturn.density = 0.687 / rho_jup
  saturn.flux = 1366/saturn.distance^2
  push!(solar_system, saturn)

  uranus = Planet()
  uranus.name = "Uranus"
  uranus.mass = M_uranus / M_jup
  uranus.distance = 19.22
  uranus.radius = 25361652 / R_jup
  uranus.period = 30660
  uranus.density = 1.27 / rho_jup
  uranus.flux = 1366/uranus.distance^2
  push!(solar_system, uranus)

  neptune = Planet()
  neptune.name = "Neptune"
  neptune.mass = M_neptune / M_jup
  neptune.distance = 30.1
  neptune.radius = 24621354 / R_jup
  neptune.period = 60225
  neptune.density = 1.64 / rho_jup
  neptune.flux = 1366/neptune.distance^2
  push!(solar_system, neptune)

  return solar_system
end


type Planet
  name::AbstractString
  mass::Float64
  gravity::Float64
  period::Float64
  eccentricity::Float64
  distance::Float64
  radius::Float64
  density::Float64
  density_e_m::Float64
  density_e_p::Float64
  temp::Float64
  stellar_temp::Float64
  stellar_mass::Float64
  stellar_rad::Float64
  flux::Float64

  function Planet()
    this = new()
    this.name = "Unknown"
    this.mass = -1.0
    this.gravity = -1.0
    this.period = -1.0
    this.eccentricity = -1.0
    this.distance = -1.0
    this.radius = -1.0
    this.density = -1.0
    this.temp = -1.0
    this.stellar_temp = -1.0
    this.stellar_mass = -1.0
    this.stellar_rad = -1.0
    this.flux = -1.0 #in W/m2
    return this
  end
end

function readdata()
  data = readdlm("exoplanet_data.txt")
  names = data[:,1] #array of system names
  stellar_temp = data[:,2:4] #[stellar temperature [K], plus error, minus error]
  stellar_mass = data[:,8:10] #mass of the star in M_sun units
  stellar_radius = data[:,11:13] #radius of the star in R_sun units
  orbital_period = data[:,20] #orbital period, [days]
  eccentricity = data[:,21:23] #[eccentricity, plus error, minus error]
  orbital_dist = data[:,24:26] #[orbit distance [AU], plus error, minus error]
  planet_masses = data[:, 27:29] #[planetary mass [MJup], plus error, minus error]
  planet_radii = data[:, 30:32] #[planetary radius [RJup], plus error, minus error]
  planet_gravity = data[:,33:35] #planetary gravity [m/s2]
  planet_density = data[:, 36:38] #[planetary density [RhoJup], plus error, minus error]
  planetary_temp = data[:, 39:41] #[planetary temp [K], plus error, minus error]
  return (names, stellar_temp, stellar_mass, stellar_radius, orbital_period, eccentricity, orbital_dist, planet_masses, planet_radii, planet_density, planetary_temp, planet_gravity)
end


function create_planet_array()
    names, stellar_temp, stellar_mass, stellar_radius, orbital_period, eccentricity, orbital_dist, planet_masses, planet_radii, planet_density, planetary_temp, planet_gravity = readdata()

    planets = []
    for i=1:length(names)
        p = Planet()
        p.name = names[i]
        p.stellar_temp = stellar_temp[i,1]
        p.stellar_mass = stellar_mass[i,1]
        p.stellar_rad = stellar_radius[i,1]
        p.gravity = planet_gravity[i,1]
        p.period = orbital_period[i]
        p.eccentricity = eccentricity[i]
        if planet_masses[i,1] == 0
          p.mass = planet_masses[i,2]
        else
          p.mass = planet_masses[i,1]
        end
        p.distance = orbital_dist[i,1]
        if planet_radii[i,1] == 0
          p.radius = planet_radii[i,2]
        else
          p.radius = planet_radii[i,1]
        end
        p.density = planet_density[i,1]
        p.density_e_p = planet_density[i,2] #positive density error
        p.density_e_m = planet_density[i,3] #negative density error
        p.temp = planetary_temp[i,1]
        if p.stellar_rad > 0 && p.stellar_temp > 0 && p.distance > 0
          p.flux = (p.stellar_rad*695700099)^2*(5.67E-8)*(p.stellar_temp)^4/(p.distance*AU)^2
        end
        push!(planets,p)
    end
    return planets
end




function findsmallplanets()
  #println("Exoplanets with known densities and mass less than 10 Earth masses.\n")
  names, stellar_temp, stellar_mass, stellar_radius, orbital_period, eccentricity, orbital_dist, planet_masses, planet_radii, planet_density, planetary_temp = readdata()

  smalls = []

  for i=1:length(names)
    emass = planet_masses[i,1]/0.00314647 #mass of earth in MJup units
    density = planet_density[i,1]*1.33 #density of jupiter in g/cm^3
    if 0.0 < emass < 10.0 && density > 0 #less than 10 earth masses
      temp = planetary_temp[i,1]
      name = names[i]
      star_temp = stellar_temp[i,1]
      dist = orbital_dist[i,1]

      p = Planet()
      p.mass = emass
      p.density = density
      p.name = name
      p.stellar_temp = star_temp
      p.temp = temp
      p.distance = dist
      p.period = orbital_period[i]
      push!(smalls, p)
      #@printf("%2d. %11s Mass=%6.3f (EM), density=%6.3f (g/cm3), planet temp=%4.0f (K), dist=%6.3f (AU), star temp=%4.0f (K)\n",
      #length(smalls), name, emass, density, temp, dist, star_temp)
    end
  end
  return smalls
end

function show_small_planets()
  planets = findsmallplanets()

  #cm = cm.get_cmap('RdYlBu')
  sc = scatter([p.mass for p in planets], [p.stellar_temp for p in planets], c=[p.density for p in planets])
  cbar = colorbar(sc)
  axis([0, 10, 2500, 6500])
  title("Super Earth Exoplanets Under 10 Earth Masses")
  cbar[:ax][:set_ylabel]("Density [g/cm3]")
  xlabel("Mass [Earth Masses]")
  ylabel("Stellar Temperature [K]")
  show()
end

function generate_rho_v_dist()
  planets = create_planet_array()

  y = []
  x = []
  z = []
  x_err_p = []
  x_err_m = []
  for p in planets
    if p.density > 0 && p.distance > 0 && p.mass > 0
      push!(y,p.flux/1366)
      push!(x,p.density*rho_jup)
      push!(z,p.mass*M_jup/M_earth)

      push!(x_err_m, p.density_e_m)
      push!(x_err_p, p.density_e_p)
    end
  end



  sc = scatter(x,y, c=z, norm=COL.LogNorm(vmin=minimum(z), vmax=maximum(z)), zorder=100)
  cbar = colorbar(sc)
  cbar[:ax][:set_ylabel]("Mass [Earth Masses]")

  errorbar(x,y, xerr=(x_err_m,x_err_p), fmt="None", marker="None", mew=0, zorder=0)

  sp = init_solar_system_planets()
  for p in sp
    plot([p.density*rho_jup], [p.flux/1366], label=p.name, marker="*", markersize=20, linestyle = "None")
  end

  legend(numpoints=1, loc="lower right")

  xscale("log")
  yscale("log")
  title("Flux vs Density")
  xlabel("Density [g/cc]")
  ylabel("Relative Stellar Flux")
  show()
end

function generate_flux_v_esc()
  planets = create_planet_array()

  y = []
  x = []
  for p in planets
    if p.flux > 0 && p.mass > 0 && p.radius > 0
      push!(y,p.flux)
      vesc = (2*G*p.mass*M_jup/(p.radius*R_jup))^0.5
      push!(x,vesc)
    end
  end

  scatter(x,y)
  xscale("log")
  yscale("log")
  title("Flux vs Escape Velocity")
  xlabel("log(Escape Velocity) [m/s]")
  ylabel("log(Stellar Flux) [W/m2]")
  show()
end

function generate_esc_v_dist()
  planets = create_planet_array()

  y = []
  x = []
  for p in planets
    if p.distance > 0 && p.mass > 0 && p.radius > 0
      vesc = (2*G*p.mass*M_jup/(p.radius*R_jup))^0.5
      push!(y,vesc)
      push!(x,p.distance)
    end
  end

  scatter(x,y)
  xscale("log")
  yscale("log")
  title("Distance vs Escape Velocity")
  xlabel("log(Distance) [AU]")
  ylabel("log(Escape Velocity) [m/s]")
  show()
end

function generate_flux_v_vesc()
  planets = create_planet_array()

  function vesc_vimp(mass, radius, period, distance)
    vesc = (2*G*mass*M_jup/(radius*R_jup))^0.5 #[m/s]
    vorb = 2*pi*distance*AU/(period*24*3600) #orbital velcity [m/s]
    venc = 0.5*vorb
    vimp = (vesc^2+venc^2)^0.5

    return (vesc, vimp)
  end

  y1 = []
  x1 = []
  z1 = []

  y2 = []
  x2 = []
  z2 = []

  y3 = []
  x3 = []
  z3 = []
  for p in planets
    if p.density > 0 && p.mass > 0 && p.radius > 0 && p.period > 0 && p.distance > 0
      vesc, vimp = vesc_vimp(p.mass, p.radius, p.period, p.distance)

      if p.density*rho_jup < 3
        push!(y1,p.flux/1366)
        push!(x1,vesc/1000)
        push!(z1,p.density*rho_jup)
      elseif p.density*rho_jup < 10
        push!(y2,p.flux/1366)
        push!(x2,vesc/1000)
        push!(z2,p.density*rho_jup)
      else
        push!(y3,p.flux/1366)
        push!(x3,vesc/1000)
        push!(z3,p.density*rho_jup)
      end

      if vesc/1000 > 175 || vesc/1000 < 14 || p.name == "Kepler-10b" || p.name == "CoRoT-07"
        annotate(p.name, xy=(vesc/1000, p.flux/1366), textcoords="data")
      end

    end
  end

  sc1 = scatter(x1,y1, c=z1, cmap="cool")#, norm=COL.LogNorm(vmin=minimum(z1), vmax=maximum(z1)))
  cbar1 = colorbar(sc1, pad = 0.003)
  cbar1[:ax][:set_ylabel]("Density [g/cc]")

  sc2 = scatter(x2,y2, c=z2, cmap="autumn")#, norm=COL.LogNorm(vmin=minimum(z2), vmax=maximum(z2)))
  cbar2 = colorbar(sc2, pad = 0.003)

  sc3 = scatter(x3,y3, c=z3, cmap="Greys")#, norm=COL.LogNorm(vmin=minimum(z3), vmax=maximum(z3)))
  cbar3 = colorbar(sc3, pad = 0.03)


  #sc = scatter(x,y, c=z, norm=COL.LogNorm(vmin=minimum(z), vmax=maximum(z)))



  #xscale("log")
  #yscale("log")
  title("Flux VS Escape Velocity")
  xlabel("Escape Velocity [km/s]")
  ylabel("Relative Stellar Heating")

  sp = init_solar_system_planets()
  for p in sp
    if p.name == "Earth" || p.name == "Mars" || p.name == "Venus"
      colr = (1,104/255,0)
      if p.name == "Mars"
        colr = "Red"
      end
      vesc, vimp = vesc_vimp(p.mass, p.radius, p.period, p.distance)
      plot([vesc/1000], [p.flux/1366], label=p.name, marker="*", markersize=10, c=colr, linestyle = "None")
      annotate(p.name, xy=(vesc/1000, p.flux/1366), textcoords="data")
    end
  end

  plot((3.0,40),(3*10.0^-2,10^4),"--")
  xlim(10.0^0.5,10^3)
  ylim(2*10.0^-1,10.0^4)

  #legend(numpoints=1, loc="upper left") #loc="upper left"

  #annotate("local max", xy=(3, 1), xytext=(0.8, 0.95))
  xscale("log")
  yscale("log")
  show()
end

function generate_vescDorb_v_rho()
  planets = create_planet_array()

  function vesc_vimp(mass, radius, period, distance)
    vesc = (2*G*mass*M_jup/(radius*R_jup))^0.5 #[m/s]
    vorb = 2*pi*distance*AU/(period*24*3600) #orbital velcity [m/s]
    venc = 0.5*vorb
    vimp = (vesc^2+venc^2)^0.5

    return (vesc, vimp)
  end

  y = []
  x = []
  z = []
  for p in planets
    if p.density > 0 && p.mass > 0 && p.radius > 0 && p.period > 0 && p.distance > 0
      vesc, vimp = vesc_vimp(p.mass, p.radius, p.period, p.distance)
      push!(y,vimp/vesc)
      push!(x,p.density*rho_jup)
      push!(z,p.flux/1366)
    end
  end


  sc = scatter(x,y, c=z, norm=COL.LogNorm(vmin=minimum(z), vmax=maximum(z)))
  cbar = colorbar(sc)
  cbar[:ax][:set_ylabel]("Stellar Flux [Earth Solar Fluxes]")
  xscale("log")
  #yscale("log")
  title("Density vs Impact Velocity / Escape Velocity")
  xlabel("Density [g/cm3]")
  ylabel("Impact/Escape Velocity")

  sp = init_solar_system_planets()
  for p in sp
    vesc, vimp = vesc_vimp(p.mass, p.radius, p.period, p.distance)
    plot([p.density*rho_jup], [vimp/vesc], label=p.name, marker="*", markersize=20, linestyle = "None")
  end

  for p in planets
    if p.name == "Kepler-10b"
      vesc, vimp = vesc_vimp(p.mass, p.radius, p.period, p.distance)
      plot([p.density*rho_jup], [vimp/vesc], label=p.name, marker="v", markersize=20, linestyle = "None")
    end
  end

  legend(numpoints=1) #loc="upper left"

  #annotate("local max", xy=(3, 1), xytext=(0.8, 0.95))

  show()
end

function generate_vescDorb_v_vesc()
  planets = create_planet_array()

  function vesc_vimp(mass, radius, period, distance)
    vesc = (2*G*mass*M_jup/(radius*R_jup))^0.5 #[m/s]
    vorb = 2*pi*distance*AU/(period*24*3600) #orbital velcity [m/s]
    venc = 0.5*vorb
    vimp = (vesc^2+venc^2)^0.5

    return (vesc, vimp)
  end

  y = []
  x = []
  z = []
  for p in planets
    if p.mass > 0 && p.radius > 0 && p.period > 0 && p.distance > 0
      vesc, vimp = vesc_vimp(p.mass, p.radius, p.period, p.distance)
      push!(y,vimp/vesc)
      push!(x,vesc)
      push!(z,p.flux/1366)


      if y[end] > 4 || p.name == "Kepler-42c"
        annotate(p.name, xy=(x[end]/1000, y[end]), textcoords="data")
      end
    end
  end


  sc = scatter(x./1000,y, c=z, norm=COL.LogNorm(vmin=minimum(z), vmax=maximum(z)))
  cbar = colorbar(sc)
  cbar[:ax][:set_ylabel]("Stellar Flux [Earth Solar Fluxes]")
  xscale("log")
  #yscale("log")
  title("Impact Velocity / Escape Velocity VS Escape Velocity")
  xlabel("Escape Velocity [km/s]")
  ylabel("Impact/Escape Velocity [km/s]")

  sp = init_solar_system_planets()
  for p in sp
    vesc, vimp = vesc_vimp(p.mass, p.radius, p.period, p.distance)
    plot([vesc/1000], [vimp/vesc], label=p.name, marker="*", markersize=10, linestyle = "None")
  end

  legend(numpoints=1) #loc="upper left"

  plot((0,1000),(5,5),"--")
  xlim(0,1000)
  #annotate("local max", xy=(3, 1), xytext=(0.8, 0.95))

  show()
end

generate_flux_v_vesc()
#generate_vescDorb_v_vesc()
#generate_rho_v_dist()
