#=
Data taken from: http://www.astro.keele.ac.uk/jkt/tepcat/tepcat.html
=#

using PyPlot

G = 6.67E-11 #m3/kg/s2
M_jup = 1.898E27 #mass of Jupiter in kg
R_jup = 69911000 #radius of Jupiter in m
AU = 149597870700 #AU in m

type Planet
  name::AbstractString
  mass::Float64
  gravity::Float64
  period::Float64
  eccentricity::Float64
  distance::Float64
  radius::Float64
  density::Float64
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
        p.mass = planet_masses[i,1]
        p.distance = orbital_dist[i,1]
        p.radius = planet_radii[i,1]
        p.density = planet_density[i,1]
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

function generate_vescDorb_v_rho()
  planets = create_planet_array()

  y = []
  x = []
  for p in planets
    if p.density > 0 && p.mass > 0 && p.radius > 0 && p.period > 0 && p.distance > 0
      vesc = (2*G*p.mass*M_jup/(p.radius*R_jup))^0.5 #[m/s]
      vimp = 2*pi*p.distance*AU/(p.period*24*3600) #impact veloctiy ~ orbital velcity [m/s]
      push!(y,vesc/vimp)
      push!(x,p.density)
    end
  end

  scatter(x,y)
  xscale("log")
  yscale("log")
  title("Density vs Impact Velocity / Escape Velocity")
  xlabel("log(Density) [kg/m3]")
  ylabel("Impact/Escape Velocity")
  show()
end

generate_vescDorb_v_rho()
