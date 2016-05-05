#=
Data taken from: http://www.astro.keele.ac.uk/jkt/tepcat/tepcat.html
=#

using PyPlot

type Planet
  name::AbstractString
  mass::Float64
  period::Float64
  eccentricity::Float64
  distance::Float64
  radius::Float64
  density::Float64
  temp::Float64
  stellar_temp::Float64

  function Planet()
    this = new()
    this.name = "Unknown"
    this.mass = -1.0
    this.period = -1.0
    this.eccentricity = -1.0
    this.distance = -1.0
    this.radius = -1.0
    this.density = -1.0
    this.temp = -1.0
    this.stellar_temp = -1.0
    return this
  end
end

function readdata()
  data = readdlm("exoplanet_data.txt")
  names = data[:,1] #array of system names
  stellar_temp = data[:,2:4] #[stellar temperature [K], plus error, minus error]
  orbital_period = data[:,20] #orbital period, [days]
  eccentricity = data[:,21:23] #[eccentricity, plus error, minus error]
  orbital_dist = data[:,24:26] #[orbit distance [AU], plus error, minus error]
  planet_masses = data[:, 27:29] #[planetary mass [MJup], plus error, minus error]
  planet_radii = data[:, 30:32] #[planetary radius [RJup], plus error, minus error]
  planet_density = data[:, 36:38] #[planetary density [RhoJup], plus error, minus error]
  planetary_temp = data[:, 39:41] #[planetary temp [K], plus error, minus error]
  return (names, stellar_temp, orbital_period, eccentricity, orbital_dist,
          planet_masses, planet_radii, planet_density, planetary_temp)
end

function findsmallplanets()
  #println("Exoplanets with known densities and mass less than 10 Earth masses.\n")
  names, stellar_temp, orbital_period, eccentricity, orbital_dist, planet_masses, planet_radii, planet_density, planetary_temp = readdata()

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
