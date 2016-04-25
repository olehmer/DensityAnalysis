
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

  return (names, stellar_temp, orbital_period, eccentricity, orbital_dist, planet_masses, planet_radii, planet_density, planetary_temp)
end

function findsmallplanets()
  println("Exoplanets with known densities and mass less than 10 Earth masses.\n")
  names, stellar_temp, orbital_period, eccentricity, orbital_dist, planet_masses, planet_radii, planet_density, planetary_temp = readdata()

  smalls = []

  for i=1:length(names)
    emass = planet_masses[i,1]/0.00314647 #mass of earth in MJup units
    density = planet_density[i,1]*1.33 #density of jupiter in g/cm^3
    if 0.0 < emass < 10.0 && density > 0 #less than 10 earth masses
      push!(smalls, i)
      temp = planetary_temp[i,1]
      name = names[i]
      @printf("%11s Mass=%6.3f (EM), density=%6.3f (g/cm3), temp=%8.3f (K)\n", name, emass, density, temp)
    end
  end

  @printf("Total of %d exoplanets known with these characteristics",length(smalls))
end

findsmallplanets()
