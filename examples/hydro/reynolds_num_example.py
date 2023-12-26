from symplyphysics.laws.hydro import reynolds_number
from symplyphysics import units, convert_to, dimensionless

rho = 1000 * units.kilogram / units.meter**3
d = 0.1 * units.meter
v = 1 * units.meter / units.second
mu = 0.00089 * (units.pressure * units.time)


result = reynolds_number.calculate_reynolds_number(d, rho, v, mu)
print(
    convert_to(result, dimensionless)
)
