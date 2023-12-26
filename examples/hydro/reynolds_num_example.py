from symplyphysics.laws.hydro import reynolds_number
from symplyphysics import units, convert_to, dimensionless, Quantity


density = Quantity(1000 * units.kilogram / units.meter**3)
diameter = Quantity(0.1 * units.meter)
velosity = Quantity(1 * units.velocity)
dynamic_viscosity = Quantity(0.000894 * units.pascal * units.second)


result = reynolds_number.calculate_reynolds_number(
    diameter_=diameter,
    density_=density,
    velocity_=velosity,
    dynamic_viscosity_=dynamic_viscosity
).evalf(6)

print(f"The formula for Reynolds number is: {reynolds_number.print_law()}")
print(
    "res=", convert_to(result, dimensionless)
)
