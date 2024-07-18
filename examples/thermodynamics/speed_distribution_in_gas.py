#!/usr/bin/env python3

from sympy import solve, symbols, integrate
from symplyphysics import Quantity, units, convert_to, quantities
from symplyphysics.laws.chemistry import molar_mass_via_mass_of_single_particle as molar_mass_law
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution

# Description
## In oxygen, whose molar mass is 0.032 kg/mol, at room temperature, what fraction of the molecules
## have speeds in the interval 599 to 601 m/s?

molar_mass = symbols("molar_mass", positive=True)
initial_speed, final_speed = symbols("initial_speed final_speed")

values = {
    molar_mass: Quantity(0.032 * units.kilogram / units.mole),
    initial_speed: Quantity(599 * units.meter / units.second),
    final_speed: Quantity(601 * units.meter / units.second),
}

molecule_mass = (solve(molar_mass_law.law,
    molar_mass_law.particle_mass)[0].subs(molar_mass_law.molar_mass, molar_mass))

speed_distribution_expr = speed_distribution.law.rhs.subs({
    speed_distribution.particle_mass: molecule_mass,
    speed_distribution.equilibrium_temperature: quantities.standard_laboratory_temperature,
})

fraction_of_molecules = integrate(
    speed_distribution_expr,
    (speed_distribution.particle_speed, initial_speed, final_speed),
)

fraction_of_molecules_ = fraction_of_molecules.subs(values)

fraction_of_molecules_value = convert_to(Quantity(fraction_of_molecules_), units.percent).evalf(3)
initial_speed_value = convert_to(values[initial_speed], units.meter / units.second)
final_speed_value = convert_to(values[final_speed], units.meter / units.second)

print(f"The fraction of oxygen molecules with speeds between {initial_speed_value} m/s and "
    f"{final_speed_value} m/s is {fraction_of_molecules_value}%.")
