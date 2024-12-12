from sympy import (Eq, solve)
from sympy.physics.units import magnetic_constant
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless)

# Description
## A solenoid is a cylindrical coil consisting of a large number of turns of wire forming a helical line.
## Energy of solenoid depends on core material, intensity of magnetic field and volume of core.

## Law is: W = mu * mu0 * H^2 * V / 2, where
## W - energy,
## mu - relative permeability of the core inside of a solenoid,
## mu0 - magnetic constant,
## H - magnetic field intensity,
## V - volume of solenoid.

# Links: Physics LibreTexts, derivable from 14.4.1 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/14%3A_Inductance/14.04%3A_Energy_in_a_Magnetic_Field>

# NOTE replace `H` with `B`?

energy = Symbol("energy", units.energy)

relative_permeability = Symbol("relative_permeability", dimensionless)
intensity = Symbol("intensity", units.current / units.length)
volume = Symbol("volume", units.volume)

law = Eq(energy, relative_permeability * magnetic_constant * intensity**2 * volume / 2)


@validate_input(relative_permeability_=relative_permeability, intensity_=intensity, volume_=volume)
@validate_output(energy)
def calculate_energy(relative_permeability_: float, intensity_: Quantity,
    volume_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({
        relative_permeability: relative_permeability_,
        intensity: intensity_,
        volume: volume_
    })
    return Quantity(result_expr)
