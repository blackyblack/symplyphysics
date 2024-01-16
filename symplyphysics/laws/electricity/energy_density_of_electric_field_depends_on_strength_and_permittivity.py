from sympy import (Eq, solve, exp)
from sympy.physics.units import electric_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

# Description
## In physics, energy density or volumic energy is the amount of energy stored in a given system or region of space per unit volume.
## It is sometimes confused with energy per unit mass which is properly called massic energy or gravimetric energy density.
## The energy of the electric field depends on the permittivity of the medium and the intensity of the electric field.

## Law is: w = (e0 * e * E^2) / 2, where
## w - the energy density of the electric field,
## e0 - electrical constant,
## e - dielectric constant of the material,
## E - the intensity of electric field.

energy_density = Symbol("energy_density", units.energy / units.volume)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

law = Eq(
    energy_density,
    (electric_constant * relative_permittivity * electric_intensity**2) / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity, electric_intensity_=electric_intensity)
@validate_output(energy_density)
def calculate_energy_density(relative_permittivity_: Quantity, electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, energy_density, dict=True)[0][energy_density]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        electric_intensity: electric_intensity_,
    })
    return Quantity(result_expr)
