from sympy import (
    Eq,
    solve,
    pi,
    sqrt,
)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## The characteristic resistance of a wave is a value determined by the ratio of the transverse component
## of the electric field strength to the transverse component of the magnetic field strength of a traveling wave.

## Law is: Z = 120 * pi * sqrt(mur / er), where
## Z - characteristic resistance,
## mur - relative permeability of the insulator material,
## er - relative permittivity of insulating material.

# TODO find link

resistance = Symbol("resistance", units.impedance)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)

impedance_vacuum = Quantity(120 * pi * units.ohm)

law = Eq(resistance, impedance_vacuum * sqrt(relative_permeability / relative_permittivity))


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability)
@validate_output(resistance)
def calculate_resistance(relative_permittivity_: float, relative_permeability_: float) -> Quantity:
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
    })
    return Quantity(result_expr)
