from sympy.physics.units import speed_of_light
from sympy import (Eq, solve, sqrt)
from symplyphysics import (units, Quantity, Symbol, print_expression, dimensionless, validate_input,
    validate_output)

# Description
## Wavespeed differs in different medium. Electromagnetic wave propagation speed depends on relative permittivity and relative permeability of medium.
## Also, these characteristics usually depend on the frequency.

# Law: Vmedium = C / sqrt(ε*μ), where
## Vmedium is speed of electromagnetic wave in medium,
## C is speed of light in vacuum (it is a fundamental constant),
## ε is relative permittivity of medium,
## μ is relative permeability of medium.

wave_speed_in_medium = Symbol("wave_speed_in_medium", units.velocity)
relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)

law = Eq(wave_speed_in_medium, speed_of_light / sqrt(relative_permittivity * relative_permeability))


def print_law() -> str:
    return print_expression(law)


@validate_input(permittivity_=relative_permittivity, permeability_=relative_permeability)
@validate_output(wave_speed_in_medium)
def calculate_wavespeed(permittivity_: float, permeability_: float) -> Quantity:
    result_expr = solve(law, wave_speed_in_medium, dict=True)[0][wave_speed_in_medium]
    wavespeed_applied = result_expr.subs({
        relative_permittivity: permittivity_,
        relative_permeability: permeability_
    })
    return Quantity(wavespeed_applied)
