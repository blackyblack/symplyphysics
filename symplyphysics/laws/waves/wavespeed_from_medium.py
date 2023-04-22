from sympy.physics.units import speed_of_light
from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, Dimensionless,
                           validate_input_symbols, validate_output_symbol)

# Description
## Wavespeed differs in different medium. Electromagnetic wave propagation speed depends on refraction factor of medium.
## Commonly refraction factor also depends on wave frequency.

# Law: Vmedium = C / n, where
## Vmedium is speed of electromagnetic wave in medium,
## C is speed of light in vacuum (it is a fundamental constant),
## n is refraction factor of medium.

wave_speed_in_medium = Symbol("wave_speed_in_medium", units.velocity)
refraction_factor = Symbol("refraction_factor", Dimensionless)

law = Eq(wave_speed_in_medium, speed_of_light / refraction_factor)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(refraction_factor_=refraction_factor)
@validate_output_symbol(wave_speed_in_medium)
def calculate_wavespeed(refraction_factor_: float) -> Quantity:
    result_expr = solve(law, wave_speed_in_medium,
                        dict=True)[0][wave_speed_in_medium]
    wavespeed_applied = result_expr.subs(refraction_factor, refraction_factor_)
    return expr_to_quantity(wavespeed_applied)
