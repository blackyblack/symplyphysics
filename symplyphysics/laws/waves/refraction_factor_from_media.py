from sympy import (Eq, solve, sqrt, symbols)
from symplyphysics import print_expression

# Description
## How media refracts electromagnetical waves depends on how this media transfers electrical and magnetical fields.

# Law: n = sqrt(epsilon * mu), where
## n is refraction factor of media,
## epsilon is relative dielectric permeability,
## mu is relative magnetic permeability.

refraction_factor = symbols("refraction_factor")
relative_dielectric_permeability = symbols("relative_dielectric_permeability")
relative_magnetic_permeability = symbols("relative_magnetic_permeability")

law = Eq(refraction_factor, sqrt(relative_dielectric_permeability * relative_magnetic_permeability))


def print() -> str:
    return print_expression(law)


def calculate_refraction_factor(relative_dielectric_permeability_: float,
    relative_magnetic_permeability_: float) -> float:
    result_expr = solve(law, refraction_factor, dict=True)[0][refraction_factor]
    factor_applied = result_expr.subs({
        relative_dielectric_permeability: relative_dielectric_permeability_,
        relative_magnetic_permeability: relative_magnetic_permeability_
    })
    return factor_applied
