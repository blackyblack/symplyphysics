from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)
from sympy.physics.units import magnetic_constant

# Description
## Magnetic induction is a physical quantity that is a force characteristic of a magnetic field, namely, a characteristic
## of its action on moving charged particles and on bodies with a magnetic moment.

## Law is: B = mu0 * mu * H, where
## B - magnetic field induction,
## mu0 - magnetic constant,
## mu - relative permeability,
## H - magnetic field intensity.

induction = Symbol("induction", units.force / units.current / units.length)

relative_permeability = Symbol("relative_permeability", dimensionless)
intensity = Symbol("intensity", units.current / units.length)

law = Eq(induction, magnetic_constant * relative_permeability * intensity)


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permeability_=relative_permeability, intensity_=intensity)
@validate_output(induction)
def calculate_induction(relative_permeability_: float, intensity_: Quantity) -> Quantity:
    result_expr = solve(law, induction, dict=True)[0][induction]
    result_expr = result_expr.subs({
        relative_permeability: relative_permeability_,
        intensity: intensity_,
    })
    return Quantity(result_expr)
