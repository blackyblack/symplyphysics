from typing import Union
from sympy import (Eq, solve, cos)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, angle_type)

# Description
## Magnetic flux is the flux of a magnetic induction vector through a certain surface.
## Law is: Φ = B * S * cos(a), where
## Φ - flux,
## B - induction,
## S - area,
## a - angle between the normal of pad and magnetic induction.

flux = Symbol("flux", units.magnetic_flux)

induction = Symbol("induction", units.magnetic_density)
area = Symbol("area", units.area)
angle = Symbol("angle", angle_type)

law = Eq(flux, induction * area * cos(angle))


def print_law() -> str:
    return print_expression(law)


@validate_input(induction_=induction, area_=area, angle_=angle)
@validate_output(flux)
def calculate_flux(induction_: Quantity, area_: Quantity, angle_: Union[float, Quantity]) -> Quantity:
    result_flux_expr = solve(law, flux, dict=True)[0][flux]
    result_expr = result_flux_expr.subs({induction: induction_, area: area_, angle: angle_})
    return Quantity(result_expr)
