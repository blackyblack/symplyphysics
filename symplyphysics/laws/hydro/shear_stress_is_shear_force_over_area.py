from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Shear stress is the component of stress coplanar with the material cross section on which it acts,
## that is, it lies in the same plane as the area it is applied to. It is opposed to normal stress,
## which arises from the force vector component perpendicular to the material cross section.

# Law: tau = F/A
## tau - shear stress
## F - force applied
## A - area of the material face parallel to the applied force vector

shear_stress = Symbol("shear_stress", units.pressure)
force_applied = clone_symbol(symbols.dynamics.force, "force_applied")
area = Symbol("area", units.area)

law = Eq(shear_stress, force_applied / area)


def print_law() -> str:
    return print_expression(law)


@validate_input(force_applied_=force_applied, area_=area)
@validate_output(shear_stress)
def calculate_shear_stress(force_applied_: Quantity, area_: Quantity) -> Quantity:
    solved = solve(law, shear_stress)[0]
    result = solved.subs({
        force_applied: force_applied_,
        area: area_,
    })
    return Quantity(result)
