from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Radiant exitance (a.k.a. radiant emittance) is the radiant flux emitted by a surface per unit area.

# M_e = d(Phi_e)/dA
## M_e - radiant exitance of a surface (e stands for "energetical" to avoid confusion with photometric quantities)
## Phi_e - radiant flux (a.k.a. radiant power) emitted
## d/dA - derivative w.r.t. area of the surface

radiant_exitance = Symbol("radiant_exitance", units.power / units.area)
radiant_flux = Function("radiant_flux", units.power)
surface_area = Symbol("surface_area", units.area)

definition = Eq(radiant_exitance, Derivative(radiant_flux(surface_area), surface_area))


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    radiant_flux_=radiant_flux,
    surface_area_=surface_area,
)
@validate_output(radiant_exitance)
def calculate_radiant_exitance(
    radiant_flux_: Quantity,
    surface_area_: Quantity,
) -> Quantity:
    """\
    Calculate radiant exitance for a surface small enough that \
    the radiant flux passed through it is constant.\
    """

    flux_function = radiant_flux_ * surface_area / surface_area_
    result = definition.rhs.subs(radiant_flux(surface_area), flux_function).doit()
    return Quantity(result)
