from sympy import (Eq, solve)
from sympy.physics.units import electric_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## An infinite plane has an arbitrary surface charge density. The tension lines are perpendicular
## to the plane and directed in both directions from the plane. Then the field strength will depend
## only on the surface dense charge.

## Law is: E = g / (2 * e0), where
## E - electric field intensity,
## g - surface charge density,
## e - electric constant.

electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

surface_charge_density = Symbol("surface_charge_density", units.charge / units.area)

law = Eq(electric_intensity, surface_charge_density / (2 * electric_constant))


def print_law() -> str:
    return print_expression(law)


@validate_input(surface_charge_density_=surface_charge_density)
@validate_output(electric_intensity)
def calculate_electric_intensity(surface_charge_density_: Quantity) -> Quantity:
    result_expr = solve(law, electric_intensity, dict=True)[0][electric_intensity]
    result_expr = result_expr.subs({
        surface_charge_density: surface_charge_density_,
    })
    return Quantity(result_expr)
