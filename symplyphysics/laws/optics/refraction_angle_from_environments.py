import numbers
from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, sin, pi, SI
)
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Dimensionless, Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## If ray of light comes from one media to another, it refracts.
## Incoming ray, refracted ray and perpendicular to media boundary are within the same plane.

## Law: n1 * sin(alpha) = n2 * sin(beta)
## Where:
## n1 is refractive index of first media.
## n2 is refractive index of second media.
## alfa is incoming angle. Angle is measured between normal vector to a media boundary plane directed outside of refracting media and incoming ray of light.
## beta is refraction angle. Angle is measured between normal vector to a media boundary plane directed inside refracting media and refracted ray of light.

# Conditions
## - alfa is in [-pi/2, pi/2] range, meaning that incoming ray comes from the outer media.
## - beta is in [-pi/2, pi/2] range, meaning that refracted ray stays in the refracting media.
## - light is monochromic, as refactive index depends on the light frequency.
## - refracting media is uniform, so refracting index does not change over ray path.

incedence_refractive_index = Symbol("incedence_refractive_index", Dimensionless)
resulting_refractive_index = Symbol("resulting_refractive_index", Dimensionless)
incedence_angle = Symbol("incedence_angle", angle_type)
refraction_angle = Symbol("refraction_angle", angle_type)

law = Eq(incedence_refractive_index * sin(incedence_angle), resulting_refractive_index * sin(refraction_angle))

def print(expr: Expr) -> str:
    symbols = [incedence_refractive_index, resulting_refractive_index, incedence_angle, refraction_angle]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(incedence_angle_=incedence_angle, incedence_refractive_index_=incedence_refractive_index, resulting_refractive_index_=resulting_refractive_index)
@validate_output_symbol(refraction_angle)
def calculate_refraction_angle(incedence_angle_: Quantity | float, incedence_refractive_index_: float, resulting_refractive_index_: float) -> Quantity:
    #HACK: sympy angles are always in radians
    incedence_angle_radians = incedence_angle_ if isinstance(incedence_angle_, numbers.Number) else incedence_angle_.scale_factor
    # Check for boundary conditions
    assert incedence_angle_radians <= pi/2
    assert incedence_angle_radians >= -pi/2
    solutions = solve(law, refraction_angle, dict=True)
    result_expr = solutions[0][refraction_angle]
    angle_applied = result_expr.subs({
        incedence_angle: incedence_angle_radians,
        incedence_refractive_index: incedence_refractive_index_,
        resulting_refractive_index: resulting_refractive_index_})
    #HACK: there are 2 solutions for refraction_angle: pi - asin() and asin(). We choose former and switch to latter
    #      if resulting angle does not pass boundary conditions.
    if(angle_applied > pi/2 or angle_applied < -pi/2):
        result_expr = solutions[1][refraction_angle]
        angle_applied = result_expr.subs({
            incedence_angle: incedence_angle_radians,
            incedence_refractive_index: incedence_refractive_index_,
            resulting_refractive_index: resulting_refractive_index_})

    # Check for boundary conditions
    assert angle_applied <= pi/2
    assert angle_applied >= -pi/2
    #HACK: angle type is automatically detected as dimensionless. Force it to angle.
    return Quantity(angle_applied * units.radian, dimension=angle_type)
