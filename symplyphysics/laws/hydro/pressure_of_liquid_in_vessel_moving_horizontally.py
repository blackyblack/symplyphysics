"""
Pressure of liquid in vessel moving horizontally
================================================

If a vessel with a liquid moves with horizontal acceleration, then the hydrostatic
pressure of the liquid depends on the density of the liquid, the acceleration of free
fall, the horizontal acceleration of the vessel and the height of liquid.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

..
    TODO: find link
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    Vector,
    add_cartesian_vectors,
    symbols,
    Quantity,
    validate_input,
    validate_output,
    vector_magnitude,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import hydrostatic_pressure_from_density_and_depth_acceleration as pressure_law

pressure = symbols.pressure
"""
:symbols:`pressure` of the liquid.
"""

density = symbols.density
"""
:symbols:`density` of the liquid.
"""

acceleration = symbols.acceleration
"""
:symbols:`acceleration` of the vessel.
"""

height = symbols.height
"""
:symbols:`height` of the liquid column.
"""

law = Eq(pressure,
    density * sqrt(quantities.acceleration_due_to_gravity**2 + acceleration**2) * height)
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via hydrostatic pressure law.
# The vessel moves horizontally and the pressure exerted by the resultant force on a surface of equal pressure,
# which is inclined with respect to the horizon, is considered. The modulus of the optimizing force will be equal
# to m * sqrt(g^2+ a^2).

# Earth free fall acceleration is directed downwards. It is equivalent to vessel acceleration
# directed vertically upwards.
_free_fall_acceleration_vector = Vector([0, quantities.acceleration_due_to_gravity])
# Horizontal vector
_vessel_acceleration_vector = Vector([acceleration, 0])
_total_acceleration = vector_magnitude(
    add_cartesian_vectors(_free_fall_acceleration_vector, _vessel_acceleration_vector))

_pressure_law_applied = pressure_law.law.subs({
    pressure_law.density: density,
    pressure_law.height: height,
    pressure_law.acceleration: _total_acceleration,
})
_pressure_derived = solve(_pressure_law_applied, pressure_law.hydrostatic_pressure,
    dict=True)[0][pressure_law.hydrostatic_pressure]

# Check if derived pressure is same as declared.
assert expr_equals(_pressure_derived, law.rhs)


@validate_input(density_liquid_=density, acceleration_=acceleration, height_=height)
@validate_output(pressure)
def calculate_pressure(density_liquid_: Quantity, acceleration_: Quantity,
    height_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_expr = result_expr.subs({
        density: density_liquid_,
        acceleration: acceleration_,
        height: height_
    })
    return Quantity(result_expr)
