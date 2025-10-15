"""
Pressure of liquid in vessel moving vertically
==============================================

If a vessel with a liquid moves with vertical acceleration, then the hydrostatic
pressure of the liquid depends on the density of the liquid, the acceleration of free
fall, the vertical acceleration of the vessel and the height of liquid.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import symbols, Quantity, validate_input, validate_output, quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import hydrostatic_pressure_via_density_height_and_acceleration as pressure_law

from symplyphysics.core.experimental.vectors import VectorNorm
from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector

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
    density * sqrt((quantities.acceleration_due_to_gravity + acceleration)**2) * height)
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via hydrostatic pressure law.
# The vessel moves vertically and the pressure exerted by the resultant force on a surface of equal pressure
# is considered. The modulus of the resulting force will be equal to m * sqrt((g+a)^2).
# If the acceleration "a" is positive, then it is directed upwards. If the acceleration "a" is negative,
# then it is directed downward.

_free_fall_acceleration_vector = CoordinateVector(
    [0, quantities.acceleration_due_to_gravity, 0],
    CARTESIAN,
)

# Vertical vector
_vessel_acceleration_vector = CoordinateVector([0, acceleration, 0], CARTESIAN)

_total_acceleration_vector = CoordinateVector.from_expr(
    _free_fall_acceleration_vector + _vessel_acceleration_vector,)
_total_acceleration = VectorNorm(_total_acceleration_vector)

_pressure_law_applied = pressure_law.law.subs({
    pressure_law.density: density,
    pressure_law.height: height,
    pressure_law.acceleration: _total_acceleration,
})

_pressure_derived = solve(
    _pressure_law_applied,
    pressure_law.hydrostatic_pressure,
    dict=True,
)[0][pressure_law.hydrostatic_pressure]

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


# UNIQUE_LAW_ID: 713
