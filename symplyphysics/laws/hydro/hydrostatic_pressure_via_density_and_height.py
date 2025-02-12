"""
Hydrostatic pressure via density and height
===========================================

Hydrostatic pressure is the pressure exerted by a fluid at equilibrium at any point of time
due to the force of gravity.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Conditions:**

#. The fluid is in statical equilibrium.
#. The only force acting on the fluid is the force of gravity.
#. The fluid is subjected to the gravity force of Earth.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Pressure#Liquid_pressure>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import hydrostatic_pressure_from_density_and_depth_acceleration as pressure_law

hydrostatic_pressure = symbols.hydrostatic_pressure
"""
:symbols:`hydrostatic_pressure` of the fluid column.
"""

density = symbols.density
"""
:symbols:`density` of the fluid.
"""

height = symbols.height
r"""
:symbols:`height` of the fluid column.
"""

law = Eq(hydrostatic_pressure, density * quantities.acceleration_due_to_gravity * height)
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via hydrostatic pressure law.

_pressure_law_applied = pressure_law.law.subs({
    pressure_law.density: density,
    pressure_law.height: height,
    pressure_law.acceleration: quantities.acceleration_due_to_gravity,
})
_pressure_derived = solve(_pressure_law_applied, pressure_law.hydrostatic_pressure,
    dict=True)[0][pressure_law.hydrostatic_pressure]

# Check if derived pressure is same as declared.
assert expr_equals(_pressure_derived, law.rhs)


@validate_input(density_=density, depth_=height)
@validate_output(hydrostatic_pressure)
def calculate_hydrostatic_pressure(density_: Quantity, depth_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, hydrostatic_pressure, dict=True)[0][hydrostatic_pressure]
    result_expr = result_pressure_expr.subs({density: density_, height: depth_})
    return Quantity(result_expr)
