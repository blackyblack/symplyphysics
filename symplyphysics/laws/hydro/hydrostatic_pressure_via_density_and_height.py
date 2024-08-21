"""
Hydrostatic pressure via density and height
===========================================

Hydrostatic pressure is the pressure exerted by a fluid at equilibrium at any point of time
due to the force of gravity.

**Conditions:**

#. The fluid is in statical equilibrium.
#. The only force acting on the fluid is the force of gravity.
"""

from sympy import Eq, solve
from symplyphysics import units, Quantity, Symbol, validate_input, validate_output
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.hydro import hydrostatic_pressure_from_density_and_depth_acceleration as pressure_law

hydrostatic_pressure = Symbol("hydrostatic_pressure", units.pressure)
"""
Hydrostatic pressure of the fluid column.

Symbol:
    :code:`p`
"""

density = Symbol("density", units.mass / units.volume)
r"""
Density of the fluid.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

height = Symbol("height", units.length)
r"""
Height of the fluid column.

Symbol:
    :code:`h`
"""

law = Eq(hydrostatic_pressure, density * units.acceleration_due_to_gravity * height)
r"""
:code:`p = rho * g * h`

Latex:
    .. math::
        p = \rho g h
"""

# This law might be derived via hydrostatic pressure law.

_pressure_law_applied = pressure_law.law.subs({
    pressure_law.density: density,
    pressure_law.depth: height,
    pressure_law.acceleration: units.acceleration_due_to_gravity,
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
