"""
Efflux speed via pressure and density
=====================================

The speed of a fluid flowing out from a small orifice can be expressed as a function
of the fluid's :doc:`hydrostatic pressure <laws.hydro.hydrostatic_pressure_via_density_and_height>`
and density.

**Conditions:**

#. The orifice is very small relative to the horizontal cross-section of the container.
#. The fluid is :ref:`ideal <ideal_fluid_def>`.

..
    TODO find link
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import hydrostatic_pressure_via_density_and_height as pressure_law
from symplyphysics.laws.hydro import efflux_speed_via_height as velocity_law

efflux_speed = Symbol("efflux_speed", units.velocity)
"""
Speed of the fluid flowing out of the pipe.

Symbol:
    :code:`v`
"""

hydrostatic_pressure = Symbol("hydrostatic_pressure", units.pressure)
"""
Hydrostatic pressure of the fluid.

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

law = Eq(efflux_speed, sqrt(2 * hydrostatic_pressure / density))
r"""
:code:`v = sqrt(2 * p / rho)`

Latex:
    .. math::
        v = \sqrt{\frac{2 p}{\rho}}
"""

# This law might be derived via "hydrostatic_pressure_from_density_and_depth" law
# and "velocity_from_height" law.

_pressure_law_applied = pressure_law.law.subs({
    pressure_law.density: density,
    pressure_law.hydrostatic_pressure: hydrostatic_pressure
})
_height_derived = solve(_pressure_law_applied, pressure_law.height,
    dict=True)[0][pressure_law.height]

_velocity_law_applied = velocity_law.law.subs({velocity_law.height: _height_derived})
_velocity_derived = solve(_velocity_law_applied, velocity_law.efflux_speed,
    dict=True)[0][velocity_law.efflux_speed]

# Check if derived efflux_speed is same as declared.
assert expr_equals(_velocity_derived, law.rhs)


@validate_input(pressure_=hydrostatic_pressure, density_=density)
@validate_output(efflux_speed)
def calculate_velocity(pressure_: Quantity, density_: Quantity) -> Quantity:
    result_expr = solve(law, efflux_speed, dict=True)[0][efflux_speed]
    result_expr = result_expr.subs({
        hydrostatic_pressure: pressure_,
        density: density_,
    })
    return Quantity(result_expr)
