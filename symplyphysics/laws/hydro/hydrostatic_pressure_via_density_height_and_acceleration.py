"""
Hydrostatic pressure via density, height and acceleration
=========================================================

Hydrostatic pressure is the pressure exerted by a fluid at equilibrium at a given point
within the fluid, due to the force of gravity.

**Conditions:**

#. The only force acting on the fluid is the gravitational force.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import symbols, Quantity, validate_input, validate_output
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import density_from_mass_volume as _density_def
from symplyphysics.laws.dynamics import (acceleration_is_force_over_mass as _newtons_second_law,
    pressure_from_force_and_area as _pressure_def)

density = symbols.density
"""
:symbols:`density` of the fluid.
"""

height = symbols.height
"""
:symbols:`height` of the fluid column.
"""

acceleration = symbols.acceleration
"""
:symbols:`acceleration` of the vessel.
"""

hydrostatic_pressure = symbols.hydrostatic_pressure
"""
:symbols:`hydrostatic_pressure` of the fluid.
"""

law = Eq(hydrostatic_pressure, density * acceleration * height)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

# For a fluid in equilibrium, the pressure inside it does not depend on the direction of
# observation. Therefore we can choose a fluid column of the given height and area oriented along
# the z-axis and calculate the pressure on the bottom side of the column due to the gravity pull of
# the column.

_pressure_eqn = _pressure_def.law.subs(_pressure_def.pressure, hydrostatic_pressure)

_force_eqn = _newtons_second_law.law.subs(_newtons_second_law.acceleration, acceleration)

_density_eqn = _density_def.definition.subs({
    _density_def.density: density,
    _density_def.mass: _newtons_second_law.mass,
    _density_def.volume: _pressure_def.area * height,
})

_pressure_derived = solve(
    (_pressure_eqn, _force_eqn, _density_eqn),
    (hydrostatic_pressure, _newtons_second_law.force, _newtons_second_law.mass),
    dict=True,
)[0][hydrostatic_pressure]

assert expr_equals(_pressure_derived, law.rhs)


@validate_input(density_=density, depth_=height, acceleration_=acceleration)
@validate_output(hydrostatic_pressure)
def calculate_hydrostatic_pressure(density_: Quantity, depth_: Quantity,
    acceleration_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, hydrostatic_pressure, dict=True)[0][hydrostatic_pressure]
    result_expr = result_pressure_expr.subs({
        density: density_,
        height: depth_,
        acceleration: acceleration_,
    })
    return Quantity(result_expr)
