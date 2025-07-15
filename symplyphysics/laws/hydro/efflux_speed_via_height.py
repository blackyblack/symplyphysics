"""
Efflux speed via height
=======================

The speed of a fluid flowing out from a small orifice can be expressed as a function of
the height of the fluid column. It is also known as the **Torricelli's law**.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Conditions:**

#. The orifice is very small relative to the horizontal cross-section of the container.
#. The fluid is :ref:`ideal <ideal_fluid_def>`.
#. The fluid is subjected to the gravity force of Earth.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Torricelli%27s_law#>`__.
"""

from sympy import Eq, solve, sqrt, dsolve, Symbol as SymSymbol
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.hydro import (
    inner_pressure_is_constant as _bernoulli_eqn,
    inner_pressure_is_sum_of_pressures as _inner_pressure_def,
    hydrostatic_pressure_via_density_and_height as _hydrostatic_law,
    dynamic_pressure_via_density_and_flow_speed as _dynamic_law,
    volume_flux_is_constant as _continuity_law,
)

efflux_speed = symbols.flow_speed
"""
:symbols:`flow_speed` of the fluid flowing out of the pipe.
"""

height = symbols.height
"""
:symbols:`height` of the fluid column above the orifice.
"""

law = Eq(efflux_speed, sqrt(2 * quantities.acceleration_due_to_gravity * height))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the Bernoulli's equation and the constancy of the volume flux
# Also refer to [this](https://en.wikipedia.org/wiki/Torricelli%27s_law#Derivation)

# 1. Bernoulli's equation

# The fluid is incompressible, therefore its density is the same at all points
_fluid_density = clone_as_symbol(symbols.density)

_flow_speed = clone_as_symbol(symbols.flow_speed)
_flow_speed_at_surface = clone_as_symbol(_flow_speed, subscript="1")
_flow_speed_at_orifice = clone_as_symbol(_flow_speed, subscript="2")

# Height is measured above an arbitrary reference point
_height = clone_as_symbol(symbols.height)
_height_at_surface = clone_as_symbol(_height, subscript="1")
_height_at_orifice = clone_as_symbol(_height, subscript="2")

_dynamic_pressure_expr = _dynamic_law.law.rhs.subs({
    _dynamic_law.density: _fluid_density,
    _dynamic_law.flow_speed: _flow_speed,
})

_hydrostatic_pressure_expr = _hydrostatic_law.law.rhs.subs({
    _hydrostatic_law.density: _fluid_density,
    _hydrostatic_law.height: _height,
})

# We assume the static pressure to be equal at all points of the fluid (in Earth's conditions it
# can be set to the atmospheric pressure):
_inner_pressure_expr = _inner_pressure_def.law.rhs.subs({
    _inner_pressure_def.dynamic_pressure: _dynamic_pressure_expr,
    _inner_pressure_def.hydrostatic_pressure: _hydrostatic_pressure_expr,
})

_inner_pressure_at_surface = _inner_pressure_expr.subs({
    _flow_speed: _flow_speed_at_surface,
    _height: _height_at_surface,
})

_inner_pressure_at_orifice = _inner_pressure_expr.subs({
    _flow_speed: _flow_speed_at_orifice,
    _height: _height_at_orifice,
})

_inner_pressure_dsolved_eqn = dsolve(
    _bernoulli_eqn.law,
    _bernoulli_eqn.inner_pressure(_bernoulli_eqn.time),
    ics={
    _bernoulli_eqn.inner_pressure(0): _inner_pressure_at_surface
    },
).subs(
    _bernoulli_eqn.inner_pressure(_bernoulli_eqn.time),
    _inner_pressure_at_orifice,
)

# 2. Continuity equation

_cross_section_at_surface = clone_as_symbol(symbols.area, subscript="1")
_cross_section_at_orifice = clone_as_symbol(symbols.area, subscript="2")

_continuity_dsolved_eqn = dsolve(
    _continuity_law.law,
    _continuity_law.tube_area(_continuity_law.time),
    ics={
    _continuity_law.tube_area(0): _cross_section_at_surface
    },
).subs({
    _continuity_law.flow_speed(0): _flow_speed_at_surface,
    _continuity_law.flow_speed(_continuity_law.time): _flow_speed_at_orifice,
    _continuity_law.tube_area(_continuity_law.time): _cross_section_at_orifice,
})

# 3. Combine with the remaining equations

_cross_section_ratio = SymSymbol("k", positive=True)
_cross_section_eqn = Eq(
    _cross_section_ratio,
    _cross_section_at_orifice / _cross_section_at_surface,
)

_height_eqn = Eq(
    height,  # thickness of the fluid layer above the orifice, per definition
    _height_at_surface - _height_at_orifice,
)

_efflux_speed_derived = solve(
    (_inner_pressure_dsolved_eqn, _continuity_dsolved_eqn, _cross_section_eqn, _height_eqn),
    (_flow_speed_at_orifice, _flow_speed_at_surface, _cross_section_at_surface, _height_at_surface),
    dict=True,
)[0][_flow_speed_at_orifice]

# Refer to the first condition of this law
_efflux_speed_derived = _efflux_speed_derived.limit(_cross_section_ratio, 0)

assert expr_equals(_efflux_speed_derived, law.rhs)


@validate_input(height_=height)
@validate_output(efflux_speed)
def calculate_velocity(height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, efflux_speed, dict=True)[0][efflux_speed]
    result_expr = result_velocity_expr.subs({height: height_})
    return Quantity(result_expr)
