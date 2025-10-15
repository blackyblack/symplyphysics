"""
Inner pressure is constant
==========================

Bernoulli's equation applied to an ideal liquid specifies that the inner pressure of the
fluid is constant at all points along a streamline.

**Conditions:**

#. The fluid must be :ref:`ideal <ideal_fluid_def>`.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Bernoulli%27s_principle#Simplified_form>`__.
"""

from sympy import Eq, dsolve, Derivative, solve
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_function,
    clone_as_symbol, quantities)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.solvers import apply

from symplyphysics.definitions import density_from_mass_volume as _density_def
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as _distance_law
from symplyphysics.laws.hydro import (
    volume_flux_is_constant as _continuity_law,
    inner_pressure_is_sum_of_pressures as _inner_pressure_def,
    dynamic_pressure_via_density_and_flow_speed as _dynamic_pressure_def,
    hydrostatic_pressure_via_density_and_height as _hydrostatic_pressure_law,
)
from symplyphysics.laws.dynamics import (
    pressure_from_force_and_area as _pressure_def,
    mechanical_work_from_force_and_distance as _work_law,
    acceleration_is_force_over_mass as _newtons_law,
    kinetic_energy_from_mass_and_speed as _kinetic_energy_law,
    total_work_is_change_in_kinetic_energy as _work_energy_law,
)

time = symbols.time
"""
:symbols:`time`.
"""

inner_pressure = clone_as_function(symbols.pressure, [time],
    display_symbol="p_inner",
    display_latex="p_\\text{inner}")
"""
Inner pressure of the fluid at a chosen point in space as a function of :attr:`~time`.
See :ref:`Inner pressure is sum of pressures`.
"""

law = Eq(Derivative(inner_pressure(time), time), 0)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law. Also see [this](https://en.wikipedia.org/wiki/Bernoulli%27s_principle#Derivations)
# Figure can be found [here](https://en.wikipedia.org/wiki/Bernoulli%27s_principle#/media/File:BernoullisLawDerivationDiagram.svg)

# The system consists of a fluid volume between the cross sections `A_1 = A(0)` and `A_2 = A(t)`.
# We let the fluid elements to flow within a time interval `dt`.

_tube_area = _continuity_law.tube_area
_flow_speed = _continuity_law.flow_speed

_density = _density_def.density
_time_change = clone_as_symbol(time, display_symbol="dt")
_mass_change = clone_as_function(symbols.mass, display_symbol="dm")

# 1. Due to mass conservation the masses of the displaced fluid are equal.

_continuity_dsolved_eqn = dsolve(
    _continuity_law.law.subs(_continuity_law.time, time),
    _tube_area(time),
    ics={_tube_area(0): _tube_area(0)},
)

_continuity_dsolved_multiplied_eqn = apply(
    _continuity_dsolved_eqn,
    lambda x: x * _flow_speed(time) * _density * _time_change,
)

_distance_change_expr = _distance_law.law.rhs.subs({
    _distance_law.initial_position: 0,
    _distance_law.speed: _flow_speed(time),
    _distance_law.time: _time_change,
})

_mass_change_original_expr = solve(
    _density_def.definition,
    _density_def.mass,
)[0].subs({
    _density_def.volume: _tube_area(time) * _distance_change_expr,
})

_mass_change_original_eqn = Eq(_mass_change(time), _mass_change_original_expr)

_initial_mass_change_eqn = _mass_change_original_eqn.subs(time, 0)
_final_mass_change_eqn = _mass_change_original_eqn  # at time `t`

_mass_change_expr = solve(
    (_continuity_dsolved_multiplied_eqn, _initial_mass_change_eqn, _final_mass_change_eqn),
    (_mass_change(time), _flow_speed(time), _flow_speed(0)),
    dict=True,
)[0][_mass_change(time)]

assert expr_equals(_mass_change_expr, _mass_change(0))  # mass change is constant!

_mass_change_eqn = Eq(_mass_change(time), _mass_change_expr)

# 2. Work of the pressure force is the work related to the pressure pushing the fluid upstream or
#    downstream. Note that the pressure work is positive when the fluid is pushed into the tube
#    (i.e. upstream), and negative when the fluid is pushed out of it (i.e. downstream).

_pressure = clone_as_function(symbols.pressure)

_pressure_force_expr = solve(_pressure_def.law, _pressure_def.force)[0].subs({
    _pressure_def.pressure: _pressure(time),
    _pressure_def.area: _tube_area(time),
})

_pressure_work_expr = _work_law.law.rhs.subs({
    _work_law.force: _pressure_force_expr,
    _work_law.distance: _distance_change_expr,
})

_inlet_pressure_work_expr = _pressure_work_expr.subs(time, 0)  # positive work, upstream

_outlet_pressure_work_expr = -1 * _pressure_work_expr  # negative work, downstream

_net_pressure_work_expr = _inlet_pressure_work_expr + _outlet_pressure_work_expr

_net_pressure_work_expr = solve(
    (
    Eq(_work_law.work, _net_pressure_work_expr),
    _initial_mass_change_eqn,
    _final_mass_change_eqn,
    _mass_change_eqn,
    ),
    (_work_law.work, _tube_area(0), _tube_area(time), _mass_change(time)),
    dict=True,
)[0][_work_law.work]

# 3. Work of the gravity force

_height = clone_as_function(symbols.height)

_gravity_force_expr = solve(
    _newtons_law.law,
    _newtons_law.force,
)[0].subs({
    _newtons_law.acceleration: quantities.acceleration_due_to_gravity,
    _newtons_law.mass: _mass_change_expr
})

_gravity_work_expr = _work_law.law.rhs.subs({
    _work_law.force: -1 * _gravity_force_expr,  # projection to the z-axis
    _work_law.distance: _height(time) - _height(0),  # final minus initial
})

# 4. Net work is the sum of pressure work and gravity work.

_net_work = _net_pressure_work_expr + _gravity_work_expr

# 5. Kinetic energy change.

_kinetic_energy_expr = _kinetic_energy_law.law.rhs.subs({
    _kinetic_energy_law.mass: _mass_change_expr,
    _kinetic_energy_law.speed: _flow_speed(time),
})

_initial_kinetic_energy_expr = _kinetic_energy_expr.subs(time, 0)
_final_kinetic_energy_expr = _kinetic_energy_expr  # at time `t`

# 6. Make use of the work—kinetic energy principle.

_work_energy_equality = _work_energy_law.law.subs({
    _work_energy_law.work: _net_work,
    _work_energy_law.kinetic_energy(_work_energy_law.time_before): _initial_kinetic_energy_expr,
    _work_energy_law.kinetic_energy(_work_energy_law.time_after): _final_kinetic_energy_expr,
})

# 7. Recall the definition of inner pressure.

_dynamic_pressure_expr = _dynamic_pressure_def.law.rhs.subs({
    _dynamic_pressure_def.density: _density,
    _dynamic_pressure_def.flow_speed: _flow_speed(time),
})

_hydrostatic_pressure_expr = _hydrostatic_pressure_law.law.rhs.subs({
    _hydrostatic_pressure_law.density: _density,
    _hydrostatic_pressure_law.height: _height(time),
})

_inner_pressure_original_expr = _inner_pressure_def.law.rhs.subs({
    _inner_pressure_def.static_pressure: _pressure(time),
    _inner_pressure_def.dynamic_pressure: _dynamic_pressure_expr,
    _inner_pressure_def.hydrostatic_pressure: _hydrostatic_pressure_expr,
})

_inner_pressure_original_eqn = Eq(inner_pressure(time), _inner_pressure_original_expr)

_initial_inner_pressure_eqn = _inner_pressure_original_eqn.subs(time, 0)
_final_inner_pressure_eqn = _inner_pressure_original_eqn  # at time `t`

_inner_pressure_expr = solve(
    (_work_energy_equality, _initial_inner_pressure_eqn, _final_inner_pressure_eqn),
    (inner_pressure(time), _height(time), _height(0)),
    dict=True,
)[0][inner_pressure(time)]

assert expr_equals(_inner_pressure_expr, inner_pressure(0))

_inner_pressure_eqn = Eq(inner_pressure(time), _inner_pressure_expr)

_inner_pressure_diff_time_eqn = apply(_inner_pressure_eqn, lambda s: s.diff(time))

assert expr_equals(law.lhs, _inner_pressure_diff_time_eqn.lhs)
assert expr_equals(law.rhs, _inner_pressure_diff_time_eqn.rhs)


@validate_input(inner_pressure_before_=inner_pressure)
@validate_output(inner_pressure)
def calculate_inner_pressure(inner_pressure_before_: Quantity) -> Quantity:
    dsolved = dsolve(law, inner_pressure(time))
    result_expr = dsolved.subs("C1", inner_pressure_before_).rhs
    return Quantity(result_expr)


# UNIQUE_LAW_ID: 729
