r"""
Rocket thrust is rocket mass times acceleration
===============================================

Assuming we are at rest relative to an inertial reference frame, we observe a rocket
through space with no gravitational or atmospheric drag forces acting on it.
The mass of the rocket changes as it burns fuel and releases the products of burning,
the total mass of the system does not change.

**Conditions:**

#. The fuel consumption rate is constant.
#. The velocities are non-relativistic.

**Notes:**

#. The quantity :math:`R v_\text{rel}` is called the **thrust of rocket engine**.

    .. _rate_note:

#. The rate :math:`R` of fuel consumption is defined as
    .. math::

        R = - \frac{d m}{d t}

    where :math:`m` is the :ref:`rocket mass <rocket_mass_definition>`.

**Links:**

#. Equation 9-87 on p. 242 of "Fundamentals of Physics" by David Halladay et al., 10th Ed.
"""

from sympy import Eq, dsolve, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import (
    momentum_after_collision_equals_to_momentum_before as momentum_conservation_law,)
from symplyphysics.definitions import (momentum_is_mass_times_speed as momentum_def, mass_flow_rate
    as flow_rate_def)
from symplyphysics.laws.kinematics import (
    classical_addition_of_velocities as galilean_law,
    speed_via_constant_acceleration_and_time as acceleration_def,
)

fuel_consumption_rate = clone_as_symbol(symbols.mass_flow_rate,
    display_symbol="R",
    display_latex="R")
"""
The rate of fuel consumption, or :symbols:`mass_flow_rate` of exhaust. See :ref:`Note <rate_note>`
for the definition.
"""

relative_speed = clone_as_symbol(symbols.speed,
    display_symbol="v_rel",
    display_latex="v_\\text{rel}")
"""
The :symbols:`speed` of the rocket relative to its products.
"""

acceleration = symbols.acceleration
"""
The :symbols:`acceleration` of the rocket.
"""

mass = symbols.mass
"""
.. _rocket_mass_definition:

The :symbols:`mass` of the rocket
"""

law = Eq(fuel_consumption_rate * relative_speed, mass * acceleration)
"""
:laws:symbol::

:laws:latex::
"""

# Derive this law from the law of conservation of momentum.

_rocket_speed = clone_as_symbol(symbols.speed)
_rocket_speed_change = clone_as_symbol(symbols.speed, display_symbol="dv")
_fuel_mass_thrusted = clone_as_symbol(symbols.mass, display_symbol="m_thrust")

_rocket_momentum_before_release = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass,
    momentum_def.speed: _rocket_speed,
})

_rocket_momentum_after_release = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass - _fuel_mass_thrusted,
    momentum_def.speed: _rocket_speed + _rocket_speed_change,
})

_products_speed = clone_as_symbol(symbols.speed, display_symbol="v_product")

_products_momentum = momentum_def.definition.rhs.subs({
    momentum_def.mass: _fuel_mass_thrusted,
    momentum_def.speed: _products_speed,
})

# Summary momentum is conserved
_final_momentum = _rocket_momentum_after_release + _products_momentum

_rocket_speed_relative_to_frame = _rocket_speed + _rocket_speed_change
_rocket_speed_relative_to_products = relative_speed
_products_speed_relative_to_frame = _products_speed

_relative_speed_eqn = galilean_law.law.subs({
    galilean_law.body_speed_in_first_frame: _rocket_speed_relative_to_frame,
    galilean_law.body_speed_in_second_frame: _rocket_speed_relative_to_products,
    galilean_law.second_frame_speed_in_first_frame: _products_speed_relative_to_frame,
})

_momentum_conservation_eqn = momentum_conservation_law.law.subs({
    momentum_conservation_law.momentum(momentum_conservation_law.initial_time):
        _rocket_momentum_before_release,
    momentum_conservation_law.momentum(momentum_conservation_law.final_time):
        _final_momentum,
})

_relative_velocity_expr = solve(
    [_relative_speed_eqn, _momentum_conservation_eqn],
    (_products_speed, relative_speed),
    dict=True,
)[0][relative_speed]

_time_change = clone_as_symbol(symbols.time, display_symbol="dt")

# solve differential equation with constant fuel_consumption_rate
_dsolved_fuel_mass = dsolve(
    flow_rate_def.definition.subs(flow_rate_def.mass_flow_rate(flow_rate_def.time),
    fuel_consumption_rate), flow_rate_def.mass(flow_rate_def.time))
_fuel_consumption_eqn = Eq(_fuel_mass_thrusted, _dsolved_fuel_mass.rhs)
# C1 is initial fuel mass thrusted
_fuel_consumption_eqn = _fuel_consumption_eqn.subs({"C1": 0, flow_rate_def.time: _time_change})

_rocket_acceleration_expr = solve(
    acceleration_def.law,
    acceleration_def.acceleration,
)[0].subs({
    acceleration_def.final_speed: _rocket_speed_change,
    acceleration_def.time: _time_change,
    acceleration_def.initial_speed: 0,
})

_relative_velocity_eqn_system = [
    Eq(relative_speed, _relative_velocity_expr),
    _fuel_consumption_eqn,
    Eq(acceleration, _rocket_acceleration_expr),
]

_relative_velocity_derived = solve(_relative_velocity_eqn_system,
    (relative_speed, _fuel_mass_thrusted, _rocket_speed_change),
    dict=True)[0][relative_speed]

_relative_velocity_from_law = solve(law, relative_speed)[0]

assert expr_equals(_relative_velocity_derived, _relative_velocity_from_law)


@validate_input(
    fuel_consumption_rate_=fuel_consumption_rate,
    rocket_mass_=mass,
    rocket_acceleration_=acceleration,
)
@validate_output(relative_speed)
def calculate_relative_velocity(
    fuel_consumption_rate_: Quantity,
    rocket_mass_: Quantity,
    rocket_acceleration_: Quantity,
) -> Quantity:
    result = solve(law, relative_speed)[0].subs({
        fuel_consumption_rate: fuel_consumption_rate_,
        mass: rocket_mass_,
        acceleration: rocket_acceleration_,
    })
    return Quantity(result)
