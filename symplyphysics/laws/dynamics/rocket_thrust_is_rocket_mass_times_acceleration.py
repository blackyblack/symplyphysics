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

#. The quantity :math:`R v_\text{rel}` is called the *thrust of rocket engine*.

    .. _rate note:

#. The rate :math:`R` of fuel consumption is defined as 
    .. math::
        R = - \frac{d M}{d t}

    where :math:`M` is the :ref:`rocket mass <rocket mass definition>`.
"""

from sympy import Eq, dsolve, solve, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
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

fuel_consumption_rate = Symbol("fuel_consumption_rate", units.mass / units.time)
r"""
The rate of fuel consumption. See :ref:`Note <rate note>` for the definition.

Symbol:
    :code:`R`
"""

relative_speed = Symbol("relative_speed", units.velocity)
r"""
The speed of the rocket relative to its products.

Symbol:
    :code:`v_rel`

Latex:
    :math:`v_\text{rel}`
"""

acceleration = symbols.kinematics.acceleration
"""
The :attr:`~symplyphysics.symbols.kinematics.acceleration` of the rocket.

Symbol:
    :code:`a`
"""

mass = symbols.basic.mass
"""
.. _rocket mass definition:

The :attr:`~symplyphysics.symbols.basic.mass` of the rocket

Symbol:
    :code:`M`
"""

law = Eq(fuel_consumption_rate * relative_speed, mass * acceleration)
r"""
:code:`R * v_rel = M * a`

Latex:
    :math:`R v_\text{rel} = M a`
"""

# Derive this law from the law of conservation of momentum.

rocket_speed = SymSymbol("rocket_speed")
rocket_speed_change = SymSymbol("rocket_speed_change")
fuel_mass_thrusted = SymSymbol("fuel_mass_thrusted")

rocket_momentum_before_release = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass,
    momentum_def.speed: rocket_speed,
})

rocket_momentum_after_release = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass - fuel_mass_thrusted,
    momentum_def.speed: rocket_speed + rocket_speed_change,
})

products_speed = SymSymbol("products_speed")

products_momentum = momentum_def.definition.rhs.subs({
    momentum_def.mass: fuel_mass_thrusted,
    momentum_def.speed: products_speed,
})

# Summary momentum is conserved
final_momentum = rocket_momentum_after_release + products_momentum

rocket_speed_relative_to_frame = rocket_speed + rocket_speed_change
rocket_speed_relative_to_products = relative_speed
products_speed_relative_to_frame = products_speed

relative_speed_eqn = galilean_law.law.subs({
    galilean_law.body_velocity_in_first_frame: rocket_speed_relative_to_frame,
    galilean_law.body_velocity_in_second_frame: rocket_speed_relative_to_products,
    galilean_law.second_frame_velocity_in_first_frame: products_speed_relative_to_frame,
})

momentum_conservation_eqn = momentum_conservation_law.law.subs({
    momentum_conservation_law.momentum(momentum_conservation_law.time_before):
        rocket_momentum_before_release,
    momentum_conservation_law.momentum(momentum_conservation_law.time_after):
        final_momentum,
})

relative_velocity_expr = solve(
    [relative_speed_eqn, momentum_conservation_eqn],
    (products_speed, relative_speed),
    dict=True,
)[0][relative_speed]

time_change = SymSymbol("time_change")

# solve differential equation with constant fuel_consumption_rate
dsolved_fuel_mass = dsolve(
    flow_rate_def.definition.subs(flow_rate_def.mass_flow_rate(flow_rate_def.time),
    fuel_consumption_rate), flow_rate_def.mass(flow_rate_def.time))
fuel_consumption_eqn = Eq(fuel_mass_thrusted, dsolved_fuel_mass.rhs)
# C1 is initial fuel mass thrusted
fuel_consumption_eqn = fuel_consumption_eqn.subs({"C1": 0, flow_rate_def.time: time_change})

rocket_acceleration_expr = solve(
    acceleration_def.law,
    acceleration_def.acceleration,
)[0].subs({
    acceleration_def.final_speed: rocket_speed_change,
    acceleration_def.time: time_change,
    acceleration_def.initial_speed: 0,
})

relative_velocity_eqn_system = [
    Eq(relative_speed, relative_velocity_expr),
    fuel_consumption_eqn,
    Eq(acceleration, rocket_acceleration_expr),
]

relative_velocity_derived = solve(relative_velocity_eqn_system,
    (relative_speed, fuel_mass_thrusted, rocket_speed_change),
    dict=True)[0][relative_speed]

relative_velocity_from_law = solve(law, relative_speed)[0]

assert expr_equals(relative_velocity_derived, relative_velocity_from_law)


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
