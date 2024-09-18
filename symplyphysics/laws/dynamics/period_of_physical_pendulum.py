"""
Period of physical pendulum
===========================

A *physical pendulum* is a pendulum with an arbitrary distribution of mass that
oscillates about a given pivot point. The period of its oscillations depends on its rotational inertia, 
mass and the distance between the pivot and the center of mass of the pendulum.
"""

from sympy import (
    Eq,
    pi,
    sqrt,
    solve,
    symbols as SymSymbols,
    Symbol as SymSymbol,
    Function as SymFunction,
    Derivative,
)
from sympy.physics.units import acceleration_due_to_gravity
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,
    angular_speed_is_angular_distance_derivative as angular_velocity_def,
    harmonic_oscillator_is_second_derivative_equation as oscillator_eqn,
    period_from_angular_frequency as period_law,
)
from symplyphysics.laws.dynamics import (
    acceleration_is_force_over_mass as newtons_second_law,
    moment_of_force_from_moment_of_inertia_and_angular_acceleration as torque_law,
    torque_via_force_and_radial_distance as torque_def,
)

period = Symbol("period", units.time, positive=True)
"""
The period of the physical pendulum.

Symbol:
    :code:`T`
"""

mass = clone_as_symbol(symbols.mass, positive=True)
"""
The :symbols:`mass` of the pendulum.
"""

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2, positive=True)
"""
The rotational inertia of the pendulum.

Symbol:
    :code:`I`
"""

distance_to_pivot = Symbol("distance_to_pivot", units.length, positive=True)
"""
The distance between the pivot and the pendulum's center of mass.

Symbol:
    :code:`h`
"""

law = Eq(
    period,
    2 * pi * sqrt(rotational_inertia / (mass * acceleration_due_to_gravity * distance_to_pivot)))
r"""
:code:`T = 2 * pi * sqrt(I / (m * g * h))`

Latex:
    .. math::
        T = 2 \pi \sqrt{\frac{I}{m g h}}
"""

# Derive from torque definition

time = SymSymbol("time")
# use symbols to avoid "not callable" warning
angle_function = SymSymbols("angle_function", cls=SymFunction)
torque = SymSymbol("torque")

gravitational_force = solve(newtons_second_law.law, newtons_second_law.force)[0].subs({
    newtons_second_law.mass: mass,
    newtons_second_law.acceleration: acceleration_due_to_gravity,
})

angular_velocity = (angular_velocity_def.definition.rhs.subs(angular_velocity_def.time,
    time).subs(angular_velocity_def.angular_distance(time), angle_function(time)))

angular_acceleration = (angular_acceleration_def.definition.rhs.subs(angular_acceleration_def.time,
    time).subs(angular_acceleration_def.angular_speed(time), angular_velocity))

# The factor of -1 indicates that gravitational force acts to reduce the angle.
torque_due_to_gravity = torque_def.law.rhs.subs({
    torque_def.force: -1 * gravitational_force,
    torque_def.radial_distance: distance_to_pivot,
    torque_def.angle_between_vectors: angle_function(time),
})

# Temporary replace function with symbol to make "series" work
angle_sym = SymSymbol("angle_sym")
torque_due_to_gravity = (torque_due_to_gravity.subs(angle_function(time),
    angle_sym).series(angle_sym, 0, 2).removeO().subs(angle_sym, angle_function(time)))

torque_due_to_acceleration = torque_law.law.rhs.subs({
    torque_law.rotational_inertia: rotational_inertia,
    torque_law.angular_acceleration: angular_acceleration,
})

diff_eqn_derived = Eq(torque_due_to_gravity, torque_due_to_acceleration)

diff_eqn_original = (oscillator_eqn.definition.subs(oscillator_eqn.time,
    time).subs(oscillator_eqn.displacement(time), angle_function(time)))

angular_velocity_expr = solve([diff_eqn_derived, diff_eqn_original],
    (Derivative(angle_function(time), (time, 2)), oscillator_eqn.angular_frequency),
    dict=True)[1][oscillator_eqn.angular_frequency]

period_derived = period_law.law.rhs.subs(period_law.angular_frequency, angular_velocity_expr)

assert expr_equals(law.rhs, period_derived)


@validate_input(
    rotational_inertia_=rotational_inertia,
    pendulum_mass_=mass,
    distance_to_pivot_=distance_to_pivot,
)
@validate_output(period)
def calculate_period(
    rotational_inertia_: Quantity,
    pendulum_mass_: Quantity,
    distance_to_pivot_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        mass: pendulum_mass_,
        distance_to_pivot: distance_to_pivot_,
    })
    return Quantity(result)
