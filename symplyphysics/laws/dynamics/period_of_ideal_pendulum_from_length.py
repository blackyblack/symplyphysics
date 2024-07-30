"""
Period of ideal pendulum from length
====================================

An ideal pendulum is an object hanging on a thread. In a gravitational field it starts oscillating after being pushed out of balance.
Period of pendulum oscillation does not depend on its mass.

**Conditions:**

#. The angle between pendulum and gravity vector is fairly small (less than 15 degrees).
#. Ideal pendulum doesn't gain or lose any energy, so there is no friction in the system.
#. The object is small enough to be considered a material point.
#. The thread is weightless and doesn't change its length.
"""

from sympy import (Derivative, Eq, Function as SymFunction, diff, sin, solve, pi, sqrt, symbols)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.geometry import planar_projection_is_cosine as projector
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as angular_velocity_law
from symplyphysics.definitions import period_from_angular_frequency as angular_frequency
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as oscillator
from symplyphysics.definitions import mechanical_energy_is_kinetic_and_potential as mechanical_energy_def
from symplyphysics.laws.conservation import mechanical_energy_is_constant as mechanical_energy_conservation

period = Symbol("period", units.time)
"""
The period of oscillations.

Symbol:
    :code:`T`
"""

length = Symbol("length", units.length)
"""
The length of the pendulum.

Symbol:
    :code:`L`
"""

law = Eq(period, 2 * pi * sqrt(length / units.acceleration_due_to_gravity))
r"""
:code:`T = 2 * pi * sqrt(L / g)`

Latex:
    .. math::
        T = 2 \pi \sqrt{\frac{L}{g}}
"""

# Derive this law from conservation of energy
## Polar coordinate system is selected for this task. Center is a fixed point of the thread.

pendulum_mass = symbols("pendulum_mass")
## Pendulum angle is angle between thread and gravity vector. In balanced position it is 0.
pendulum_angle = symbols("pendulum_angle", cls=SymFunction)
time = symbols("time")

## Pendulum oscillation is cyclic transfer of energy from kinetic to potential. To set oscillation up we have to input some energy. Usually it is done by biasing the pendulum to some angle and letting it go.
## Biasing the pendulum is giving to it some amount of potential energy.

pendulum_height_before = length
pendulum_height_after = projector.law.subs({
    projector.vector_length: length,
    projector.vector_angle: pendulum_angle(time)
}).rhs
amount_of_potential_energy = potential_energy.law.subs({
    potential_energy.mass: pendulum_mass,
    potential_energy.height: (pendulum_height_before - pendulum_height_after)
}).rhs

## Kinetic energy of the pendulum is:
## pendulum_mass * (length * angular_velocity)**2 / 2

linear_velocity = angular_velocity_law.law.subs({
    angular_velocity_law.curve_radius: length,
    angular_velocity_law.angular_velocity: Derivative(pendulum_angle(time), time),
}).rhs
amount_of_kinetic_energy = kinetic_energy.law.subs({
    kinetic_energy.mass: pendulum_mass,
    kinetic_energy.speed: linear_velocity
}).rhs

mechanical_energy = mechanical_energy_def.definition.subs({
    mechanical_energy_def.kinetic_energy: amount_of_kinetic_energy,
    mechanical_energy_def.potential_energy: amount_of_potential_energy
}).rhs

## Total mechanical energy for pendulum is constant

conserved_energy_eq = mechanical_energy_conservation.law.subs(mechanical_energy_conservation.time,
    time)

## Differentiate both sides of equation.
## Derivative of constant mechanical energy will be zero, so is the left side of this equation.
total_energy_diff_eq = Eq(Derivative(mechanical_energy_conservation.mechanical_energy(time), time),
    diff(mechanical_energy, time))

## We do not replace it with zero, but solve system of equations instead
total_energy_diff_solved = solve([total_energy_diff_eq, conserved_energy_eq],
    (Derivative(pendulum_angle(time),
    (time, 2)), Derivative(mechanical_energy_conservation.mechanical_energy(time), time)),
    dict=True)[0][Derivative(pendulum_angle(time), (time, 2))]
## Now we've found the solution for second order derivative of angle function over time
total_energy_diff_solved_eq = Eq(Derivative(pendulum_angle(time), (time, 2)),
    total_energy_diff_solved)

#NOTE: large displacement angle (over 15 degrees) gives quite a complex solution for the differential equation.

# For small angles, sin(pendulum_angle) can be reduced to pendulum_angle
small_angle_harmonic_oscillation_eq = total_energy_diff_solved_eq.subs(
    sin(pendulum_angle(time)), pendulum_angle(time))

# Will result in harmonic oscillator equation:
## Derivative(pendulum_angle(time), (time, 2)) = -free_fall_acceleration / length * pendulum_angle(time)
oscillator_eq = oscillator.definition.subs(oscillator.time, time)
oscillator_eq = oscillator_eq.subs(oscillator.displacement(time), pendulum_angle(time))
angular_frequency_solved = solve([oscillator_eq, small_angle_harmonic_oscillation_eq],
    (oscillator.angular_frequency, pendulum_angle(time)),
    dict=True)[0][oscillator.angular_frequency]

## Check that expected period matches our law.
## Square roots fail to compare with each other. Raise both parts to power of 2 before checking for equality.
oscillation_period_derived = angular_frequency.law.subs(angular_frequency.angular_frequency,
    angular_frequency_solved).rhs
assert expr_equals(oscillation_period_derived**2, law.rhs**2)


@validate_input(pendulum_length_=length)
@validate_output(period)
def calculate_period(pendulum_length_: Quantity) -> Quantity:
    solved = solve(law, period, dict=True)[0][period]
    result_expr = solved.subs(length, pendulum_length_)
    return Quantity(result_expr)
