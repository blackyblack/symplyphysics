"""
Period of spring from mass
==========================

Mass on spring is a system of object with mass :math:`m` and spring with stiffness :math:`k`.
It starts oscillating after being pushed out of balance.

**Conditions:**

#. The spring does not gain or lose energy. There is no friction or inelastic deformation.
#. The object is small enough to be considered a material point.
#. The spring is weightless.
"""

from sympy import (Derivative, Eq, Function as SymFunction, diff, solve, pi, sqrt, symbols as
    SymSymbols, simplify)
from symplyphysics import (Quantity, units, Symbol, validate_input, validate_output, symbols)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import potential_energy_from_deformation as spring_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy
from symplyphysics.definitions import speed_is_distance_derivative as velocity_def
from symplyphysics.definitions import period_from_angular_frequency as period_definition
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as oscillator

period = Symbol("period", units.time)
"""
The period of spring oscillations.

Symbol:
    :code:`T`
"""

mass = symbols.basic.mass
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the object attached to the spring.

Symbol:
    :code:`m`
"""

stiffness = Symbol("stiffness", units.force / units.length)
"""
Spring's stiffness, or spring constant.

Symbol:
    :code:`k`
"""

law = Eq(period, 2 * pi * sqrt(mass / stiffness))
r"""
:code:`T = 2 * pi * sqrt(m / k)`

Latex:
    .. math::
        T = 2 \pi \sqrt{\frac{m}{k}}
"""

# Derive this law from conservation of energy

## Spring displacement is distance between current and balanced positions.
spring_displacement = SymSymbols("spring_displacement", cls=SymFunction)
time = SymSymbols("time")

## Spring oscillation is cyclic transfer of energy from kinetic to potential. To set oscillation up we have to input some energy. Usually it is done by biasing the spring and letting it go.
## Biasing the spring is giving to it some amount of potential energy.

amount_of_potential_energy = spring_energy.law.subs({
    spring_energy.stiffness: stiffness,
    spring_energy.displacement: spring_displacement(time)
}).rhs

## Kinetic energy of the pendulum is:
## object_mass * (linear_velocity)**2 / 2
velocity_def_eq = velocity_def.definition.subs(velocity_def.time, time)
linear_velocity = velocity_def_eq.subs(velocity_def.distance(time), spring_displacement(time)).rhs
amount_of_kinetic_energy = kinetic_energy.law.subs({
    kinetic_energy.mass: mass,
    kinetic_energy.speed: linear_velocity
}).rhs

## Total energy is constant and any of it's derivatives is 0.
total_energy = SymSymbols("total_energy", constant=True)
total_energy_eq = Eq(total_energy, amount_of_kinetic_energy + amount_of_potential_energy)

## Differentiate twice both sides of equation
total_energy_diff_eq = Eq(diff(total_energy_eq.lhs, time), diff(total_energy_eq.rhs, time))

## The second derivative of displacement is acceleration
spring_acceleration_derived_from_energy = solve(total_energy_diff_eq,
    Derivative(spring_displacement(time), (time, 2)),
    dict=True)[0][Derivative(spring_displacement(time), (time, 2))]
spring_acceleration_diff_eq = Eq(Derivative(spring_displacement(time), (time, 2)),
    spring_acceleration_derived_from_energy)

oscillator_eq = oscillator.definition.subs(oscillator.time, time)
oscillator_eq = oscillator_eq.subs(oscillator.displacement(time), spring_displacement(time))
angular_frequency_solved = simplify(
    solve([oscillator_eq, spring_acceleration_diff_eq],
    (oscillator.angular_frequency, spring_displacement(time)),
    dict=True)[0][oscillator.angular_frequency])

# 6. Derive period from frequency
period_law = period_definition.law.subs(period_definition.angular_frequency,
    angular_frequency_solved)
period_solved = solve(period_law, period_definition.period, dict=True)[0][period_definition.period]
## Square roots fail to compare with each other. Raise both parts to power of 2 before checking for equality.
assert expr_equals(period_solved**2, law.rhs**2)


@validate_input(stiffness_=stiffness, object_mass_=mass)
@validate_output(period)
def calculate_period(stiffness_: Quantity, object_mass_: Quantity) -> Quantity:
    solved = solve(law, period, dict=True)[0][period]
    result_expr = solved.subs({stiffness: stiffness_, mass: object_mass_})
    return Quantity(result_expr)
