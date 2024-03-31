from sympy import (Derivative, Eq, Function as SymFunction, diff, solve, pi, sqrt, symbols as
    SymSymbols, simplify)
from symplyphysics import (Quantity, units, Symbol, print_expression, validate_input,
    validate_output, symbols)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import potential_energy_from_deformation as spring_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy
from symplyphysics.definitions import velocity_is_movement_derivative as velocity_def
from symplyphysics.laws.kinematic import period_from_angular_frequency as period_definition
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as oscillator

# Description
## Mass on spring is a system of object with mass m and spring with elasticity k. It starts oscillating after been pushed out of balance.

# Law: T = 2 * pi * sqrt(m/k)
# Where:
## T is period of oscillation,
## m is object mass,
## k is spring elastic coefficient.

# Conditions:
## - Spring doesn't accept or loose any energy. No any friction or inelastic deformations.
## - Object is small.
## - Another end of a spring is not moving in current coordinate system.
## - Spring is weightless.

spring_elasticity = Symbol("spring_elasticity", units.force / units.length)
oscillation_period = Symbol("oscillation_period", units.time)

law = Eq(oscillation_period, 2 * pi * sqrt(symbols.basic.mass / spring_elasticity))


def print_law() -> str:
    return print_expression(law)


# Derive this law from conservation of energy

## Spring displacement is distance between current and balanced positions.
spring_displacement = SymSymbols("spring_displacement", cls=SymFunction)
time = SymSymbols("time")

## Spring oscillation is cyclic transfer of energy from kinetic to potential. To set oscillation up we have to input some energy. Usually it is done by biasing the spring and letting it go.
## Biasing the spring is giving to it some amount of potential energy.

amount_of_potential_energy = spring_energy.law.subs({
    spring_energy.elastic_koefficient: spring_elasticity,
    spring_energy.deformation: spring_displacement(time)
}).rhs

## Kinetic energy of the pendulum is:
## object_mass * (linear_velocity)**2 / 2
velocity_def_eq = velocity_def.definition.subs(velocity_def.moving_time, time)
linear_velocity = velocity_def_eq.subs(velocity_def.movement(time), spring_displacement(time)).rhs
amount_of_kinetic_energy = kinetic_energy.law.subs({
    kinetic_energy.symbols.basic.mass: symbols.basic.mass,
    kinetic_energy.body_velocity: linear_velocity
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
oscillator_eq = oscillator_eq.subs(oscillator.displacement_function(time),
    spring_displacement(time))
angular_frequency_solved = simplify(
    solve([oscillator_eq, spring_acceleration_diff_eq],
    (oscillator.angular_frequency, spring_displacement(time)),
    dict=True)[0][oscillator.angular_frequency])

# 6. Derive period from frequency
period_law = period_definition.law.subs(period_definition.circular_frequency,
    angular_frequency_solved)
period_solved = solve(period_law, period_definition.period, dict=True)[0][period_definition.period]
## Square roots fail to compare with each other. Raise both parts to power of 2 before checking for equality.
assert expr_equals(period_solved**2, law.rhs**2)


@validate_input(spring_elasticity_=spring_elasticity, object_mass_=symbols.basic.mass)
@validate_output(oscillation_period)
def calculate_period(spring_elasticity_: Quantity, object_mass_: Quantity) -> Quantity:
    solved = solve(law, oscillation_period, dict=True)[0][oscillation_period]
    result_expr = solved.subs({
        spring_elasticity: spring_elasticity_,
        symbols.basic.mass: object_mass_
    })
    return Quantity(result_expr)
