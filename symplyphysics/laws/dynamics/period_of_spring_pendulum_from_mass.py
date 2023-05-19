from sympy import (Derivative, Eq, Function as SymFunction, cos, diff, sin, solve, pi, sqrt, symbols)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import potential_energy_from_deformation as spring_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy
from symplyphysics.laws.dynamics import spring_reaction_from_deformation as hook_law

# Description
## Spring pendulum is a system of object with mass m and spring with elastice k. It starts oscillating after been pushed out of balance.

# Law: T = 2 * pi * sqrt(m/k)
# Where:
## T is period of oscillation,
## m is pendulum mass,
## k is pendulum elastic coefficient.

# Conditions:
## - Pendulum doesn't accept or loose any energy. No any friction or inelastic deformations.
## - Object is small.
## - Another end of a spring is not moving in current coordinate system.
## - Spring is weightless.

pendulum_mass = Symbol("pendulum_mass", units.mass)
pendulum_elastice = Symbol("pendulum_elastice", units.force / units.length)
oscillation_period = Symbol("oscillation_period", units.time)

law = Eq(oscillation_period, 2 * pi * sqrt(pendulum_mass / pendulum_elastice))

# Derive this law from conservation of energy

## Pendulum displacement is distance between current and balanced positions.
pendulum_displacement = symbols("pendulum_displacement", cls=SymFunction)
time = symbols("time")

## Pendulum oscillation is cyclic transfer of energy from kinetic to potential. To set oscillation up we have to input some energy. Usually it is done by biasing the pendulum and letting it go.
## Biasing the pendulum is giving to it some amount of potential energy.

amount_of_potential_energy = spring_energy.law.subs({    
    spring_energy.elastic_koefficient: pendulum_elastice,
    spring_energy.deformation: pendulum_displacement
}).rhs

## Kinetic energy of the pendulum is:
## pendulum_mass * (linear_velocity)**2 / 2
linear_velocity = symbols("linear_velocity", cls=SymFunction)
amount_of_kinetic_energy = kinetic_energy.law.subs({
    kinetic_energy.body_mass: pendulum_mass,
    kinetic_energy.body_velocity: linear_velocity
}).rhs

## Total energy is constant and any of it's derivatives is 0.
total_energy = symbols("total_energy", constant=True)
total_energy_eq = Eq(total_energy, amount_of_kinetic_energy + amount_of_potential_energy)

## Differentiate twice both sides of equation
total_energy_diff_eq = Eq(0, diff(total_energy_eq.rhs, (time, 2)))

## The second derivative of displacement is acceleration
pendulum_acceleration_derived_from_energy = solve(total_energy_diff_eq, Derivative(pendulum_displacement(time), (time, 2)), dict=True)[0][Derivative(pendulum_displacement(time), (time, 2))]

## On the other hand the pendulum acceleration is caused by spring reaction force and Newton's the second law.
pendulum_acceleration_derived_from_force = hook_law.law.rhs.subs({hook_law.law.elastic_coefficient: pendulum_elastice, hook_law.law.deformation: pendulum_displacement}) / pendulum_mass

main_equation = Eq(pendulum_acceleration_derived_from_energy, pendulum_acceleration_derived_from_force)
displacement_function = solve(main_equation, pendulum_displacement(time), dict = True)[0][pendulum_displacement(time)]



'''
# For small angles, sin(pendulum_angle) can be reduced to pendulum_angle
small_angle_harmonic_oscillation = total_energy_diff_solved.subs(sin(pendulum_angle(time)), pendulum_angle(time))

# Will result in harmonic oscillator equation:
## Derivative(pendulum_angle(time), (time, 2)) = -free_fall_acceleration / pendulum_length * pendulum_angle(time)
harmonic_oscillation_eq = Eq(Derivative(pendulum_angle(time), (time, 2)), small_angle_harmonic_oscillation)

## dsolve() gives us solution in exponential form - we are looking for the solution as trigonometric function

maximum_angle = symbols("maximum_angle")
initial_phase = symbols("initial_phase")
oscillation_angular_frequency = sqrt(free_fall_acceleration / pendulum_length)
angle_function_eq = Eq(pendulum_angle(time), maximum_angle * cos(oscillation_angular_frequency * time + initial_phase))
dsolved = harmonic_oscillation_eq.subs(pendulum_angle(time), angle_function_eq.rhs)
assert expr_equals(dsolved.lhs, dsolved.rhs)

## There are many solutions for harmonic_oscillation_eq. Add condition, that at initial point of time (time = 0)
## pendulum is at its highest point (pendulum_angle(time) = maximum_angle).
## Let's prove that initial phase of cosine function (angle_function_eq) should be zero.

initial_pendulum_condition = Eq(pendulum_angle(0), maximum_angle)
angle_function_at_zero_time_eq = angle_function_eq.subs(time, 0)
## Initial phase solutions have period of 2*pi. Take first solution.
initial_phase_solved = solve([angle_function_at_zero_time_eq, initial_pendulum_condition], (maximum_angle, initial_phase), dict=True)[0][initial_phase]
assert expr_equals(initial_phase_solved, 0)

## Check that expected period matches our law.
## Square roots fail to compare with each other. Raise both parts to power of 2 before checking for equality.
oscillation_period_derived = angular_frequency.law.subs(angular_frequency.circular_frequency, oscillation_angular_frequency).rhs
assert expr_equals(oscillation_period_derived**2, law.rhs**2)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(pendulum_length_=pendulum_length)
@validate_output_symbol(oscillation_period)
def calculate_period(pendulum_length_: Quantity) -> Quantity:
    solved = solve(law, oscillation_period, dict=True)[0][oscillation_period]
    result_expr = solved.subs(pendulum_length, pendulum_length_)
    return expr_to_quantity(result_expr)
'''