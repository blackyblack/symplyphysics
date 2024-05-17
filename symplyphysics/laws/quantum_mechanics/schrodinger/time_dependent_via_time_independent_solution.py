from sympy import Eq, sqrt, exp, I, solve
from sympy.physics.units import hbar
from symplyphysics import (
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    units,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.conditions.quantum_mechanics import (
    one_dimensional_wave_function_is_normalized as condition,)
from symplyphysics.laws.quantum_mechanics.schrodinger import (
    time_independent_equation_in_one_dimension as stationary_eqn,
    time_dependent_equation_in_one_dimension as general_eqn,
)

# Description
## When the potential energy is independent of time, the solution of the time-dependent Schrödinger
## equation can be constructed using the solution of the time-independent equation.

# Law: Psi(x, t) = psi(x) * exp(-1 * (i / hbar) * E * t)
## Psi(x, t) - solution of time-dependent Schrödinger equation
## psi(x) - solution of time-independent Schrödinger equation
## x - position, spatial variable
## t - time
## i - imaginary unit
## hbar - reduced Planck constant
## E - particle energy

# Conditions
## - The potential energy is a function independent of time.
## - This law applies for 1-dimensional systems.

time_dependent_wave_function = Function("time_dependent_wave_function", 1 / sqrt(units.length))
time_independent_wave_function = Function("time_independent_wave_function", 1 / sqrt(units.length))
particle_energy = Symbol("particle_energy", units.energy, positive=True)
position = Symbol("position", units.length, real=True)
time = Symbol("time", units.time, real=True)

law = Eq(
    time_dependent_wave_function(position, time),
    time_independent_wave_function(position) * exp(-1 * (I / hbar) * particle_energy * time),
)

# Show that the normalization condition still holds

_condition = condition.normalization_condition.subs({
    condition.position: position,
    condition.time: time,
})

_time_independent_condition = _condition.replace(
    condition.wave_function,
    lambda position_, _time: time_independent_wave_function(position_),
)

_time_dependent_condition = _condition.replace(
    condition.wave_function, lambda position_, time_: law.rhs.subs({
    position: position_,
    time: time_,
    }))

# The latter condition actually reduces to the former one
assert expr_equals(_time_independent_condition.lhs, _time_dependent_condition.lhs)
assert expr_equals(_time_independent_condition.rhs, _time_dependent_condition.rhs)

# Prove this law

_stationary_eqn = stationary_eqn.law.subs({
    stationary_eqn.position: position,
    stationary_eqn.particle_energy: particle_energy,
}).replace(stationary_eqn.wave_function, time_independent_wave_function)

# Substitute the time-dependent wave function defined in this law into the general Schrödinger equation

_general_eqn = general_eqn.law.subs({
    general_eqn.position: position,
    general_eqn.time: time,
    general_eqn.particle_mass: stationary_eqn.particle_mass,
}).replace(
    general_eqn.wave_function,
    lambda position_, time_: law.rhs.subs(position, position_),
).replace(
    general_eqn.potential_energy,
    lambda position_, time_: stationary_eqn.potential_energy(position_),
).doit()

# Solve for the stationary wave function in both equations and compare them

_wave_function_of_stationary_eqn = solve(
    _stationary_eqn,
    time_independent_wave_function(position),
)[0]

_wave_function_of_general_eqn = solve(_general_eqn, time_independent_wave_function(position))[0]

# The two expressions happen to be equal, which means that the wave function proposed in this
# law is indeed the solution to the general Schödinger equation.

assert expr_equals(_wave_function_of_stationary_eqn, _wave_function_of_general_eqn)


@validate_input(
    time_independent_wave_function_value_=time_independent_wave_function,
    particle_energy_=particle_energy,
    time_=time,
)
@validate_output(time_dependent_wave_function)
def calculate_time_dependent_wave_function_value(
    time_independent_wave_function_value_: Quantity,
    particle_energy_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        time_independent_wave_function(position): time_independent_wave_function_value_,
        particle_energy: particle_energy_,
        time: time_,
    })
    return Quantity(result)
