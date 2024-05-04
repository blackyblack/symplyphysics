from sympy import Eq, sqrt, exp, I
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
    one_dimensional_wave_function_is_normalized as condition,
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
    lambda position_, time_: time_independent_wave_function(position_),
)

_time_dependent_condition = _condition.replace(
    condition.wave_function,
    lambda position_, time_: law.rhs.subs({
        position: position_,
        time: time_,
    })
)

# The latter condition actually reduces to the former one
assert expr_equals(_time_independent_condition.lhs, _time_dependent_condition.lhs)
assert expr_equals(_time_independent_condition.rhs, _time_dependent_condition.rhs)

# TODO: prove this law


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
