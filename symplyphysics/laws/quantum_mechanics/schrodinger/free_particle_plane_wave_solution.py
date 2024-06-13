from sympy import Eq, exp, I, S, dsolve
from sympy.physics.units import hbar
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    units,
    dimensionless,
    convert_to,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_via_momentum
from symplyphysics.laws.quantum_mechanics.schrodinger import (
    time_independent_equation_in_one_dimension as time_independent_eqn,
    time_dependent_via_time_independent_solution as time_dependent_law,
)

# Description
## In the absence of a potential field (U = 0) the wave function can be represented in the form of
## a [plane wave](https://en.wikipedia.org/wiki/Plane_wave).

# Law: Psi = exp((i / hbar) * (p * x - E * t))
## Psi - wave function
## x - position
## t - time
## hbar - reduced Planck constant
## p - momentum of quantum particle
## E - energy of quantum particle

# Notes
## - The energy and momentum of a quantum particle are related by the equation `E = p**2 / (2 * m)`
##   where `m` is the mass of the quantum particle.
## - The wave function must be normalized, i.e. the integral of the square of its absolute value
##   over the whole range of the spatial variable must equal 1, but the integral of the complex expontential
##   diverges if taken over the real line. This is not a problem, though, because the states described by
##   such a wave function would never be infinite (they would not be defined over the whole real line).
##   Moreover, the correct solution can be represented by a linear combination of planar waves, which can be
##   made convergent.

wave_function = Symbol("wave_function", dimensionless)
particle_momentum = Symbol("particle_momentum", units.momentum)
particle_energy = Symbol("particle_energy", units.energy)
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(wave_function, exp((I / hbar) * (particle_momentum * position - particle_energy * time)))

# Derive from Schrödinger equation

_time_independent_eqn = time_independent_eqn.law.subs(
    time_independent_eqn.position, position,
).subs({
    time_independent_eqn.potential_energy(position): 0,
    time_independent_eqn.particle_energy: kinetic_energy_via_momentum.law.rhs.subs({
        kinetic_energy_via_momentum.mass: time_independent_eqn.particle_mass,
        kinetic_energy_via_momentum.momentum: particle_momentum,
    }),
})

_time_independent_solution = dsolve(
    _time_independent_eqn,
    time_independent_eqn.wave_function(position),
).rhs.subs({
    "C1": 0,
    "C2": 1,
})

# The solution is a linear combination of exponents corresponding to positive and negative values of particle's momentum.
# If we can prove this law for any real value of momentum, it would be correct for any linear combination of solutions.
# The magnitude of the constants before the exponents does not affect the solution either. Therefore we can set the constant
# at the term with negative exponent to 0, and the constant at the term with positive exponent to 1.

_time_dependent_solution = time_dependent_law.law.rhs.replace(
    time_dependent_law.time_independent_wave_function,
    lambda _: _time_independent_solution,
).subs({
    time_dependent_law.particle_energy: particle_energy,
    time_dependent_law.time: time,
})

assert expr_equals(_time_dependent_solution, law.rhs)


@validate_input(
    particle_momentum_=particle_momentum,
    particle_energy_=particle_energy,
    position_=position,
    time_=time,
)
@validate_output(wave_function)
def calculate_wave_function_value(
    particle_momentum_: Quantity,
    particle_energy_: Quantity,
    position_: Quantity,
    time_: Quantity,
) -> complex:
    result = law.rhs.subs({
        particle_momentum: particle_momentum_,
        particle_energy: particle_energy_,
        position: position_,
        time: time_,
    })
    return complex(convert_to(result, S.One))
