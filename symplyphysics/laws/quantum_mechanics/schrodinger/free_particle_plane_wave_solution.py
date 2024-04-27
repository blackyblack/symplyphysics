from sympy import Eq, exp, I
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
## - This solution is divergent in terms of normalization over the whole real line.

wave_function = Symbol("wave_function", dimensionless)
particle_momentum = Symbol("particle_momentum", units.momentum)
particle_energy = Symbol("particle_energy", units.energy)
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(
    wave_function,
    exp((I / hbar) * (particle_momentum * position - particle_energy * time))
)

# TODO: derive from Schrödinger equation


@validate_input(
    particle_momentum_=particle_momentum,
    particle_energy_=particle_energy,
    position_=position,
    time_=time,
)
@validate_output(wave_function)
def calculate_wave_function(
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
    return complex(convert_to(result, 1))
