from sympy import Eq, sqrt, Derivative
from sympy.physics.units import hbar
from symplyphysics import (
    units,
    Symbol,
    Function,
    symbols,
    clone_symbol,
    angle_type,
)

# Description
## The Schrödinger equation for the quantum simple harmonic oscillator governs the wave function
## of the quantum oscillator.

# Law: (-hbar**2 / (2 * m)) * d^2(psi(x))/dx^2 + (1 / 2) * m * w**2 * x**2 * psi(x) = E * psi(x)
## psi(x) - wave function of oscillating particle
## x - position
## hbar - reduced Planck constant
## m - particle mass
## w - angular frequency of oscillations
## E - energy eigenvalue of the equation

wave_function = Function("wave_function", 1 / sqrt(units.length))
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")
particle_energy = Symbol("particle_energy", units.energy)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
position = Symbol("position", units.length)

law = Eq(
    (-1 * hbar**2 / (2 * particle_mass)) * Derivative(wave_function(position), position, 2)
    + (particle_mass * angular_frequency**2 / 2) * position**2 * wave_function(position),
    particle_energy * wave_function(position),
)
