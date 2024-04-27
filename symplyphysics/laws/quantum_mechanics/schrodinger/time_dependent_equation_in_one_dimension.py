from sympy import Eq, sqrt, Derivative, I
from sympy.physics.units import hbar
from symplyphysics import (
    units,
    Symbol,
    Function,
    symbols,
    clone_symbol,
)

# Description
## The Schrödinger equation is a linear partial differential equation that governs the wave
## function of a quantum-mechanical system. This law describes the general case of a time-dependent
## potential and a time-dependent wave function.

# Notes
## - This law works in the case of 1 spatial dimension. To use it for the 3-dimensional space
##   replace the spatial second derivative with the Laplace operator.

wave_function = Function("wave_function", 1 / sqrt(units.length))
potential_energy = Function("potential_energy", units.energy)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")
position = Symbol("position", units.length)
time = Symbol("time", units.time)

law = Eq(
    -1 * hbar**2 / (2 * particle_mass) * Derivative(wave_function(position, time), position, 2)
    + potential_energy(position, time) * wave_function(position, time),
    I * hbar * Derivative(wave_function(position, time), time)
)
