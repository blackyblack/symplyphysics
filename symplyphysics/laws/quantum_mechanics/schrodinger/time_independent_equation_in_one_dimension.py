from sympy import Eq, sqrt, Derivative
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
## function of a quantum-mechanical system.

# Law: (-hbar**2 / (2 * m)) * d^2(psi(x))/dx^2 + U(x) * psi(x) = E * psi(x)
## psi(x) - wave function of particle
## x - position
## hbar - reduced Planck constant
## m - particle mass
## U(x) - potential energy
## E - energy of particle, eigenvalue of the equation

# Conditions
## - This law works in the case of 1 spatial dimension. To use it for the 3-dimensional space
##   replace the second derivative with the Laplace operator.

wave_function = Function("wave_function", 1 / sqrt(units.length))
position = Symbol("position", units.length)
potential_energy = Function("potential_energy", units.energy)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")
particle_energy = Symbol("particle_energy", units.energy)

law = Eq(
    (-1 * hbar**2 / (2 * particle_mass)) * Derivative(wave_function(position), position, 2)
    + potential_energy(position) * wave_function(position),
    particle_energy * wave_function(position),
)

# The solutions of the Schrödinger equation vary drastically with the kind of potential energy
# function used and also depend on the boundary conditions.
