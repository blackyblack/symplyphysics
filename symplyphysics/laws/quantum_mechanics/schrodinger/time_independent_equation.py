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
## psi(x) - wavefunction of particle
## x - position
## hbar - reduced Planck constant
## m - particle mass
## U(x) - potential energy
## E - energy of particle, eigenvalue of the equation

wavefunction = Function("wavefunction", 1 / sqrt(units.length), real=1)
position = Symbol("position", units.length, real=1)
potential_energy = Function("potential_energy", units.energy, real=1)
particle_mass = clone_symbol(symbols.basic.mass, "particle^mass", positive=1)
particle_energy = Symbol("particle_energy", units.energy, positive=1)

law = Eq(
    (-1 * hbar**2 / (2 * particle_mass)) * Derivative(wavefunction(position), position, 2)
    + potential_energy(position) * wavefunction(position),
    particle_energy * wavefunction(position),
)

# The solutions of the Schrödinger equation vary drastically with the kind of potential energy
# function is used and depend on the boundary conditions.
