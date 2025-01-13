"""
Energy levels of harmonic oscillator
====================================

As opposed to the classical harmonic oscillator, the energy levels of a quantum harmonic oscillator
are quantized, meaning that its energy take a value out of a discrete range. These energy levels
are equidistant, i.e. the difference between successive energy levels is the same for all levels.

**Notation:**

#. :quantity_notation:`hbar`.

**Notes**

#. This means that the energy of a quantum oscillator cannot be zero and the lowest it can be
   is the zero-point energy :math:`E_0 = \\hbar \\omega / 2`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Hamiltonian_and_energy_eigenstates>`__.
"""

from sympy import Eq, Rational, solve, symbols as sym_symbols, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.quantum_mechanics.harmonic_oscillator import equation, wave_eigenfunctions

energy_level = clone_as_symbol(symbols.energy, subscript="n")
"""
Energy of the level corresponding to the :attr:`~mode_number`.
"""

mode_number = symbols.nonnegative_number
"""
Quantum number of oscillator, which is any non-negative integer (:math:`0, 1, 2, \\dots`).
See :symbols:`nonnegative_number`.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the oscillator.
"""

law = Eq(
    energy_level,
    (mode_number + Rational(1, 2)) * quantities.hbar * angular_frequency,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from Schrödinger equation and wave eigenfunction expressions

_position = equation.position
_mass = equation.particle_mass

_schrodinger_eqn = equation.law.subs({
    equation.angular_frequency: angular_frequency,
    equation.particle_energy: energy_level,
})

_eigenfunction_expr = wave_eigenfunctions.law.rhs.subs({
    wave_eigenfunctions.mode_number: mode_number,
    wave_eigenfunctions.oscillator_mass: _mass,
    wave_eigenfunctions.angular_frequency: angular_frequency,
    wave_eigenfunctions.position: _position,
})

_schrodinger_eqn = _schrodinger_eqn.replace(
    equation.wave_function,
    lambda position_: _eigenfunction_expr.subs(_position, position_),
).doit()

_reduced_position = sym_symbols("_reduced_position")

_reduced_position_eqn = Eq(
    _reduced_position,
    _position * sqrt(angular_frequency * _mass / quantities.hbar),
)

_energy_expr = solve(
    (_schrodinger_eqn, _reduced_position_eqn),
    (energy_level, _position),
    dict=True,
)[0][energy_level]

_hermite = wave_eigenfunctions.hermite

# See [this](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
_hermite_recurrence_relation = Eq(
    _hermite(mode_number + 1, _reduced_position),
    2 * _reduced_position * _hermite(mode_number, _reduced_position) -
    2 * mode_number * _hermite(mode_number - 1, _reduced_position),
)

_energy_expr = solve(
    (
    Eq(energy_level, _energy_expr),
    _hermite_recurrence_relation.subs(mode_number, mode_number - 1),
    ),
    (
    energy_level,
    _hermite(mode_number - 2, _reduced_position),
    ),
    dict=True,
)[0][energy_level].simplify()

assert expr_equals(_energy_expr, law.rhs)


@validate_input(
    mode_number_=mode_number,
    angular_frequency_=angular_frequency,
)
@validate_output(energy_level)
def calculate_energy_level(
    mode_number_: int,
    angular_frequency_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        mode_number: mode_number_,
        angular_frequency: angular_frequency_,
    })
    return Quantity(result)
