"""
Wave eigenfunctions of quantum harmonic oscillator
==================================================

The time-independent Schrödinger equation describing the wave function of the quantum oscillator
can be solved to get the corresponding wave eigenfunctions and energy eigenvalues of the Hamiltonian
operator. Each eigenfunction describes a stationary state of the quantum mechanical system with the
corresponding energy value (eigenvalue of the Hamiltonian). The combination of all eigenfunctions and
eigenvalues represent the energy states allowed.

**Notation:**

#. :quantity_notation:`hbar`.
#. :math:`H_n` (:code:`hermite`) is the :math:`n`-th physicists' Hermite polynomial.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Hamiltonian_and_energy_eigenstates>`__.
"""

from sympy import Eq, Rational, sqrt, exp, factorial, pi
from sympy.functions.special.polynomials import hermite
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.quantities import hbar

wave_function = symbols.wave_function
"""
:math:`n`-th wave eigenfunction (solution of the time-independent Schrödinger equation).
See :symbols:`wave_function`.
"""

mode_number = symbols.nonnegative_number
"""
Mode number. See :symbols:`nonnegative_number`.
"""

oscillator_mass = symbols.mass
"""
:symbols:`mass` of the oscillator.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the oscillator.
"""

position = symbols.position
"""
:symbols:`position` of the oscillator.
"""

law = Eq(wave_function,
    (1 / (sqrt(2**mode_number * factorial(mode_number)))) * ((oscillator_mass * angular_frequency) /
    (pi * hbar))**Rational(1, 4) * exp(-1 * (oscillator_mass * angular_frequency) /
    (2 * hbar) * position**2) * hermite(mode_number,
    sqrt(oscillator_mass * angular_frequency / hbar) * position))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    mode_number_=mode_number,
    oscillator_mass_=oscillator_mass,
    angular_frequency_=angular_frequency,
    position_=position,
)
@validate_output(wave_function)
def calculate_wave_function_value(
    mode_number_: int,
    oscillator_mass_: Quantity,
    angular_frequency_: Quantity,
    position_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        mode_number: mode_number_,
        oscillator_mass: oscillator_mass_,
        angular_frequency: angular_frequency_,
        position: position_,
    })
    return Quantity(result)
