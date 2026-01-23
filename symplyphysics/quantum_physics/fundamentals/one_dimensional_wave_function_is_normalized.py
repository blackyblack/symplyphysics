"""
One-dimensional wave function is normalized
===========================================

For the wave function to be physically acceptable, it needs to be normalized, i.e. the integral of
the square of its absolute value must converge to one. The physical meaning of this is that the
particle, whose distribution in space is described by the wave function, must exists *somewhere*
in space.

**Links:**

#. `Physics LibreTexts, formula 3.6.3 <https://chem.libretexts.org/Courses/Pacific_Union_College/Quantum_Chemistry/03%3A_The_Schrodinger_Equation_and_a_Particle_in_a_Box/3.06%3A_Wavefunctions_Must_Be_Normalized>`__.
"""

from sympy import Eq, Integral, S
from symplyphysics import symbols, clone_as_function

position = symbols.position
"""
:symbols:`position` in the 1D space.
"""

time = symbols.time
"""
:symbols:`time`.
"""

wave_function = clone_as_function(symbols.wave_function, [position, time])
"""
:symbols:`wave_function` as a function of :attr:`~position` and :attr:`~time`.
"""

normalization_condition = Eq(
    Integral(abs(wave_function(position, time))**2, (position, S.NegativeInfinity, S.Infinity)),
    1,
)
"""
:laws:symbol::

:laws:latex::
"""
