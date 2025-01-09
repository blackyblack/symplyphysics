"""
Expectation value of observable
===============================

Roughly speaking, an observable is a measurable property of a physical system, such as spin, position,
energy, momentum, etc. In terms of Quantum Mechanics, every physical observable corresponds to an operator
which acts on the wave function.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Operator_(physics)#Expectation_values_of_operators_on_%CE%A8>`__.
"""

from sympy import Eq, Integral, conjugate, S, Expr
from symplyphysics import clone_as_symbol, clone_as_function, symbols
from symplyphysics.core.type_aliases import WaveFunction, Observable

mean_observable_value = clone_as_symbol(
    symbols.any_quantity,
    display_symbol="avg(O)",
    display_latex="\\langle O \\rangle",
)
"""
The mean value of the :attr:`~observable` operator. See :symbols:`any_quantity`.
"""

position = symbols.position
"""
:symbols:`position`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

wave_function = clone_as_function(symbols.wave_function, [position, time])
"""
:symbols:`wave_function` as a function of :attr:`~position` and :attr:`~time`.
"""

observable = clone_as_function(
    symbols.any_quantity,
    [wave_function(position, time)], # should technically be a function of unapplied wave function
    display_symbol="O",
    display_latex="O",
)
"""
Observable operator as a function of :attr:`~wave_function`. See :symbols:`any_quantity`.
"""

law = Eq(
    mean_observable_value,
    Integral(
    conjugate(wave_function(position, time)) * observable(wave_function(position, time)),
    (position, S.NegativeInfinity, S.Infinity)),
)
r"""
:laws:symbol::

..
    The conjugation operator isn't recognized yet by the LaTeX printer.

.. math::
    
    \langle O \rangle = \int \limits_{-\infty}^{\infty} \psi^* (x, t) O[\psi](x, t) dx
"""


def calculate_mean_observable_value(
    observable_: Observable,
    wave_function_: WaveFunction,
) -> Expr:
    result = law.rhs.replace(
        wave_function,
        wave_function_,
    ).replace(
        observable,
        lambda _wave_function: observable_(wave_function_)(position, time),
    ).doit()
    return result
