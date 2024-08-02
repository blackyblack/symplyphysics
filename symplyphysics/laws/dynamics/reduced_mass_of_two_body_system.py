"""
Reduced mass of a two-body system
=================================

Reduced mass is effective inertial mass in a system with two or more particles when they
are interacting with each other. This allowes the two-body problem to be solved as if it
were a one-body problem.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

reduced_mass = clone_symbol(symbols.basic.mass, "reduced_mass")
r"""
The reduced :attr:`~symplyphysics.symbols.basic.mass` of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

first_mass = clone_symbol(symbols.basic.mass, "first_mass")
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the first body.

Symbol:
    :code:`m1`

Latex:
    :math:`m_1`
"""

second_mass = clone_symbol(symbols.basic.mass, "second_mass")
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the second body.

Symbol:
    :code:`m2`

Latex:
    :math:`m_2`
"""

law = Eq(reduced_mass, 1 / (1 / first_mass + 1 / second_mass))
r"""
:code:`mu = 1 / (1 / m1 + 1 / m2)`

Latex:
    .. math::
        \mu = \frac{1}{\frac{1}{m_1} + \frac{1}{m_2}}
"""


@validate_input(
    first_mass_=first_mass,
    second_mass_=second_mass,
)
@validate_output(reduced_mass)
def calculate_reduced_mass(
    first_mass_: Quantity,
    second_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_mass: first_mass_,
        second_mass: second_mass_,
    })
    return Quantity(result)
