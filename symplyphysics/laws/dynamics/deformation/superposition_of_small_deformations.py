"""
Superposition of small deformations
===================================

If the deformations are small, the elastic properties of bodies stay the same under deformations.
If several forces act on a body, its total deformation can be computed from the deformations of
the individual forces, however these deformations must be small.

**Conditions:**

#. The deformations are small.

**Links:**

#. `Wikipedia, general information can be found here <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`__.
"""

from sympy import Eq
from symplyphysics import (
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

total_strain = symbols.strain
"""
Total :symbols:`strain` of the body.
"""

first_strain = clone_as_symbol(symbols.strain, display_symbol="e_1", display_latex="e_1")
"""
:symbols:`strain` caused by one force.
"""

second_strain = clone_as_symbol(symbols.strain, display_symbol="e_2", display_latex="e_2")
"""
:symbols:`strain` caused by another force.
"""

law = Eq(total_strain, first_strain + second_strain)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    first_strain_=first_strain,
    second_strain_=second_strain,
)
@validate_output(total_strain)
def calculate_total_strain(
    first_strain_: float,
    second_strain_: float,
) -> float:
    result = law.rhs.subs({
        first_strain: first_strain_,
        second_strain: second_strain_,
    })
    return convert_to_float(result)
