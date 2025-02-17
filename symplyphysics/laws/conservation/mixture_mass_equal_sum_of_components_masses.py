"""
Mixture mass equals sum of component masses
===========================================

The mass of a mixture of fluids is equal to the sum of the masses of the components of
the mixture.

**Conditions:**

#. Mixture in a closed impenetrable volume, that is, its molecules cannot leave it and
   they are always inside;
#. Mass is not transformed to energy, for example due to annihilation.

**Links:**

#. `Engineering LibreTexts, Composition on a Mass Basis <https://eng.libretexts.org/Bookshelves/Introductory_Engineering/Basic_Engineering_Science_-_A_Systems_Accounting_and_Modeling_Approach_(Richards)/03%3A_Conservation_of_Mass/3.04%3A_Mixture_Composition>`__.

..
    TODO: fix file name
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    IndexedSum,
    global_index,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

mixture_mass = symbols.mass
"""
:symbols:`mass` of the mixture.
"""

component_mass = clone_as_indexed(symbols.mass)
"""
:symbols:`mass` of the :math:`i`-th component.
"""

law = Eq(mixture_mass, IndexedSum(component_mass[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(masses_of_components_=component_mass)
@validate_output(mixture_mass)
def calculate_mass_of_mixture(masses_of_components_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(masses_of_components_)))
    masses_of_components_law = law.subs(global_index, local_index)
    masses_of_components_law = masses_of_components_law.doit()
    solved = solve(masses_of_components_law, mixture_mass, dict=True)[0][mixture_mass]
    for i, v in enumerate(masses_of_components_):
        solved = solved.subs(component_mass[i + 1], v)
    return Quantity(solved)
