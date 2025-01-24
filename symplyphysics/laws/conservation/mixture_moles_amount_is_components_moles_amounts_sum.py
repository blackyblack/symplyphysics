"""
Amount of mixture is sum of amounts of components
=================================================

The total amount of the mixture is equal to the sum of the amounts of each of the components.

**Conditions:**

#. Mixture in a closed impenetrable volume, that is, molecules/atoms cannot leave it and
   they are always inside;
#. Mass is not transformed to energy, for example due to annihilation.

**Links:**

#. `Engineering LibreTexts, Composition on a Molar Basis <https://eng.libretexts.org/Bookshelves/Introductory_Engineering/Basic_Engineering_Science_-_A_Systems_Accounting_and_Modeling_Approach_(Richards)/03%3A_Conservation_of_Mass/3.04%3A_Mixture_Composition>`__.

..
    TODO fix file name
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    SumIndexed,
    global_index,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

amount_of_mixture = symbols.amount_of_substance
"""
:symbols:`amount_of_substance` in the mixture.
"""

amount_of_component = clone_as_indexed(symbols.amount_of_substance)
"""
:symbols:`amount_of_substance` in the :math:`i`-th component.
"""

law = Eq(amount_of_mixture, SumIndexed(amount_of_component[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(amount_of_components_=amount_of_component)
@validate_output(amount_of_mixture)
def calculate_moles_count_of_mixture(amount_of_components_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(amount_of_components_)))
    moles_counts_of_mixture_law = law.subs(global_index, local_index)
    moles_counts_of_mixture_law = moles_counts_of_mixture_law.doit()
    solved = solve(moles_counts_of_mixture_law, amount_of_mixture,
        dict=True)[0][amount_of_mixture]
    for i, v in enumerate(amount_of_components_):
        solved = solved.subs(amount_of_component[i + 1], v)
    return Quantity(solved)
