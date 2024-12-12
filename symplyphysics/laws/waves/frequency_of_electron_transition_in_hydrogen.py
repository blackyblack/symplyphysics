"""
Frequency of electron transition in Hydrogen
============================================

The *Rydberg formula* relates the frequency of light emitted or absorbed during
electron transition to the levels between which the transition is occurring.

**Links:**

#. `Wikipedia, equivalent formula in terms of wavelength <https://phys.libretexts.org/Bookshelves/University_Physics/Book%3A_Introductory_Physics_-_Building_Models_to_Describe_Our_World_(Martin_Neary_Rinaldo_and_Woodman)/14%3A_Waves/14.07%3A_Standing_waves>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
    clone_as_symbol,
)

transition_frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of light emitted or absorbed during electron transition in the Hydrogen atom.
"""

lower_principal_quantum_number = clone_as_symbol(
    symbols.positive_number,
    display_symbol="n_1",
    display_latex="n_1",
)
"""
The principal quantum number of the lower transitional level. See :symbols:`positive_number`.
"""

higher_principal_quantum_number = clone_as_symbol(
    symbols.positive_number,
    display_symbol="n_2",
    display_latex="n_2",
)
"""
The principal quantum number of the higher transitional level.  See :symbols:`positive_number`.
"""

law = Eq(
    transition_frequency,
    quantities.rydberg_frequency * ((1 / lower_principal_quantum_number**2) -
    (1 / higher_principal_quantum_number**2)))
"""
.. only:: comment

    For now the auto-generation doesn't handle powers correctly here.

:code:`f = R_H * (1 / n_1^2 - 1 / n_2^2)`

:laws:latex::
"""


@validate_input(
    number_level_to_=lower_principal_quantum_number,
    number_level_from_=higher_principal_quantum_number,
)
@validate_output(transition_frequency)
def calculate_transition_frequency(number_level_to_: int, number_level_from_: int) -> Quantity:
    result_expr = solve(law, transition_frequency, dict=True)[0][transition_frequency]
    transition_frequency_applied = result_expr.subs({
        lower_principal_quantum_number: number_level_to_,
        higher_principal_quantum_number: number_level_from_
    })
    return Quantity(transition_frequency_applied)
