"""
Frequency of electron transition in Hydrogen
============================================

The *Rydberg formula* relates the frequency of light emitted or absorbed during
electron transition to the levels between which the transition is occurring.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, dimensionless, validate_input, validate_output)

transition_frequency = Symbol("transition_frequency", units.frequency)
r"""
Frequency of light emitted or absorbed during electron transition in the Hydrogen atom.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

lower_principal_quantum_number = Symbol("lower_principal_quantum_number", dimensionless)
r"""
The principal quantum number of the lower transitional level.

Symbol:
    :code:`n_1`

Latex:
    :math:`n_1`
"""

higher_principal_quantum_number = Symbol("higher_principal_quantum_number", dimensionless)
"""
The principal quantum number of the higher transitional level.

Symbol:
    :code:`n_2`

Latex:
    :math:`n_2`
"""

rydberg_constant = Quantity(3.2898419602500e15 * units.hertz,
    display_symbol="R_H",
    display_latex="R_\\text{H}")
"""
Rydberg frequency constant for Hydrogen.
"""

law = Eq(
    transition_frequency,
    rydberg_constant * ((1 / lower_principal_quantum_number**2) -
    (1 / higher_principal_quantum_number**2)))
r"""
:code:`nu = R_H * (1 / n_1^2 - 1 / n_2^2)`

Latex:
    .. math::
        \nu = R_\text{H} \left( \frac{1}{n_1^2} - \frac{1}{n_2^2} \right)
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
