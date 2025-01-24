"""
Electrochemical equivalent from molar mass and valence
======================================================

Faraday's second law of electrolysis. The equivalent mass of a substance in general in
chemistry is its molar mass divided by an integer depending on the chemical reaction in
which the substance participates; in this case, the equivalent is the molar mass of the
substance formed during ion discharge divided by the sum of the ion charges (measured
in elementary units), resulting in a molecule or atom of the substance.

**Links:**

#. `Wikipedia, derivable <https://en.wikipedia.org/wiki/Electrochemical_equivalent>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

equivalent = symbols.electrochemical_equivalent
"""
:symbols:`electrochemical_equivalent`.
"""

molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass`.
"""

valence = symbols.valence
"""
:symbols:`valence`.
"""

law = Eq(equivalent, molar_mass / (quantities.faraday_constant * valence))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(molar_mass_=molar_mass, valence_=valence)
@validate_output(equivalent)
def calculate_equivalent(molar_mass_: Quantity, valence_: int) -> Quantity:
    if not isinstance(valence_, int):
        raise ValueError("valence_ must be an integer.")
    if valence_ <= 0:
        raise ValueError("valence_ must be greater than 0.")

    result_expr = solve(law, equivalent, dict=True)[0][equivalent]
    result_expr = result_expr.subs({
        molar_mass: molar_mass_,
        valence: valence_,
    })
    return Quantity(result_expr)
