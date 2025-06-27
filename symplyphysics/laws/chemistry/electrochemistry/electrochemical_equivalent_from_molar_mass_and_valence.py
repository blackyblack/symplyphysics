"""
Electrochemical equivalent from molar mass and valence
======================================================

**Faraday's second law of electrolysis** states that if the same amount of electricity is passed
through different electrolytes, the masses of ions deposited at the electrodes are directly
proportional to their chemical equivalents.

**Notations:**

#. :quantity_notation:`faraday_constant`.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Faraday%27s_laws_of_electrolysis#Derivation>`__.

#. `BYJY'S <https://byjus.com/chemistry/laws-of-electrolysis/>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

electrochemical_equivalent = symbols.electrochemical_equivalent
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

law = Eq(electrochemical_equivalent, molar_mass / (quantities.faraday_constant * valence))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(molar_mass_=molar_mass, valence_=valence)
@validate_output(electrochemical_equivalent)
def calculate_equivalent(molar_mass_: Quantity, valence_: int) -> Quantity:
    if not isinstance(valence_, int):
        raise ValueError("valence_ must be an integer.")
    if valence_ <= 0:
        raise ValueError("valence_ must be greater than 0.")

    result_expr = solve(law, electrochemical_equivalent, dict=True)[0][electrochemical_equivalent]
    result_expr = result_expr.subs({
        molar_mass: molar_mass_,
        valence: valence_,
    })
    return Quantity(result_expr)
