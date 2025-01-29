"""
Mass of film deposited during electrolysis
==========================================

Electrolysis is a physico-chemical process consisting in the release of components of
dissolved substances or other substances on the electrodes, which are the result of
secondary reactions on the electrodes, which occurs when an electric current passes
through a solution or melt of an electrolyte.

**Notation:**

#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`avogadro_constant`.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, SymbolNew, dimensionless
from symplyphysics.quantities import elementary_charge, avogadro_constant

film_mass = symbols.mass
"""
:symbols:`mass` of the film deposited.
"""

current = symbols.current
"""
:symbols:`current` in the electrolyte.
"""

molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass` of deposited metal.
"""

current_output = SymbolNew("B", dimensionless)
"""
Current output is the fraction of electric current spent on the passage of the target electrochemical reaction.
"""

valence = symbols.valence
"""
:symbols:`valence` of the deposited metal.
"""

time = symbols.time
"""
:symbols:`time` of electrolysis.
"""

law = Eq(
    film_mass,
    current * molar_mass * current_output * time /
    (valence * elementary_charge * avogadro_constant),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(current_=current,
    molar_mass_=molar_mass,
    current_output_=current_output,
    valence_=valence,
    time_=time)
@validate_output(film_mass)
def calculate_mass_of_film(current_: Quantity, molar_mass_: Quantity, current_output_: float,
    valence_: int, time_: Quantity) -> Quantity:
    result_expr = solve(law, film_mass, dict=True)[0][film_mass]
    result_expr = result_expr.subs({
        current: current_,
        molar_mass: molar_mass_,
        current_output: current_output_,
        valence: valence_,
        time: time_,
    })
    return Quantity(result_expr)
