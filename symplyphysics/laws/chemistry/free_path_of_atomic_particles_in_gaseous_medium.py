"""
Mean free path of particles in gaseous medium
=============================================

The atoms of the target material evaporate and move towards the substrate inside the
magnetron. At the same time, it collides with gas atoms. The free path length is the
distance that a traveling atom travels between two collisions.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Notes:**

#. Assuming the model of spherical gas molecules, :math:`\\sigma = pi * d^2`, where
   :math:`\\sigma` is the cross section and :math:`d` is the molecule diameter.

**Links:**

#. `Wikipedia, the fourth formula <https://en.wikipedia.org/wiki/Mean_free_path#Kinetic_theory_of_gases>`__.
#. `Chemistry LibreTexts, "27.6.4. Mean Free Path" <https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Physical_Chemistry_(LibreTexts)/27%3A_The_Kinetic_Theory_of_Gases/27.06%3A_Mean_Free_Path>`__.

..
    NOTE: remove the mention of a magnetron from the description?
    TODO: move to `magnetron` folder?
"""

from sympy import Eq, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

mean_free_path = symbols.mean_free_path
"""
:symbols:`mean_free_path` of particle.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

cross_section = symbols.cross_section
"""
:symbols:`cross_section` of the interaction between the particle and the gas.
"""

law = Eq(mean_free_path,
    quantities.boltzmann_constant * temperature / (sqrt(2) * pressure * cross_section))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(pressure_=pressure,
    temperature_=temperature,
    cross_sectional_area_of_interaction_=cross_section)
@validate_output(mean_free_path)
def calculate_free_path_length(pressure_: Quantity, temperature_: Quantity,
    cross_sectional_area_of_interaction_: Quantity) -> Quantity:
    result_expr = solve(law, mean_free_path, dict=True)[0][mean_free_path]
    result_expr = result_expr.subs({
        pressure: pressure_,
        temperature: temperature_,
        cross_section: cross_sectional_area_of_interaction_,
    })
    return Quantity(result_expr)
