"""
Interaction cross section in elastic interaction model
======================================================

The Sutherland formula adjusts the formula of the cross section of interacting gas molecules
taking into account their pairwise interactions. See :ref:`Effective cross section <effective_cross_section>`.

**Conditions:**

#. The sizes of the target and field molecules are equal.

#. Applicable to gases under moderate pressure.

**Links:**

#. `Wikipedia <https://ru.wikipedia.org/wiki/%D0%93%D0%B0%D0%B7%D0%BE%D0%BA%D0%B8%D0%BD%D0%B5%D1%82%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%BE%D0%B5_%D1%8D%D1%84%D1%84%D0%B5%D0%BA%D1%82%D0%B8%D0%B2%D0%BD%D0%BE%D0%B5_%D1%81%D0%B5%D1%87%D0%B5%D0%BD%D0%B8%D0%B5_%D0%BC%D0%BE%D0%BB%D0%B5%D0%BA%D1%83%D0%BB%D1%8B>`__.

..
    TODO: find English link
    NOTE: a more proper law would replace `particle_diameter` with `test_radius + target_radius`
"""

from sympy import Eq, solve, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols

cross_sectional_area_of_interaction = symbols.cross_section
"""
:symbols:`cross_section` of interaction of particles.
"""

particle_diameter = symbols.diameter
"""
:symbols:`diameter` of a gas particle.
"""

sutherland_constant = symbols.sutherland_constant
"""
:symbols:`sutherland_constant` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

law = Eq(cross_sectional_area_of_interaction,
    pi * particle_diameter**2 * (1 + sutherland_constant / temperature))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(diameter_of_atom_=particle_diameter,
    sutherland_constant_=sutherland_constant,
    temperature_=temperature)
@validate_output(cross_sectional_area_of_interaction)
def calculate_cross_sectional_area_of_interaction(diameter_of_atom_: Quantity,
    sutherland_constant_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, cross_sectional_area_of_interaction,
        dict=True)[0][cross_sectional_area_of_interaction]
    result_expr = result_expr.subs({
        particle_diameter: diameter_of_atom_,
        sutherland_constant: sutherland_constant_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
