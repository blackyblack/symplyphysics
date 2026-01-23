"""
Macroscopic cross section from mean free path
=============================================

Macroscopic cross section is the inverse of mean free path of neutrons. As an example,
the greater the absorption cross section of a material is, the shorter the mean free
path of neutrons is and the closer to the surface the neutron absorption will occur.

**Links:**

#. `NuclearPower, section "Mean Free Path" <https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-engineering-fundamentals/neutron-nuclear-reactions/macroscopic-cross-section/>`__.

..
    TODO: fix file name
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

mean_free_path = symbols.mean_free_path
"""
:symbols:`mean_free_path` of neutrons.
"""

macroscopic_cross_section = symbols.macroscopic_cross_section
"""
:symbols:`macroscopic_cross_section`.
"""

law = Eq(macroscopic_cross_section, 1 / mean_free_path)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mean_free_path_=mean_free_path)
@validate_output(macroscopic_cross_section)
def calculate_cross_section(mean_free_path_: Quantity) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_cross_section,
        dict=True)[0][macroscopic_cross_section]
    result_expr = result_cross_section_expr.subs(mean_free_path, mean_free_path_)
    return Quantity(result_expr)
