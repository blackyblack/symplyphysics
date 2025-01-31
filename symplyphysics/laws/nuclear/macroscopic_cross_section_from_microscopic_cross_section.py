"""
Macroscopic cross section from microscopic cross section and number density
===========================================================================

Macroscopic cross section is equal to the product of the microscopic cross section, i.e.
the effective target area where the nucleus interacts with an incident neutron, and the
atomic number density (concentration).

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-engineering-fundamentals/neutron-nuclear-reactions/macroscopic-cross-section/>`__.
#. `ScienceDirect <https://www.sciencedirect.com/topics/engineering/macroscopic-cross-section>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

microscopic_cross_section = symbols.cross_section
"""
Microscopic :symbols:`cross_section`.
"""

number_density = symbols.number_density
"""
:symbols:`number_density` of particles.
"""

macroscopic_cross_section = symbols.macroscopic_cross_section
"""
:symbols:`macroscopic_cross_section`.
"""

law = Eq(macroscopic_cross_section, microscopic_cross_section * number_density)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(microscopic_cross_section_=microscopic_cross_section,
    atomic_number_density_=number_density)
@validate_output(macroscopic_cross_section)
def calculate_cross_section(microscopic_cross_section_: Quantity,
    atomic_number_density_: Quantity) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_cross_section,
        dict=True)[0][macroscopic_cross_section]
    result_expr = result_cross_section_expr.subs({
        microscopic_cross_section: microscopic_cross_section_,
        number_density: atomic_number_density_
    })
    return Quantity(result_expr)
