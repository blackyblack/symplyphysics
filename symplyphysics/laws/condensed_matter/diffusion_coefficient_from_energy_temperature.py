"""
Diffusion coefficient from energy and temperature
=================================================

The diffusion coefficient is a quantitative characteristic of the diffusion rate, equal
to the amount of matter passing per unit time through a section of a unit area as a
result of the thermal motion of molecules.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The material is a solid.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mass_diffusivity#Solids>`__.
"""

from sympy import (Eq, solve, exp)
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    quantities,
)

diffusion_coefficient = symbols.diffusion_coefficient
"""
:symbols:`diffusion_coefficient`.
"""

energy = clone_as_symbol(symbols.energy, subscript="\\text{A}")
"""
Activation :symbols:`energy` of the dopant.
"""

maximum_diffusion_coefficient = clone_as_symbol(symbols.diffusion_coefficient, subscript="0")
"""
Maximum :symbols:`diffusion_coefficient`, e.g. at infinite temperature.
"""

temperature = symbols.temperature
"""
:symbols:`temperature`.
"""

law = Eq(diffusion_coefficient,
    maximum_diffusion_coefficient * exp(-energy / (quantities.boltzmann_constant * temperature)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(energy_=energy, diffusion_constant_=maximum_diffusion_coefficient, temperature_=temperature)
@validate_output(diffusion_coefficient)
def calculate_diffusion_coefficient(energy_: Quantity, diffusion_constant_: Quantity,
    temperature_: Quantity) -> Quantity:
    result_expr = solve(law, diffusion_coefficient, dict=True)[0][diffusion_coefficient]
    result_expr = result_expr.subs({
        energy: energy_,
        maximum_diffusion_coefficient: diffusion_constant_,
        temperature: temperature_
    })
    return Quantity(result_expr)
