r"""
Current density in thermionic emission per Richardson
=====================================================

Thermionic emission is the liberation of electrons from an electrode by virtue of its temperature.
This occurs because the thermal energy given to the charge carrier overcomes the work function of
the material. This formula for emission current density was proposed by Owen Williams Richardson.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.
"""

from sympy import (Eq, solve, exp)
from symplyphysics import (
    symbols,
    units,
    quantities,
    Quantity,
    validate_input,
    validate_output,
)

current_density = symbols.current_density
"""
Emission :symbols:`current_density`.
"""

work_function = symbols.work_function
"""
:symbols:`work_function` of the material.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the metal.
"""

richardson_constant = Quantity(
    120.17 * (units.ampere / units.kelvin**2 / units.centimeter**2),
    display_symbol="a")
"""
Constant of proportionality proposed by Richardson.
"""

law = Eq(
    current_density,
    richardson_constant * temperature**2 * exp(-1 * work_function / (quantities.boltzmann_constant * temperature)))
r"""
:laws:symbol::

:laws:latex::
"""


@validate_input(thermodynamic_work_=work_function, temperature_=temperature)
@validate_output(current_density)
def calculate_current(thermodynamic_work_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, current_density, dict=True)[0][current_density]
    result_expr = result_expr.subs({
        work_function: thermodynamic_work_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
