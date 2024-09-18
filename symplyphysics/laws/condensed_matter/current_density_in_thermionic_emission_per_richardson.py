r"""
Current density in thermionic emission per Richardson
=====================================================

Thermionic emission is the liberation of electrons from an electrode by virtue of its temperature.
This occurs because the thermal energy given to the charge carrier overcomes the work function of
the material. This formula for emission current density was proposed by Owen Williams Richardson.

**Notation:**

#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.
"""

from sympy import (Eq, solve, exp)
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

current_density = Symbol("current_density", units.current / units.area)
"""
Emission current density.

Symbol:
    :code:`j`
"""

work_function = Symbol("work_function", units.energy)
"""
Work function of the material.

Symbol: 
    :code:`W`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the metal.
"""

richardson_constant = Quantity(120.17 * (units.ampere / units.kelvin**2 / units.centimeter**2), display_symbol="A")
"""
Constant of proportionality proposed by Richardson.
"""

law = Eq(
    current_density,
    richardson_constant * temperature**2 * exp(-1 * work_function / (units.boltzmann_constant * temperature)))
r"""
:code:`j = A * T^2 * exp(-1 * W / (k_B * T))`

Latex:
    .. math::
        j = A T^2 \exp \left( - \frac{W}{k_\text{B} T} \right)
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
