r"""
Pressure of ideal gas from height and temperature
=================================================

The *barometric formula* determines the dependence of the pressure or density of a gas on the height in the gravity field.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.
#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is ideal.
#. The gas is in a uniform gravity field.
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
    clone_as_symbol,
)

final_pressure = symbols.pressure
"""
:symbols:`pressure` of the gas at final height.
"""

initial_pressure = clone_as_symbol(symbols.pressure, subscript="0")
"""
:symbols:`pressure` of the gas at initial height.
"""

molecular_mass = symbols.mass
"""
:symbols:`mass` of a single gas molecule.
"""

height_change = clone_as_symbol(
    symbols.height,
    display_symbol="Delta(h)",
    display_latex="\\Delta h",
)
"""
Change in :symbols:`height` between :attr:`~initial_pressure` and :attr:`~final_pressure`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

law = Eq(
    final_pressure,
    initial_pressure *
    exp(-quantities.acceleration_due_to_gravity * molecular_mass * height_change /
    (quantities.boltzmann_constant * temperature)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(initial_pressure_=initial_pressure,
    molecular_mass_=molecular_mass,
    height_change_=height_change,
    temperature_=temperature)
@validate_output(final_pressure)
def calculate_final_pressure(initial_pressure_: Quantity, molecular_mass_: Quantity,
    height_change_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, final_pressure, dict=True)[0][final_pressure]
    result_final_pressure = result_expr.subs({
        initial_pressure: initial_pressure_,
        molecular_mass: molecular_mass_,
        height_change: height_change_,
        temperature: temperature_
    })
    return Quantity(result_final_pressure)
