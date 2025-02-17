"""
Coefficient of volumetric ionization of neutral molecules by electrons
======================================================================

At a certain voltage, the gas discharge becomes independent. To find this voltage, it is
necessary to know the (volumetric) ionization coefficient. And it, in turn, depends on
the energy distribution of electrons and can be approximated by the expression below.

**Links:**

#. `StudFiles, formula 1.12 <https://studfile.net/preview/3079348/page:3/>`__.

..
    TODO: find English link
    TODO: move to `ionization` folder?
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

ionization_coefficient = symbols.ionization_coefficient
"""
:symbols:`ionization_coefficient` of the gas.
"""

first_constant = Symbol("A", 1 / units.length / units.pressure)
"""
The first gas constant used in this model.
"""

second_constant = Symbol("B", units.voltage / units.length / units.pressure)
"""
The second gas constant used in this model.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` of the gas.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

law = Eq(ionization_coefficient,
    first_constant * pressure * exp(-second_constant * pressure / electric_field_strength))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(first_constant_of_gas_=first_constant,
    second_constant_=second_constant,
    pressure_=pressure,
    electric_intensity_=electric_field_strength)
@validate_output(ionization_coefficient)
def calculate_ionization_coefficient(first_constant_of_gas_: Quantity, second_constant_: Quantity,
    pressure_: Quantity, electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, ionization_coefficient, dict=True)[0][ionization_coefficient]
    result_expr = result_expr.subs({
        first_constant: first_constant_of_gas_,
        second_constant: second_constant_,
        pressure: pressure_,
        electric_field_strength: electric_intensity_,
    })
    return Quantity(result_expr)
