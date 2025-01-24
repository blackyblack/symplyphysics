"""
Corona discharge current from voltage
=====================================

The corona discharge is an independent discharge in a relatively dense gas. If an
electric field is applied to two electrodes between which there is a gas gap, then a
corona discharge occurs at a certain potential difference between the electrodes.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

current = symbols.current
"""
Corona discharge :symbols:`current`.
"""

experimental_coefficient = SymbolNew("A", units.charge / (units.area * units.voltage))
"""
Coefficient of the gas which is experimentally determined.
"""

mobility = symbols.mobility
"""
:symbols:`mobility` of charge carriers.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` applied.
"""

discharge_voltage = clone_as_symbol(symbols.voltage, subscript="0")
"""
:symbols:`voltage` at which the corona discharge occurs.
"""

law = Eq(
    current,
    experimental_coefficient * mobility * voltage *
    (voltage - discharge_voltage))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(gas_coefficient_=experimental_coefficient,
    mobility_of_charged_particles_=mobility,
    voltage_=voltage,
    corona_discharge_occurrence_voltage_=discharge_voltage)
@validate_output(current)
def calculate_current(gas_coefficient_: Quantity, mobility_of_charged_particles_: Quantity,
    voltage_: Quantity, corona_discharge_occurrence_voltage_: Quantity) -> Quantity:
    result_expr = solve(law, current, dict=True)[0][current]
    result_expr = result_expr.subs({
        experimental_coefficient: gas_coefficient_,
        mobility: mobility_of_charged_particles_,
        voltage: voltage_,
        discharge_voltage: corona_discharge_occurrence_voltage_,
    })
    return Quantity(result_expr)
