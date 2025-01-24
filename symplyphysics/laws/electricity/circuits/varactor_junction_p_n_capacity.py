"""
Capacitance of p-n varactor junction
====================================

The principle of operation of the variator is based on the dependence of the barrier
capacitance of the p-n junction on the voltage value. Knowing the doping coefficient,
the material parameter and the capacity without bias voltage, it is possible to
calculate the capacity of the varactor at a given voltage.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol,
)

capacitance = symbols.capacitance
"""
Barrier :symbols:`capacitance` of the p-n junction.
"""

capacitance_without_bias_voltage = clone_as_symbol(symbols.capacitance, subscript="0")
"""
Barrier :symbols:`capacitance` of the p-n junction without the bias voltage.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` applied.
"""

material_parameter = clone_as_symbol(symbols.voltage, subscript="0")
"""
A certain constant having the dimension of :symbols:`voltage` which depends on the exact
material.
"""

doping_coefficient = SymbolNew("y", dimensionless)
"""
Doping coefficient.
"""

law = Eq(
    capacitance, capacitance_without_bias_voltage /
    (1 - voltage / material_parameter)**doping_coefficient)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    junction_capacitance_without_bias_voltage_=capacitance_without_bias_voltage,
    voltage_=voltage,
    material_parameter_=material_parameter,
    doping_coefficient_=doping_coefficient)
@validate_output(capacitance)
def calculate_junction_capacitance(junction_capacitance_without_bias_voltage_: Quantity,
    voltage_: Quantity, material_parameter_: Quantity, doping_coefficient_: float) -> Quantity:
    result_expr = solve(law, capacitance, dict=True)[0][capacitance]
    result_expr = result_expr.subs({
        capacitance_without_bias_voltage: junction_capacitance_without_bias_voltage_,
        voltage: voltage_,
        material_parameter: material_parameter_,
        doping_coefficient: doping_coefficient_,
    })
    return Quantity(result_expr)
