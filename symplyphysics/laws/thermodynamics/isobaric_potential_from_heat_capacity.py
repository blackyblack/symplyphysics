"""
Isobaric potential from heat capacity
=====================================

The isobaric potential of a reaction is a value whose change during a chemical reaction
is equal to the change in the internal energy of the system. The isobaric potential
shows how much of the total internal energy of the system can be used for chemical
transformations. Thermal effect of reaction is change of enthalpy of the system. The
standard state is the state at a temperature of :math:`298 \\, \\text{K}` and a total
pressure of :math:`1 \\, \\text{atm}`, as well as at a fixed composition of the system.

**Notation:**

#. :quantity_notation:`standard_laboratory_temperature`.

**Conditions:**

#. We neglect the temperature dependence of the heat capacities.
#. Pressure is held constant during the process.
#. Temperature is held constant during the process.

..
    TODO: find link
"""

from sympy import (Eq, solve, log)
from symplyphysics import (
    quantities,
    symbols,
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    clone_as_symbol,
)

standard_molar_gibbs_energy_change = SymbolNew("Delta(G_m)", units.energy / units.amount_of_substance, display_latex="\\Delta G_\\text{m}")
"""
Standard change of :symbols:`gibbs_energy` per unit :symbols:`amount_of_substance`.
"""

standard_molar_enthalpy_change = SymbolNew("Delta(H_m)", units.energy / units.amount_of_substance, display_latex="\\Delta H_\\text{m}")
"""
Standard change of :symbols:`enthalpy` per unit :symbols:`amount_of_substance`.
"""

standard_molar_entropy_change = SymbolNew("Delta(S_m)", units.energy / units.amount_of_substance / units.temperature, display_latex="\\Delta S_\\text{m}")
"""
Standard change of :symbols:`entropy` per unit :symbols:`amount_of_substance`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature`.
"""

standard_molar_heat_capacity_change = clone_as_symbol(symbols.molar_heat_capacity, display_symbol="Delta(c_m)", display_latex="\\Delta c_\\text{m}")
"""
Standard change of :symbols:`molar_heat_capacity`.
"""

law = Eq(
    standard_molar_gibbs_energy_change, standard_molar_enthalpy_change -
    temperature * standard_molar_entropy_change - standard_molar_heat_capacity_change * temperature *
    (log(temperature / quantities.standard_laboratory_temperature) +
    (quantities.standard_laboratory_temperature / temperature) - 1))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(standard_thermal_effect_=standard_molar_enthalpy_change,
    standard_change_entropy_=standard_molar_entropy_change,
    temperature_=temperature,
    standard_change_heat_capacity_=standard_molar_heat_capacity_change)
@validate_output(standard_molar_gibbs_energy_change)
def calculate_standard_change_isobaric_potential(
        standard_thermal_effect_: Quantity, standard_change_entropy_: Quantity,
        temperature_: Quantity, standard_change_heat_capacity_: Quantity) -> Quantity:
    result_expr = solve(law, standard_molar_gibbs_energy_change,
        dict=True)[0][standard_molar_gibbs_energy_change]
    result_expr = result_expr.subs({
        standard_molar_enthalpy_change: standard_thermal_effect_,
        standard_molar_entropy_change: standard_change_entropy_,
        temperature: temperature_,
        standard_molar_heat_capacity_change: standard_change_heat_capacity_,
    })
    return Quantity(result_expr)
