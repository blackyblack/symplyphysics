"""
Isobaric potential of temperature dependent heat capacity
=========================================================

The isobaric potential of a reaction is a value whose change during a chemical reaction
is equal to the change in the internal energy of the system. The isobaric potential
shows how much of the total internal energy of the system can be used for chemical
transformations. Thermal effect of reaction is enthalpy of the system. The heat capacity
coefficients are tabular values for the reaction. They are used to express the
dependence of heat capacity on temperature.

**Notation:**

#. :quantity_notation:`standard_laboratory_temperature`.

**Conditions:**

#. We take into account the dependence of heat capacity on temperature according to the
   Temkin-Schwarzman formula.
#. The process is isobaric-isothermal.

..
    TODO: find link
    TODO: fix file and documentation name
    TODO: find reference to the 'Temkin-Schwarzman formula'
"""

from sympy import (Eq, solve, log)
from symplyphysics import (
    quantities,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

molar_gibbs_energy_change = Symbol("Delta(G_m)", units.energy / units.amount_of_substance, display_latex="\\Delta G_\\text{m}")
"""
:symbols:`gibbs_energy` change, or isobaric potential, per unit :symbols:`amount_of_substance`.
"""

molar_enthalpy_change = Symbol("Delta(H_m)", units.energy / units.amount_of_substance, display_latex="\\Delta H_\\text{m}")
"""
:symbols:`enthalpy` change, or thermal effect, per unit :symbols:`amount_of_substance`.
"""

molar_entropy = Symbol("S_m", units.energy / units.amount_of_substance / units.temperature, display_latex="S_\\text{m}")
"""
:symbols:`entropy` per unit :symbols:`amount_of_substance`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature`.
"""

first_capacity_coefficient = Symbol("a", units.energy / units.amount_of_substance / units.temperature)
"""
First tabular coefficient of heat capacity.
"""

second_capacity_coefficient = Symbol("b", units.energy / units.amount_of_substance / units.temperature**2)
"""
Second tabular coefficient of heat capacity.
"""

third_capacity_coefficient = Symbol("c", units.energy * units.temperature / units.amount_of_substance)
"""
Third tabular coefficient of heat capacity.
"""

law = Eq(
    molar_gibbs_energy_change, molar_enthalpy_change - temperature * molar_entropy - temperature *
    (first_capacity_coefficient * (log(quantities.standard_laboratory_temperature / temperature) +
    (quantities.standard_laboratory_temperature / temperature) - 1) + second_capacity_coefficient *
    ((temperature / 2) + (quantities.standard_laboratory_temperature**2 /
    (2 * temperature)) - quantities.standard_laboratory_temperature) + third_capacity_coefficient *
    ((temperature**(-2) / 2) - (quantities.standard_laboratory_temperature**(-1) / temperature) +
    (quantities.standard_laboratory_temperature**(-2) / 2))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    thermal_effect_=molar_enthalpy_change,
    entropy_=molar_entropy,
    temperature_=temperature,
    coefficient_capacity_1_=first_capacity_coefficient,
    coefficient_capacity_2_=second_capacity_coefficient,
    coefficient_capacity_3_=third_capacity_coefficient,
)
@validate_output(molar_gibbs_energy_change)
def calculate_isobaric_potential(thermal_effect_: Quantity, entropy_: Quantity,
    temperature_: Quantity, coefficient_capacity_1_: Quantity,
    coefficient_capacity_2_: Quantity, coefficient_capacity_3_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    result_expr = solve(law, molar_gibbs_energy_change, dict=True)[0][molar_gibbs_energy_change]
    result_expr = result_expr.subs({
        molar_enthalpy_change: thermal_effect_,
        molar_entropy: entropy_,
        temperature: temperature_,
        first_capacity_coefficient: coefficient_capacity_1_,
        second_capacity_coefficient: coefficient_capacity_2_,
        third_capacity_coefficient: coefficient_capacity_3_,
    })
    return Quantity(result_expr)
