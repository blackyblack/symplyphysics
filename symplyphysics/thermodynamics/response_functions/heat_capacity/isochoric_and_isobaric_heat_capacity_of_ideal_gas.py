"""
Isochoric and isobaric heat capacities of ideal gas
===================================================

The **Mayer's relation** is the relation between heat capacity at constant pressure and that at
constant volume in the case of an ideal gas.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. Gas is ideal.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Mayer%27s_relation>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.thermodynamics.equations_of_state.ideal_gas import ideal_gas_equation
from symplyphysics.thermodynamics.response_functions.heat_capacity import isochoric_and_isobaric_heat_capacities_of_homogeneous_substance as general_mayers_relation
from symplyphysics.thermodynamics.response_functions.thermal_expansion import volumetric_coefficient_of_thermal_expansion as expansion_coefficient_def
from symplyphysics.thermodynamics.response_functions.compressibility import thermodynamic_compressibility as compressibility_def

isobaric_heat_capacity = clone_as_symbol(symbols.heat_capacity, subscript="p")
"""
:symbols:`heat_capacity` of gas at constant :symbols:`pressure`.
"""

isochoric_heat_capacity = clone_as_symbol(symbols.heat_capacity, subscript="V")
"""
:symbols:`heat_capacity` of gas at constant :symbols:`volume`.
"""

amount_of_substance = symbols.amount_of_substance
"""
:symbols:`amount_of_substance`.
"""

law = Eq(
    isobaric_heat_capacity - isochoric_heat_capacity,
    amount_of_substance * quantities.molar_gas_constant,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the Mayer's relation for a homogeneous substance by evaluating the
# thermal expansion coefficient and the isothermal compressibility on the ideal gas
# equation of state.

## Express the volume of the ideal gas as a function of temperature and pressure.
_ideal_gas_volume = solve(ideal_gas_equation.law, ideal_gas_equation.volume)[0]

## Evaluate the volumetric thermal expansion coefficient of the ideal gas. Note that
## the pressure is held constant during the expansion, as required by the definition.
_expansion_coefficient = expansion_coefficient_def.law.rhs.subs(
    expansion_coefficient_def.volume(expansion_coefficient_def.temperature,
    expansion_coefficient_def.parameters),
    _ideal_gas_volume,
).doit()

## Evaluate the isothermal compressibility of the ideal gas. Note that the temperature
## is held constant, so the compressibility is indeed the isothermal one.
_compressibility = compressibility_def.law.rhs.subs(
    compressibility_def.volume(compressibility_def.pressure, compressibility_def.parameters),
    _ideal_gas_volume,
).doit()

## Substitute the ideal gas volume and response functions into the general Mayer's relation.
_heat_capacity_difference = general_mayers_relation.law.rhs.subs({
    general_mayers_relation.volume: _ideal_gas_volume,
    general_mayers_relation.temperature: ideal_gas_equation.temperature,
    general_mayers_relation.thermal_expansion_coefficient: _expansion_coefficient,
    general_mayers_relation.isothermal_compressibility: _compressibility,
})

assert expr_equals(_heat_capacity_difference, law.rhs)


@validate_input(
    isobaric_heat_capacity_=isobaric_heat_capacity,
    amount_of_substance_=amount_of_substance,
)
@validate_output(isochoric_heat_capacity)
def calculate_isochoric_heat_capacity(
    isobaric_heat_capacity_: Quantity,
    amount_of_substance_: Quantity,
) -> Quantity:
    expr = solve(law, isochoric_heat_capacity)[0]
    result = expr.subs({
        isobaric_heat_capacity: isobaric_heat_capacity_,
        amount_of_substance: amount_of_substance_,
    })
    return Quantity(result)
