from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The internal energy of an ideal gas consisting of molecules whose energy is all kinetic
## depends only on the degrees of freedom of the molecule and the temperature of the gas.
## From that can be derived the expression of the isochoric heat capacity of ideal gases.

# Law: C_V = (f / 2) * R
## C_V - isochoric molar heat capacity
## f - degrees of freedom of gas molecules
## R - molar gas constant

# Note
## - For applications, see [internal energy law](./internal_energy_of_ideal_gas_is_proportional_to_temperature)
## - f = 3 for monatomic molecules,
##   f = 5 for diatomic molecules,
##   f = 6 for non-linear polyatomic molecules

# Conditions
## - This is the classical theory of heat capacity of gases, for more accurate represention refer to
##   the quantum theory, which accounts for the "freezing" of degrees of freedom and other phenomena.

isochoric_molar_heat_capacity = Symbol(
    "isochoric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance),
)
degrees_of_freedom = Symbol("degrees_of_freedom", integer=True)

law = Eq(isochoric_molar_heat_capacity, (degrees_of_freedom / 2) * units.molar_gas_constant)


@validate_input(degrees_of_freedom_=degrees_of_freedom)
@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity(degrees_of_freedom_: int) -> Quantity:
    result = law.rhs.subs(degrees_of_freedom, degrees_of_freedom_)
    return Quantity(result)
