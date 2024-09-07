"""
Linear charge density from charge and length
============================================

*Linear charge density* is the quantity of charge per unit length, at any point on a line charge distribution.
See :doc:`laws.quantities.quantity_is_linear_density_times_length` for a general version of this law.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, SymbolNew, validate_input, validate_output)

linear_charge_density = SymbolNew("lambda", units.charge / units.length, display_latex="\\lambda")
"""
Linear charge density of the object.
"""

charge = SymbolNew("q", units.charge)
"""
Charge accumulated in the object.
"""

length = SymbolNew("l", units.length)
"""
Length of the object.
"""

definition = Eq(linear_charge_density, charge / length)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_=charge, length_=length)
@validate_output(linear_charge_density)
def calculate_linear_charge_density(charge_: Quantity, length_: Quantity) -> Quantity:
    result_expr = solve(definition, linear_charge_density, dict=True)[0][linear_charge_density]
    result_linear_charge_density = result_expr.subs({
        charge: charge_,
        length: length_,
    })
    return Quantity(result_linear_charge_density)
