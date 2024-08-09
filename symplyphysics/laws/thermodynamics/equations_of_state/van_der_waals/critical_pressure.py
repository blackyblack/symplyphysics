"""
Critical pressure
=================

Critical parameters of the van der Waals equation of state are such value of volume, pressure, and
temperature at which the isotherm has an inflection point whose tangent at that point is zero, i.e.
the first and second derivatives of pressure with respect to volume at constant temperature are zero.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

critical_pressure = Symbol("critical_pressure", units.pressure)
r"""
Critical pressure of the van der Waals fluid.

Symbol:
    :code:`p_c`

Latex:
    :math:`p_\text{c}`
"""

attractive_forces_parameter = Symbol(
    "attractive_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
"""
Parameter of the van der Waals equation denoting the magnitude of attractive
forces between gas molecules.

Symbol:
    :code:`a`
"""

excluded_volume_parameter = Symbol(
    "excluded_volume_parameter",
    units.volume / units.amount_of_substance,
)
"""
Parameter of the van der Waals equation denoting an excluded molar volume
due to a finite size of molecules.

Symbol:
    :code:`b`
"""

law = Eq(
    critical_pressure,
    attractive_forces_parameter / (27 * excluded_volume_parameter**2),
)
r"""
:code:`p_c = a / (27 * b^2)`

Latex:
    .. math::
        p_\text{c} = \frac{a}{27 b^2}
"""


@validate_input(
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
)
@validate_output(critical_pressure)
def calculate_critical_pressure(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    pressure_ = law.rhs.subs({
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(pressure_)
