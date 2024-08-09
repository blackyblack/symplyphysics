r"""
Critical temperature
====================

Critical parameters of the van der Waals equation of state are such values of volume, pressure, and
temperature at which the isotherm has an inflection point whose tangent at that point is zero, i.e.
the first and second derivatives of pressure with respect to volume at constant temperature are zero.

**Notation:**

#. :math:`R` is the molar gas constant.
"""

from sympy import Eq
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

critical_temperature = clone_symbol(symbols.thermodynamics.temperature, "critical_temperature")
r"""
Critical :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the van der Waals fluid.

Symbol:
    :code:`T_c`

Latex:
    :math:`T_\text{c}`
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
    critical_temperature,
    (8 * attractive_forces_parameter) / (27 * units.molar_gas_constant * excluded_volume_parameter),
)
r"""
:code:`T_c = 8 * a / (27 * R * b)`

Latex:
    .. math::
        T_\text{c} = \frac{8 a}{27 R b}
"""


@validate_input(
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
)
@validate_output(critical_temperature)
def calculate_critical_temperature(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    temperature_ = law.rhs.subs({
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(temperature_)
