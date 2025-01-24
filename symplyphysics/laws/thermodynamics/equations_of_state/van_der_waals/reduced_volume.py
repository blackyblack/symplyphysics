r"""
Reduced volume
==============

See :ref:`vdw_reduced_units_def`.

**Note:**

#. Specific or molar volumes can be used in the right-hand side of the law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

reduced_volume = SymbolNew("V_r", dimensionless)
"""
Reduced :symbols:`volume` of the van der Waals fluid.
"""

volume = symbols.volume
"""
:symbols:`volume` of the van der Waals fluid.
"""

critical_volume = clone_as_symbol(symbols.volume, subscript="\\text{c}")
"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.critical_molar_volume`
and :symbols:`volume`.
"""

law = Eq(reduced_volume, volume / critical_volume)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    volume_=volume,
    critical_volume_=critical_volume,
)
@validate_output(reduced_volume)
def calculate_reduced_volume(
    volume_: Quantity,
    critical_volume_: Quantity,
) -> float:
    result = law.rhs.subs({
        volume: volume_,
        critical_volume: critical_volume_,
    })
    return convert_to_float(result)
