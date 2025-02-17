"""
Molar internal energy
=====================

If the equation of state is known, the internal energy of a substance can be found
as a function of volume at constant temperature.

**Conditions:**

#. The fluid is homogeneous and in a single phase state.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Internal_energy_and_specific_heat_at_constant_volume>`__.
"""

from sympy import Eq, Integral
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

molar_internal_energy = Symbol("u_m", units.energy / units.amount_of_substance, display_latex="u_\\text{m}")
"""
:symbols:`internal_energy` of the van der Waals fluid per unit :symbols:`amount_of_substance`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the van der Waals fluid.
"""

isochoric_molar_heat_capacity = clone_as_function(symbols.molar_heat_capacity, [temperature], display_symbol="c_Vm", display_latex="c_{V, \\text{m}}")
"""
:symbols:`molar_heat_capacity` at constant :symbols:`volume` as a function of :attr:`~temperature`.
"""

attractive_forces_parameter = symbols.attractive_forces_parameter
"""
:symbols:`attractive_forces_parameter`.
"""

molar_volume = symbols.molar_volume
"""
:symbols:`molar_volume`.
"""

law = Eq(
    molar_internal_energy,
    Integral(isochoric_molar_heat_capacity(temperature), temperature) -
    attractive_forces_parameter / molar_volume,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    isochoric_molar_heat_capacity_=isochoric_molar_heat_capacity,
    temperature_=temperature,
    bonding_forces_parameter_=attractive_forces_parameter,
    molar_volume_=molar_volume,
)
@validate_output(molar_internal_energy)
def calculate_internal_energy(
    isochoric_molar_heat_capacity_: Quantity,
    temperature_: Quantity,
    bonding_forces_parameter_: Quantity,
    molar_volume_: Quantity,
) -> Quantity:
    # Note that internal energy is only known up to a constant term
    # Isochoric heat capacity is assumed to be a constant independent of temperature

    isochoric_molar_heat_capacity_function = isochoric_molar_heat_capacity_
    result = law.rhs.subs(isochoric_molar_heat_capacity(temperature),
        isochoric_molar_heat_capacity_function).doit().subs({
        temperature: temperature_,
        attractive_forces_parameter: bonding_forces_parameter_,
        molar_volume: molar_volume_,
        })
    return Quantity(result)
