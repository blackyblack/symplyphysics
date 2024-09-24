r"""
Change in entropy of ideal gas from volume and temperature
==========================================================

*Entropy* is a property of a thermodynamic system that expresses the direction or outcome of spontaneous changes in the system.
The term was introduced to explain the relationship of the internal energy that is available or unavailable for transformations in form of heat and work.
Entropy predicts that certain processes are irreversible or impossible, despite not violating the conservation of energy.
The definition of entropy is central to the establishment of the second law of thermodynamics, which states that the entropy of isolated systems cannot decrease with time,
as they always tend to arrive at a state of thermodynamic equilibrium, where the entropy is highest.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The gas is ideal.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols,
    clone_as_symbol, quantities)
from symplyphysics.core.functions import log

entropy_change = Symbol("entropy_change", units.energy / units.temperature)
"""
Entropy change during the transition between the two states.

Symbol:
    :code:`S`
"""

mass = symbols.mass
"""
:symbols:`mass` of the gas.
"""

molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
"""
Molar mass, or molecular weight, of the gas.

Symbol:
    :code:`M`
"""

molar_isochoric_heat_capacity = Symbol(
    "molar_isochoric_heat_capacity", units.energy / (units.temperature * units.amount_of_substance))
r"""
Heat capacity at constant volume per unit amount of substance.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

final_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T1",
    display_latex="T_1")
"""
:symbols:`temperature` of the final state.
"""

initial_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T0",
    display_latex="T_0")
"""
:symbols:`temperature` of the initial state.
"""

final_volume = Symbol("final_volume", units.volume)
"""
Volume of the final state.

Symbol:
    :code:`V1`

Latex:
    :math:`V_1`
"""

initial_volume = Symbol("initial_volume", units.volume)
"""
Volume of the initial state.

Symbol:
    :code:`V0`

Latex:
    :math:`V_0`
"""

law = Eq(entropy_change, (mass / molar_mass) *
    ((molar_isochoric_heat_capacity * log(final_temperature / initial_temperature)) +
    (quantities.molar_gas_constant * log(final_volume / initial_volume))))
r"""
:code:`S = (m / M) * (C_V * log(T1 / T0) + R * log(V1 / V0))`

Latex:
    .. math::
        S = \frac{m}{M} \left( C_V \log \frac{T_1}{T_0} + R \log \frac{V_1}{V_0} \right)
"""


@validate_input(mass_=mass,
    molar_mass_=molar_mass,
    molar_heat_capacity_=molar_isochoric_heat_capacity,
    temperature_limits=initial_temperature,
    volume_limits=initial_volume)
@validate_output(entropy_change)
def calculate_entropy_change(mass_: Quantity, molar_mass_: Quantity, molar_heat_capacity_: Quantity,
    temperature_limits: tuple[Quantity, Quantity], volume_limits: tuple[Quantity,
    Quantity]) -> Quantity:
    initial_temperature_, final_temperature_ = temperature_limits
    start_volume_, final_volume_ = volume_limits
    result_expr = solve(law, entropy_change, dict=True)[0][entropy_change]
    result_entropy_change = result_expr.subs({
        mass: mass_,
        molar_mass: molar_mass_,
        molar_isochoric_heat_capacity: molar_heat_capacity_,
        initial_temperature: initial_temperature_,
        final_temperature: final_temperature_,
        initial_volume: start_volume_,
        final_volume: final_volume_,
    })
    return Quantity(result_entropy_change)
