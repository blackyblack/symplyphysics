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

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Entropy#Entropy_change_formulas_for_simple_processes>`__.

..
    TODO refactor `mass / molar_mass` into `amount_of_substance`
    TODO find other link
"""

from sympy import Eq, solve, log
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)

entropy_change = symbols.entropy
"""
:symbols:`entropy` change during the transition between the two states.
"""

mass = symbols.mass
"""
:symbols:`mass` of the gas.
"""

molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass`, or molecular weight, of the gas.
"""

molar_isochoric_heat_capacity = clone_as_symbol(
    symbols.molar_heat_capacity,
    display_symbol="c_Vm",
    display_latex="c_{V, \\text{m}}",
)
"""
:symbols:`molar_heat_capacity` at constant volume.
"""

final_temperature = clone_as_symbol(symbols.temperature, subscript="1")
"""
:symbols:`temperature` of the final state.
"""

initial_temperature = clone_as_symbol(symbols.temperature, subscript="0")
"""
:symbols:`temperature` of the initial state.
"""

final_volume = clone_as_symbol(symbols.volume, subscript="1")
"""
:symbols:`volume` of the final state.
"""

initial_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
:symbols:`volume` of the initial state.
"""

law = Eq(entropy_change, (mass / molar_mass) *
    ((molar_isochoric_heat_capacity * log(final_temperature / initial_temperature)) +
    (quantities.molar_gas_constant * log(final_volume / initial_volume))))
"""
:laws:symbol::

:laws:latex::
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
