r"""
Pressure of ideal gas from height and temperature
=================================================

The *barometric formula* determines the dependence of the pressure or density of a gas on the height in the gravity field.

**Notation:**

#. :math:`g` is the acceleration due to gravity.
#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Conditions:**

#. The gas is ideal.
#. The gas is in a uniform gravity field.
"""

from sympy import (Eq, solve, exp)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input, validate_output,
    clone_symbol)

final_pressure = Symbol("final_pressure", units.pressure)
"""
Pressure of the gas at final height.

Symbol:
    :code:`p`
"""

initial_pressure = Symbol("initial_pressure", units.pressure)
"""
Pressure of the gas at initial height.

Symbol:
    :code:`p0`

Latex:
    :math:`p_0`
"""

molecular_mass = symbols.mass
"""
:attr:`~symplyphysics.symbols.mass` of a single gas molecule.
"""

height_change = Symbol("height_change", units.length)
r"""
Change in height between :math:`p_0` and :math:`p`.

Symbol:
    :code:`dh`

Latex:
    :math:`\Delta h`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the gas.
"""

law = Eq(
    final_pressure,
    initial_pressure * exp(-units.acceleration_due_to_gravity * molecular_mass * height_change /
    (units.boltzmann_constant * temperature)))
r"""
:code:`p = p0 * exp(-1 * g * m * dh / (k_B * T))`

Latex:
    .. math::
        p = p_0 \exp \left( - \frac{g m \Delta h}{k_\text{B} T} \right)
"""


@validate_input(initial_pressure_=initial_pressure,
    molecular_mass_=molecular_mass,
    height_change_=height_change,
    temperature_=temperature)
@validate_output(final_pressure)
def calculate_final_pressure(initial_pressure_: Quantity, molecular_mass_: Quantity,
    height_change_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, final_pressure, dict=True)[0][final_pressure]
    result_final_pressure = result_expr.subs({
        initial_pressure: initial_pressure_,
        molecular_mass: molecular_mass_,
        height_change: height_change_,
        temperature: temperature_
    })
    return Quantity(result_final_pressure)
