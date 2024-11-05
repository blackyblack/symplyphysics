"""
Internal energy of ideal gas via temperature
============================================

Internal energy of an ideal gas is the sum of the kinetic energy of all of its molecules.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The gas is ideal.

..
    TODO replace `mass/molar_mass` with `amount_of_substance`
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

internal_energy = symbols.internal_energy
"""
:symbols:`internal_energy` of the gas.
"""

mass = symbols.mass
"""
:symbols:`mass` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass` of the gas.
"""

law = Eq(internal_energy, 3 * mass * quantities.molar_gas_constant * temperature / (2 * molar_mass))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_of_gas_=mass, temperature_=temperature, mole_mass_=molar_mass)
@validate_output(internal_energy)
def calculate_inner_energy(mass_of_gas_: Quantity, temperature_: Quantity,
    mole_mass_: Quantity) -> Quantity:
    solved = solve(law, internal_energy, dict=True)[0][internal_energy]
    result_expr = solved.subs({
        mass: mass_of_gas_,
        temperature: temperature_,
        molar_mass: mole_mass_
    })
    return Quantity(result_expr)
