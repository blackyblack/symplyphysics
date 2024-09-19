r"""
Internal energy of ideal gas via temperature
============================================

Internal energy of an ideal gas is the sum of the kinetic energy of all of its molecules.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The gas is ideal.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols, quantities)

internal_energy = Symbol("energy", units.energy)
"""
Internal energy of the gas.

Symbol:
    :code:`U`
"""

mass = symbols.mass
"""
:symbols:`mass` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
"""
Mass of gas per unit amount of substance.

Symbol:
    :code:`M`
"""

law = Eq(internal_energy, 1.5 * mass * quantities.molar_gas_constant * temperature / molar_mass)
r"""
:code:`U = (3 / 2) * (m / M) * R * T`

Latex:
    .. math::
        U = \frac{3}{2} \frac{m}{M} R T
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
