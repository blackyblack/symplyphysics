"""
Reaction equilibrium constant via standard Gibbs energy
=======================================================

The equilibrium constant is a value that determines for a given chemical reaction the ratio between
thermodynamic activities (or, depending on the conditions of the reaction, partial pressures,
concentrations or fugitives) of reactants and products in a state of chemical equilibrium.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The term "standard" applies to a solution of an infinite dilution and of a hypothetical standard
   concentration, typically 1 mol/kg.
"""

from sympy import Eq, solve, exp
from symplyphysics import (
    symbols,
    units,
    quantities,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

equilibrium_constant = Symbol("equilibrium_constant", dimensionless)
"""
Equilibrium constant of the reaction.

Symbol:
    :code:`K`
"""

reaction_standard_gibbs_energy = Symbol("reaction_standard_gibbs_energy",
    units.energy / units.amount_of_substance)
r"""
Reaction standard Gibbs energy, which is the sum of the standard Gibbs energies of the reaction
products minus that of reactants.

Symbol:
    :code:`Delta(G)`

Latex:
    :math:`\Delta G`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

law = Eq(equilibrium_constant,
    exp(-reaction_standard_gibbs_energy / (quantities.molar_gas_constant * temperature)))
r"""
:code:`K = exp(-1 * Delta(G) / (R * T))`

Latex:
    .. math::
        K = \exp \left( - \frac{\Delta G}{R T} \right)
"""


@validate_input(standard_change_isobaric_potential_=reaction_standard_gibbs_energy,
    temperature_=temperature)
@validate_output(equilibrium_constant)
def calculate_equilibrium_constant(standard_change_isobaric_potential_: Quantity,
    temperature_: Quantity) -> float:
    result_expr = solve(law, equilibrium_constant, dict=True)[0][equilibrium_constant]
    result_expr = result_expr.subs({
        reaction_standard_gibbs_energy: standard_change_isobaric_potential_,
        temperature: temperature_
    })
    return convert_to_float(result_expr)
