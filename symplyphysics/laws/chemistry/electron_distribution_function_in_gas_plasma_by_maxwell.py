from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float
)
from sympy.physics.units import elementary_charge
from symplyphysics.core.symbols.probability import Probability

# Description
## In a gas discharge, electrons have a wide range of energies, which is described by the electron energy distribution function.
## Electrons in a gas-discharge plasma acquire their energy under the action of an electric field. Energy consumption occurs due to
## elastic and, especially, inelastic collisions with atoms. In addition, energy exchange between electrons is also possible in plasma.
## Depending on the relationship between all these factors, different electron energy distributions are established. Under equilibrium
## conditions, the Maxwell distribution is most common.

## Law is: f = 2.07 * ((q * U)^0.5 / We^1.5) * exp(-1.55 * q * U / We), where
## f - the value of the electron distribution function,
## q - elementary charge,
## U - voltage between electrodes,
## We - electron energy.

value_of_distribution_function = Symbol("value_of_distribution_function", dimensionless)

voltage_between_electrodes = Symbol("voltage_between_electrodes", units.voltage)
electron_energy = Symbol("electron_energy", units.energy)

constant_energy = Quantity(2.07 * units.electronvolt)
expression_1 = constant_energy * (elementary_charge * voltage_between_electrodes)**0.5 / electron_energy**1.5
expression_2 = exp(-1.55 * elementary_charge * voltage_between_electrodes / electron_energy)
law = Eq(value_of_distribution_function, expression_1 * expression_2)


@validate_input(voltage_between_electrodes_=voltage_between_electrodes,
    electron_energy_=electron_energy)
@validate_output(value_of_distribution_function)
def calculate_value_of_distribution_function(voltage_between_electrodes_: Quantity,
    electron_energy_: Quantity) -> Probability:
    result_expr = solve(law, value_of_distribution_function,
        dict=True)[0][value_of_distribution_function]
    result_expr = result_expr.subs({
        voltage_between_electrodes: voltage_between_electrodes_,
        electron_energy: electron_energy_,
    })
    return Probability(convert_to_float(result_expr))
