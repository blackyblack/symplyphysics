from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description

# Law: E = F / q0
## E - electric field
## F - electrostatic force
## q0 - test_charge

electric_field = Symbol("electric_field", units.force / units.charge)
electrostatic_force = clone_symbol(symbols.dynamics.force)
test_charge = Symbol("test_charge", units.charge)

law = Eq(electric_field, electrostatic_force / test_charge)


def print_law() -> str:
    return print_expression(law)


@validate_input(electrostatic_force_=electrostatic_force, test_charge_=test_charge)
@validate_output(electric_field)
def calculate_electric_field(electrostatic_force_: Quantity, test_charge_: Quantity) -> Quantity:
    result = solve(law, electric_field)[0]
    result_field = result.subs({
        electrostatic_force: electrostatic_force_,
        test_charge: test_charge_,
    })
    return Quantity(result_field)
